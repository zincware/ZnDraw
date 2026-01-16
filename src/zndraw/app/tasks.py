import itertools
import logging
from datetime import datetime, timezone
from pathlib import Path

import ase.db
import ase.io
import requests
import znh5md

# Tqdm is removed from the generator as it won't render in a Celery worker.
# from tqdm import tqdm
from celery import shared_task

from zndraw.app.file_browser import FORMAT_BACKENDS

log = logging.getLogger(__name__)


def batch_generator(iterable, size):
    """Yields successive n-sized chunks from any iterable."""
    it = iter(iterable)
    while True:
        chunk = list(itertools.islice(it, size))
        if not chunk:
            return
        yield chunk


def calculate_adaptive_resolution(num_particles: int) -> int:
    """Calculate adaptive sphere resolution based on particle count.

    Reduces resolution for large systems to maintain performance.

    Parameters
    ----------
    num_particles : int
        Number of particles in the system

    Returns
    -------
    int
        Resolution value between 6 and 16

    Examples
    --------
    >>> calculate_adaptive_resolution(100)
    16
    >>> calculate_adaptive_resolution(1000)
    14
    >>> calculate_adaptive_resolution(10000)
    10
    >>> calculate_adaptive_resolution(100000)
    8
    """
    if num_particles < 1000:
        return 16
    elif num_particles < 5000:
        return 14
    elif num_particles < 10000:
        return 12
    elif num_particles < 25000:
        return 10
    elif num_particles <= 100000:
        return 8
    else:
        return 6


def upload_frames_to_room(
    vis,
    frame_iterator,
    batch_size: int = 10,
    total_expected_frames: int | None = None,
    progress_description: str = "Loading frames",
) -> tuple[int, int]:
    """Upload frames to a ZnDraw room with progress tracking.

    Parameters
    ----------
    vis : ZnDraw
        Target ZnDraw instance to upload frames to
    frame_iterator : Iterator[ase.Atoms]
        Iterator yielding ASE Atoms objects
    batch_size : int
        Number of frames per batch
    total_expected_frames : int | None
        Expected total frames for progress calculation
    progress_description : str
        Description shown in progress tracker

    Returns
    -------
    tuple[int, int]
        (loaded_frame_count, max_particles)
    """
    loaded_frame_count = 0
    max_particles = 0

    with vis.progress_tracker(progress_description) as task:
        with vis.get_lock(msg="Uploading frames"):
            for batch in batch_generator(frame_iterator, batch_size):
                for atoms in batch:
                    max_particles = max(max_particles, len(atoms))
                vis._extend(batch)
                loaded_frame_count += len(batch)

                if total_expected_frames and total_expected_frames > 0:
                    progress = (loaded_frame_count / total_expected_frames) * 100
                    task.update(progress=min(progress, 99))

    return loaded_frame_count, max_particles


def apply_adaptive_resolution(vis, max_particles: int) -> None:
    """Apply adaptive resolution based on particle count.

    Reduces sphere resolution for large systems to maintain performance.

    Parameters
    ----------
    vis : ZnDraw
        Target ZnDraw instance
    max_particles : int
        Maximum particle count across all frames
    """
    if max_particles <= 0:
        return

    adaptive_resolution = calculate_adaptive_resolution(max_particles)
    if adaptive_resolution >= 16:
        return

    try:
        from zndraw.geometries import Sphere

        current_particles = vis.geometries["particles"]
        updated_particles = Sphere(
            **{**current_particles.model_dump(), "resolution": adaptive_resolution}
        )
        vis.geometries["particles"] = updated_particles
        log.info(
            f"Adaptive resolution: {max_particles} particles -> resolution {adaptive_resolution}"
        )
        vis.log(
            f"Adjusted render quality for {max_particles} particles (resolution: {adaptive_resolution})"
        )
    except Exception as e:
        log.warning(f"Failed to apply adaptive resolution: {e}")


@shared_task
def read_file(
    file: str,
    room: str,
    server_url: str,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    make_default: bool = False,
    batch_size: int = 10,
    root_path: str | None = None,
    cleanup_after: bool = False,
    description: str | None = None,
) -> None:
    from zndraw import ZnDraw
    from zndraw.server_manager import wait_for_server_ready

    # Wait for server to be ready before connecting
    # This handles race conditions during server startup when Celery task
    # is dispatched before the HTTP server is fully listening
    if not wait_for_server_ready(server_url, timeout=30.0):
        log.error(f"Server at {server_url} not ready after 30 seconds")
        raise RuntimeError(f"Server at {server_url} not ready after 30 seconds")

    file_path = Path(file)
    # Use custom description if provided, otherwise default to file path
    room_description = description if description is not None else f"{file}"
    vis = ZnDraw(
        room=room, url=server_url, user="uploader", description=room_description
    )
    if not file_path.exists():
        vis.log(f"File {file} does not exist.")
        return

    # Determine backend based on file extension
    ext = file_path.suffix.lstrip(".").lower()
    backends = FORMAT_BACKENDS.get(ext, ["ASE"])  # Default to ASE for unknown formats
    backend_name = backends[0]

    # Format slice information for user-friendly display
    def format_slice_info(start, stop, step):
        """Format slice parameters for user-friendly display."""
        if start is None and stop is None and step is None:
            return "all frames"

        parts = []
        if start is not None:
            parts.append(f"from frame {start}")
        if stop is not None:
            parts.append(f"to frame {stop - 1}")
        if step is not None and step != 1:
            parts.append(f"every {step} frames")

        return " ".join(parts) if parts else "all frames"

    slice_info = format_slice_info(start, stop, step)

    if ext in FORMAT_BACKENDS:
        vis.log(f"Reading {slice_info} from {file} using {backend_name}...")
    else:
        vis.log(
            f"Reading {slice_info} from {file} (unknown format, attempting with ASE)..."
        )

    # Build progress description
    progress_description = f"Loading {file_path.name}"
    if slice_info != "all frames":
        progress_description = f"Loading {file_path.name} ({slice_info})"

    try:
        frame_iterator = None
        total_expected_frames = None

        if backend_name == "ZnH5MD":
            if step is not None and step <= 0:
                vis.log("❌ Error: Step must be a positive integer (e.g., 1, 2, 5)")
                raise ValueError("Step must be a positive integer for H5MD files.")
            io = znh5md.IO(file_path)
            n_frames = len(io)

            if start is not None and start >= n_frames:
                vis.log(
                    f"❌ Error: Start frame ({start}) exceeds file length ({n_frames} frames)"
                )
                raise ValueError(f"Start frame {start} exceeds file length")

            if stop is not None and stop > n_frames:
                vis.log(
                    f"⚠️ Warning: Stop frame ({stop}) exceeds file length ({n_frames}), using end of file"
                )

            s = slice(start, stop, step)
            _start, _stop, _step = s.indices(n_frames)
            total_expected_frames = (_stop - _start + _step - 1) // _step
            frame_iterator = itertools.islice(io, _start, _stop, _step)

        elif backend_name == "ASE-DB":
            vis.log("Connecting to ASE database...")
            db = ase.db.connect(file_path)
            n_rows = db.count()
            vis.log(f"Database contains {n_rows} structures")

            if n_rows == 0:
                vis.log("⚠️ Warning: Database is empty")
                return

            if step is not None and step <= 0:
                vis.log("❌ Error: Step must be a positive integer (e.g., 1, 2, 5)")
                raise ValueError("Step must be a positive integer for database files.")

            _start = start if start is not None else 0
            _stop = stop if stop is not None else n_rows
            _step = step if step is not None else 1

            if _start >= n_rows:
                vis.log(
                    f"❌ Error: Start row ({_start}) exceeds database size ({n_rows} rows)"
                )
                raise ValueError(f"Start row {_start} exceeds database size")

            if _stop > n_rows:
                vis.log(
                    f"⚠️ Warning: Stop row ({_stop}) exceeds database size ({n_rows}), using end of database"
                )
                _stop = n_rows

            selected_ids = list(range(_start + 1, _stop + 1, _step))
            vis.log(f"Selecting {len(selected_ids)} structures from database")
            total_expected_frames = len(selected_ids)

            def db_iterator():
                for row_id in selected_ids:
                    try:
                        row = db.get(id=row_id)
                        yield row.toatoms()
                    except KeyError:
                        vis.log(
                            f"⚠️ Warning: Row ID {row_id} not found in database, skipping"
                        )
                        continue

            frame_iterator = db_iterator()

        else:
            # Use ASE for all other formats
            start_str = str(start) if start is not None else ""
            stop_str = str(stop) if stop is not None else ""
            step_str = str(step) if step is not None else ""
            index_str = f"{start_str}:{stop_str}:{step_str}"
            frame_iterator = ase.io.iread(file_path, index=index_str)

        # Upload frames using shared helper
        if frame_iterator:
            log.info(f"Processing frames from {file_path} in batches of {batch_size}")
            loaded_frame_count, max_particles = upload_frames_to_room(
                vis,
                frame_iterator,
                batch_size=batch_size,
                total_expected_frames=total_expected_frames,
                progress_description=progress_description,
            )
        else:
            loaded_frame_count, max_particles = 0, 0

    except Exception as e:
        log.exception(f"An error occurred while reading file {file_path}")
        vis.log(f"❌ Error reading file {file_path}: {e}")
        raise

    # Report success
    if slice_info != "all frames":
        vis.log(f"✓ Successfully loaded {loaded_frame_count} frames ({slice_info})")
    else:
        vis.log(f"✓ Successfully loaded {loaded_frame_count} frames")

    # Apply adaptive resolution
    if loaded_frame_count > 0:
        apply_adaptive_resolution(vis, max_particles)

    # Store file metadata
    try:
        stat = file_path.stat()
        # Get relative path from root_path if provided
        if root_path:
            try:
                relative_path = file_path.relative_to(Path(root_path).resolve())
            except ValueError:
                # If not relative to root_path, use absolute path
                relative_path = file_path.absolute()
        else:
            # Fallback: try relative to cwd
            try:
                relative_path = file_path.relative_to(Path.cwd())
            except ValueError:
                # If not relative to cwd, use absolute path
                relative_path = file_path.absolute()

        metadata = {
            "relative_file_path": str(relative_path),
            "file_size": str(stat.st_size),
            "file_modified_timestamp": datetime.fromtimestamp(
                stat.st_mtime, tz=timezone.utc
            ).isoformat(),
        }

        # Store metadata via REST API
        requests.post(
            f"{server_url}/api/rooms/{room}/metadata", json=metadata
        ).raise_for_status()
        log.debug(f"Stored metadata for room {room}: {metadata}")
    except Exception as e:
        log.error(f"Failed to store metadata for room {room}: {e}")
        # Don't fail the whole task if metadata storage fails

    # Set as default room if requested
    if make_default:
        try:
            # Use authentication headers from the ZnDraw API manager
            headers = vis.api._get_headers()
            requests.put(
                f"{server_url}/api/rooms/default",
                json={"roomId": room},
                headers=headers,
            ).raise_for_status()
            log.debug(f"Set room {room} as default")
        except requests.RequestException as e:
            log.error(f"Failed to set default room {room}: {e}")
            vis.log("Failed to set room as default.")

    # Lock template room if configured
    from zndraw.config import get_config

    config = get_config()
    if config.lock_template_room:
        try:
            # Use authentication headers from the ZnDraw API manager
            headers = vis.api._get_headers()
            requests.patch(
                f"{server_url}/api/rooms/{room}",
                json={"locked": True},
                headers=headers,
            ).raise_for_status()
            log.debug(f"Locked template room {room} (admin-locked)")
            vis.log("✓ Room locked by admin (template room)")
        except requests.RequestException as e:
            log.error(f"Failed to lock template room {room}: {e}")
            vis.log("Failed to lock template room.")

    vis.disconnect()

    # Cleanup temporary file if requested
    if cleanup_after and file_path.exists():
        try:
            import os

            os.remove(file_path)
            log.debug(f"Cleaned up temporary file: {file_path}")
        except Exception as e:
            log.error(f"Failed to cleanup temp file {file_path}: {e}")


@shared_task(bind=True)
def celery_job_worker(self, job_data: dict, server_url: str):
    """Celery worker that executes a server-side extension job.

    This task is triggered when a job is submitted for a Celery extension.
    The job has already been created and assigned to this worker - this task
    just needs to execute it following the standard workflow:
    1. Fetch job details via GET /api/jobs/{job_id}
    2. Update status to 'processing'
    3. Execute extension
    4. Update status to 'completed' or 'failed'

    This uses the same execution path as remote workers (job_executor.execute_job_for_worker).

    Args:
        self: Celery task instance (bound)
        job_data: Job data dict with jobId, room, category, extension, etc.
        server_url: The ZnDraw server URL
    """
    from zndraw.job_executor import execute_job_for_worker

    worker_id = f"celery:{self.request.id}"
    job_id = job_data.get("id")

    log.debug(
        f"Celery worker {worker_id} starting job {job_id}: "
        f"{job_data.get('category')}/{job_data.get('extension')} in room {job_data.get('room')}"
    )

    # Use shared job executor (same code as remote workers)
    try:
        execute_job_for_worker(
            job_id=job_id,
            server_url=server_url,
            worker_id=worker_id,
        )
        log.debug(f"Celery worker {worker_id} completed job {job_id}")
    except Exception as e:
        log.error(
            f"Celery worker {worker_id} failed to execute job {job_id}: {e}",
            exc_info=True,
        )
