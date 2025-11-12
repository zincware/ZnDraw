import itertools
import logging
import traceback
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

    # Use vis.lock context manager to lock room during upload
    # Initialize tracking variables
    loaded_frame_count = 0  # Track number of frames loaded
    max_particles = 0  # Track maximum particle count across all frames

    # with vis.lock(msg="Uploading ..."):
    # TODO: vis._extend() which does extend without aquire/release 
    # so we can do it here?
    try:
        # Create progress tracker with slice info
        progress_description_text = f"Loading {file_path.name}"
        if slice_info != "all frames":
            progress_description_text = f"Loading {file_path.name} ({slice_info})"

        with vis.progress_tracker(progress_description_text) as task_desc:
            frame_iterator = None
            total_expected_frames = None  # Will be set if we know the frame count

            if backend_name == "ZnH5MD":
                # Use ZnH5MD for H5/H5MD files
                if step is not None and step <= 0:
                    vis.log("❌ Error: Step must be a positive integer (e.g., 1, 2, 5)")
                    raise ValueError("Step must be a positive integer for H5MD files.")
                io = znh5md.IO(file_path)

                n_frames = len(io)

                # Validate slice parameters against file length
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
                # The 'indices' method converts the slice into a (start, stop, step)
                # tuple of non-negative integers that can be used with islice.
                _start, _stop, _step = s.indices(n_frames)

                # Calculate expected frames for progress tracking
                total_expected_frames = (_stop - _start + _step - 1) // _step

                frame_iterator = itertools.islice(io, _start, _stop, _step)
            elif backend_name == "ASE-DB":
                # Use ASE database connection for database files
                # Supports: SQLite (.db), JSON (.json), PostgreSQL, MySQL, ASELMDB (.aselmdb)
                # Connection strings for PostgreSQL/MySQL: postgresql://..., mysql://...
                vis.log("Connecting to ASE database...")
                db = ase.db.connect(file_path)
                n_rows = db.count()
                vis.log(f"Database contains {n_rows} structures")

                if n_rows == 0:
                    vis.log("⚠️ Warning: Database is empty")
                    return

                # Validate step
                if step is not None and step <= 0:
                    vis.log("❌ Error: Step must be a positive integer (e.g., 1, 2, 5)")
                    raise ValueError(
                        "Step must be a positive integer for database files."
                    )

                # Validate and adjust slice parameters
                # User provides 0-indexed, database uses 1-indexed IDs internally
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

                # Create list of row IDs to retrieve (convert 0-indexed to 1-indexed)
                # Database IDs are 1-indexed: first row is ID=1
                selected_ids = list(range(_start + 1, _stop + 1, _step))
                vis.log(f"Selecting {len(selected_ids)} structures from database")

                # Track expected frames for progress
                total_expected_frames = len(selected_ids)

                # Create iterator over selected database rows
                def db_iterator():
                    for row_id in selected_ids:
                        try:
                            row = db.get(id=row_id)
                            yield row.toatoms()
                        except KeyError:
                            # Row ID might not exist (gaps in database)
                            vis.log(
                                f"⚠️ Warning: Row ID {row_id} not found in database, skipping"
                            )
                            continue

                frame_iterator = db_iterator()
            else:
                # Use ASE for all other formats (known and unknown)
                # --- This logic is correct for ASE's string-based index ---
                # It properly handles None by creating empty strings, e.g., ":-1:"
                start_str = str(start) if start is not None else ""
                stop_str = str(stop) if stop is not None else ""
                step_str = str(step) if step is not None else ""
                index_str = f"{start_str}:{stop_str}:{step_str}"

                # Use ase.io.iread() with the correctly formatted index string.
                # This may fail for unsupported formats, which will be caught below
                frame_iterator = ase.io.iread(file_path, index=index_str)

            # Now, the batching logic is the same for both file types
            if frame_iterator:
                # A simple log message is often better for background tasks.
                log.info(
                    f"Processing frames from {file_path} in batches of {batch_size}"
                )

                for batch in batch_generator(frame_iterator, batch_size):
                    # Track max particle count
                    for atoms in batch:
                        max_particles = max(max_particles, len(atoms))
                    vis.extend(batch)
                    loaded_frame_count += len(batch)

                    # Update task progress if we know total frames
                    if total_expected_frames and total_expected_frames > 0:
                        progress = (loaded_frame_count / total_expected_frames) * 100
                        # Cap at 99 until complete (will auto-complete on context exit)
                        task_desc.update(progress=min(progress, 99))

    except Exception as e:
        # Log the full exception for better debugging
        log.exception(f"An error occurred while reading file {file_path}")
        vis.log(f"❌ Error reading file {file_path}: {e}")
        raise  # Re-raise to exit the context manager properly

    # Report success with frame count
    if slice_info != "all frames":
        vis.log(f"✓ Successfully loaded {loaded_frame_count} frames ({slice_info})")
    else:
        vis.log(f"✓ Successfully loaded {loaded_frame_count} frames")

    # Apply adaptive resolution based on particle count
    if loaded_frame_count > 0 and max_particles > 0:
        adaptive_resolution = calculate_adaptive_resolution(max_particles)
        # Only update if resolution should be reduced
        if adaptive_resolution < 16:
            try:
                from zndraw.geometries import Sphere

                # Get current particles geometry
                current_particles = vis.geometries["particles"]
                # Create new Sphere with updated resolution (Pydantic models are frozen)
                updated_particles = Sphere(
                    **{
                        **current_particles.model_dump(),
                        "resolution": adaptive_resolution,
                    }
                )
                # Write back to server
                vis.geometries["particles"] = updated_particles
                log.info(
                    f"Adaptive resolution: {max_particles} particles -> resolution {adaptive_resolution}"
                )
                vis.log(
                    f"Adjusted render quality for {max_particles} particles (resolution: {adaptive_resolution})"
                )
            except Exception as e:
                log.warning(f"Failed to apply adaptive resolution: {e}")

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
        log.info(f"Stored metadata for room {room}: {metadata}")
    except Exception as e:
        log.error(f"Failed to store metadata for room {room}: {e}")
        # Don't fail the whole task if metadata storage fails

    # Set as default room if requested
    if make_default:
        try:
            requests.put(
                f"{server_url}/api/rooms/default", json={"roomId": room}
            ).raise_for_status()
            log.info(f"Set room {room} as default")
        except requests.RequestException as e:
            log.error(f"Failed to set default room {room}: {e}")
            vis.log("Failed to set room as default.")

    vis.disconnect()

    # Cleanup temporary file if requested
    if cleanup_after and file_path.exists():
        try:
            import os

            os.remove(file_path)
            log.info(f"Cleaned up temporary file: {file_path}")
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

    log.info(
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
        log.info(f"Celery worker {worker_id} completed job {job_id}")
    except Exception as e:
        log.error(
            f"Celery worker {worker_id} failed to execute job {job_id}: {e}",
            exc_info=True,
        )
