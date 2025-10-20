import itertools
import logging
from datetime import datetime, timezone
from pathlib import Path

import ase.io
import requests
import znh5md
import traceback


# Tqdm is removed from the generator as it won't render in a Celery worker.
# from tqdm import tqdm
from celery import shared_task
from zndraw.connectivity import add_connectivity
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


@shared_task
def read_file(
    file: str,
    room: str,
    server_url: str = "http://localhost:5000",
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    make_default: bool = False,
    batch_size: int = 10,
    root_path: str | None = None,
) -> None:
    from zndraw import ZnDraw

    file_path = Path(file)
    vis = ZnDraw(room=room, url=server_url, user="uploader", description=f"{file}")
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
        vis.log(f"Reading {slice_info} from {file} (unknown format, attempting with ASE)...")

    # Use vis.lock context manager to lock room during upload
    with vis.lock:
        try:
            frame_iterator = None
            loaded_frame_count = 0  # Track number of frames loaded

            if backend_name == "ZnH5MD":
                # Use ZnH5MD for H5/H5MD files
                if step is not None and step <= 0:
                    vis.log("❌ Error: Step must be a positive integer (e.g., 1, 2, 5)")
                    raise ValueError("Step must be a positive integer for H5MD files.")
                io = znh5md.IO(file_path)

                n_frames = len(io)

                # Validate slice parameters against file length
                if start is not None and start >= n_frames:
                    vis.log(f"❌ Error: Start frame ({start}) exceeds file length ({n_frames} frames)")
                    raise ValueError(f"Start frame {start} exceeds file length")

                if stop is not None and stop > n_frames:
                    vis.log(f"⚠️ Warning: Stop frame ({stop}) exceeds file length ({n_frames}), using end of file")

                s = slice(start, stop, step)
                # The 'indices' method converts the slice into a (start, stop, step)
                # tuple of non-negative integers that can be used with islice.
                _start, _stop, _step = s.indices(n_frames)

                frame_iterator = itertools.islice(io, _start, _stop, _step)
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
                log.info(f"Processing frames from {file_path} in batches of {batch_size}")

                for batch in batch_generator(frame_iterator, batch_size):
                    # compute connectivity for each frame in the batch, if reasonable sized
                    for atoms in batch:
                        if len(atoms) < 1000:
                            add_connectivity(atoms)
                    vis.extend(batch)
                    loaded_frame_count += len(batch)

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
            "file_modified_timestamp": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
        }
        
        # Store metadata via REST API
        requests.post(
            f"{server_url}/api/rooms/{room}/metadata",
            json=metadata
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


@shared_task(bind=True)
def celery_job_worker(self, room: str, server_url: str = "http://localhost:5000"):
    """Celery worker that polls for server-side extension jobs.

    This task continuously polls the /jobs/next endpoint looking for jobs
    with provider="celery". It executes the extension and marks the job
    as completed or failed.

    Args:
        self: Celery task instance (bound)
        room: The room ID to poll for jobs
        server_url: The ZnDraw server URL
    """
    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings
    from zndraw.zndraw import ZnDraw

    worker_id = f"celery:{self.request.id}"

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
        "analysis": analysis,
    }
    log.info(f"Celery worker {worker_id} starting to poll for jobs in room {room}")

    response = requests.post(
        f"{server_url}/api/rooms/{room}/jobs/next",
        json={"workerId": worker_id},
    )

    if response.status_code == 400:
        # No jobs available or worker not idle
        error_msg = response.json().get("error", "")
        if "No jobs available" in error_msg:
            log.debug(f"No jobs available, worker {worker_id}.")
            return
        else:
            log.warning(f"Worker {worker_id} got 400: {error_msg}")
            return

    if response.status_code != 200:
        log.error(f"Failed to fetch next job: {response.status_code} {response.text}")
        return

    job_data = response.json()
    job_id = job_data.get("jobId")
    category = job_data.get("category")
    extension = job_data.get("extension")
    data = job_data.get("data", {})

    log.info(f"Worker {worker_id} picked up job {job_id}: {category}/{extension}")

    try:
        # Validate category and extension
        if category not in category_map:
            raise ValueError(f"Unknown category: {category}")

        if extension not in category_map[category]:
            raise ValueError(
                f"Unknown extension '{extension}' in category '{category}'"
            )

        # Get the extension class and instantiate it with the provided data
        ext_class = category_map[category][extension]
        instance = ext_class(**data)

        # Create a ZnDraw client connected to the specific room
        vis = ZnDraw(room=room, url=server_url, user=worker_id)

        # Run the extension
        instance.run(vis)

        log.info(f"Worker {worker_id} successfully completed job {job_id}")

        # Mark job as completed
        requests.post(
            f"{server_url}/api/rooms/{room}/jobs/{job_id}/complete",
            json={"result": {}, "workerId": worker_id},
        )

    except Exception as e:
        vis.log(f"Error executing job {job_id}: {e} \n{traceback.format_exc()}")
        log.error(
            f"Worker {worker_id} error executing job {job_id}: {e}", exc_info=True
        )
        requests.post(
            f"{server_url}/api/rooms/{room}/jobs/{job_id}/fail",
            json={"error": str(e), "workerId": worker_id},
        )
    finally:
        try:
            vis.disconnect()
        except Exception:
            pass
