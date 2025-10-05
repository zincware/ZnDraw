import itertools
import logging
from pathlib import Path

import ase.io
import requests
import znh5md

# Tqdm is removed from the generator as it won't render in a Celery worker.
# from tqdm import tqdm
from celery import shared_task
from zndraw.connectivity import add_connectivity

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
) -> None:
    from zndraw import ZnDraw

    file_path = Path(file)
    vis = ZnDraw(room=room, url=server_url, user="uploader")
    if not file_path.exists():
        vis.log(f"File {file} does not exist.")
        return
    vis.log(f"Reading file {file}...")
    try:
        frame_iterator = None
        if file_path.suffix in [".h5", ".h5md"]:
            if step is not None and step <= 0:
                vis.log("Step must be a positive integer for H5MD files.")
                raise ValueError("Step must be a positive integer for H5MD files.")
            io = znh5md.IO(file_path)

            n_frames = len(io)
            s = slice(start, stop, step)
            # The 'indices' method converts the slice into a (start, stop, step)
            # tuple of non-negative integers that can be used with islice.
            _start, _stop, _step = s.indices(n_frames)

            frame_iterator = itertools.islice(io, _start, _stop, _step)
        else:
            # --- This logic is correct for ASE's string-based index ---
            # It properly handles None by creating empty strings, e.g., ":-1:"
            start_str = str(start) if start is not None else ""
            stop_str = str(stop) if stop is not None else ""
            step_str = str(step) if step is not None else ""
            index_str = f"{start_str}:{stop_str}:{step_str}"

            # Use ase.io.iread() with the correctly formatted index string.
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

    except Exception as e:
        # Log the full exception for better debugging
        log.exception(f"An error occurred while reading file {file_path}")
        vis.log(f"Error reading file {file_path}: {e}")
        return

    vis.log(f"Finished reading file {file}.")
    # promote to template
    try:
        requests.post(
            f"{server_url}/api/rooms/{room}/promote",
            json={"name": file, "description": f"Data uploaded from file {file}"},
        ).raise_for_status()
        if make_default:
            requests.put(
                f"{server_url}/api/templates/default", json={"template_id": room}
            ).raise_for_status()
    except requests.RequestException as e:
        log.error(f"Failed to promote template for room {room}: {e}")
        vis.log("Failed to promote data to template.")

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
        vis.log(f"Error executing job {job_id}: {e}")
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
