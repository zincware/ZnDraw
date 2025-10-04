import logging
from pathlib import Path

import ase.io
import requests
import znh5md
from celery import shared_task

log = logging.getLogger(__name__)


@shared_task
def read_file(
    file: str,
    room: str,
    server_url: str = "http://localhost:5000",
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    make_default: bool = False,
) -> None:
    from zndraw import ZnDraw

    file_path = Path(file)
    vis = ZnDraw(room=room, url=server_url, user="uploader")
    if not file_path.exists():
        vis.log(f"File {file} does not exist.")
        return
    vis.log(f"Reading file {file}...")
    if file_path.suffix in [".h5", ".h5md"]:
        io = znh5md.IO(file)
        if step is not None and step < 0:
            frames = io[:][start:stop:step]
        else:
            frames = io[start:stop:step]
        vis.extend(frames)
    else:
        try:
            frames = ase.io.read(file, index=slice(start, stop, step))
            if not isinstance(frames, list):
                frames = [frames]
            vis.extend(frames)
        except Exception as e:
            vis.log(f"Error reading file {file}: {e}")
            return

    vis.log(f"Finished reading file {file}.")
    # promote to template
    requests.post(
        f"{server_url}/api/rooms/{room}/promote",
        json={"name": file, "description": f"Data uploaded from file {file}"},
    ).raise_for_status()
    if make_default:
        requests.put(
            f"{server_url}/api/templates/default", json={"template_id": room}
        ).raise_for_status()

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
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings
    from zndraw.zndraw import ZnDraw

    worker_id = f"celery:{self.request.id}"

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
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


@shared_task
def start_celery_workers(
    room: str, num_workers: int = 1, server_url: str = "http://localhost:5000"
):
    """Start multiple Celery workers for a room.

    Args:
        room: The room ID to poll for jobs
        num_workers: Number of worker tasks to start
        server_url: The ZnDraw server URL
    """
    for i in range(num_workers):
        celery_job_worker.delay(room, server_url)
        log.info(f"Started Celery worker {i + 1}/{num_workers} for room {room}")
