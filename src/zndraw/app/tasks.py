from celery import shared_task
import znh5md
from tqdm import tqdm
from zndraw.utils import update_colors_and_radii, atoms_to_dict
import logging
import requests
import time

log = logging.getLogger(__name__)


@shared_task
def read_file() -> None:
    from zndraw import Client

    client = Client(room="testroom", url="http://localhost:5000")
    client.connect()
    io = znh5md.IO("/Users/fzills/tools/zndraw-communication-testing/structures.h5")
    frames = io[:]
    for atoms in tqdm(frames, desc="Uploading frames"):
        update_colors_and_radii(atoms)
        # atoms = atoms.repeat((4, 4, 4))
        client.append(atoms_to_dict(atoms))


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
    from zndraw.zndraw import ZnDraw
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

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
            raise ValueError(f"Unknown extension '{extension}' in category '{category}'")

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
        log.error(f"Worker {worker_id} error executing job {job_id}: {e}", exc_info=True)
        requests.post(
            f"{server_url}/api/rooms/{room}/jobs/{job_id}/fail",
            json={"error": str(e), "workerId": worker_id},
        )



@shared_task
def start_celery_workers(room: str, num_workers: int = 1, server_url: str = "http://localhost:5000"):
    """Start multiple Celery workers for a room.

    Args:
        room: The room ID to poll for jobs
        num_workers: Number of worker tasks to start
        server_url: The ZnDraw server URL
    """
    for i in range(num_workers):
        celery_job_worker.delay(room, server_url)
        log.info(f"Started Celery worker {i+1}/{num_workers} for room {room}")
