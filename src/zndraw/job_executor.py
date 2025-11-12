"""Shared job execution logic for both Celery and remote workers.

This module provides a unified execution path for extension jobs,
ensuring consistent behavior across different worker types.
"""

import logging
import traceback

import requests

from zndraw.app.job_manager import JobStatus

log = logging.getLogger(__name__)


def _fetch_job_details(server_url: str, job_id: str) -> dict:
    """Fetch job details from server via GET /api/jobs/{job_id}.

    Parameters
    ----------
    server_url : str
        ZnDraw server URL
    job_id : str
        Job ID

    Returns
    -------
    dict
        Job details with keys: jobId, room, category, extension, data, public, status

    Raises
    ------
    requests.HTTPError
        If job not found or request fails
    """
    url = f"{server_url}/api/jobs/{job_id}"
    response = requests.get(url)
    response.raise_for_status()

    job_data = response.json()
    log.debug(f"Fetched job {job_id} details: {job_data['category']}/{job_data['extension']}")

    return job_data


def execute_job_for_worker(
    job_id: str,
    server_url: str,
    worker_id: str,
    run_kwargs: dict | None = None,
    extension_registry: dict[str, dict] | None = None,
) -> None:
    """Execute an extension job following the new REST-based workflow.

    This function implements the complete job execution workflow:
    1. Fetches job details via GET /api/jobs/{job_id}
    2. Creates ZnDraw instance in target room
    3. Gets extension class from registry or built-ins
    4. Updates job status to 'processing' via vis.api.update_job_status()
    5. Executes extension.run()
    6. Updates job status to 'completed' or 'failed'
    7. Cleans up resources

    Parameters
    ----------
    job_id : str
        Unique job identifier
    server_url : str
        ZnDraw server URL
    worker_id : str
        Worker identifier (e.g., 'celery:task_id' or socket sid)
    run_kwargs : dict | None
        Additional keyword arguments to pass to extension.run().
        Merged with run_kwargs from extension registration (passed kwargs take precedence).
    extension_registry : dict[str, dict] | None
        Extension registry from parent ZnDraw (for socket workers).
        Maps "category:name" to _ExtensionStore dict with 'extension' and 'run_kwargs'.
        If None, uses built-in extensions only (for Celery workers).
    """
    from zndraw.api_manager import APIManager
    from zndraw.zndraw import ZnDraw

    vis = None
    room = None

    try:
        # Step 1: Fetch job details via GET /api/jobs/{job_id}
        job_data = _fetch_job_details(server_url, job_id)

        room = job_data["room"]
        category = job_data["category"]
        extension = job_data["extension"]
        data = job_data["data"]
        public = job_data["public"]

        log.info(
            f"Worker {worker_id} fetched job {job_id}: {category}/{extension} in room {room}"
        )

        # Step 2: Create ZnDraw instance in target room
        vis = ZnDraw(room=room, url=server_url, user=worker_id)

        # Step 3: Get extension class from registry or built-ins
        if extension_registry is not None:
            # Socket worker: use parent's registered extensions
            ext_key = f"{category}:{extension}"
            ext_store = extension_registry.get(ext_key)
            if ext_store is None:
                raise ValueError(
                    f"Extension '{extension}' not in registry. "
                    f"Registered: {list(extension_registry.keys())}"
                )
            # Extract extension class and run_kwargs from _ExtensionStore
            extension_class = ext_store["extension"]
            stored_run_kwargs = ext_store.get("run_kwargs") or {}
            # Merge: stored kwargs + passed kwargs (passed takes precedence)
            run_kwargs = {**stored_run_kwargs, **(run_kwargs or {})}
        else:
            # Celery worker: use built-in extensions only
            extension_class = _get_extension_class(category, extension, public)

        # Step 4: Update job status to 'processing' via vis.api.update_job_status()
        vis.api.update_job_status(
            job_id=job_id,
            status=JobStatus.PROCESSING,
            worker_id=worker_id,
        )

        log.info(
            f"Worker {worker_id} starting job {job_id}: {category}/{extension}"
        )

        # Step 5: Execute extension
        try:
            instance = extension_class(**data)
            instance.run(vis, **(run_kwargs or {}))

            log.info(f"Worker {worker_id} completed job {job_id}")

            # Step 6a: Mark as completed
            vis.api.update_job_status(
                job_id=job_id,
                status=JobStatus.COMPLETED,
                worker_id=worker_id,
            )

        except Exception as e:
            # Step 6b: Mark as failed
            error_msg = f"Error executing job {job_id}: {e}\n{traceback.format_exc()}"
            vis.log(error_msg)
            log.error(f"Worker {worker_id} failed job {job_id}: {e}", exc_info=True)

            vis.api.update_job_status(
                job_id=job_id,
                status=JobStatus.FAILED,
                worker_id=worker_id,
                error=str(e),
            )

    except Exception as e:
        # Handle errors during setup (vis creation, extension validation, etc.)
        error_msg = f"Error setting up job {job_id}: {e}\n{traceback.format_exc()}"
        log.error(f"Worker {worker_id} setup error for job {job_id}: {e}", exc_info=True)

        # Try to log error to vis if available
        if vis is not None:
            try:
                vis.log(error_msg)
            except Exception:
                pass

        # Mark job as failed using APIManager directly (since vis might not be fully initialized)
        if room is not None:
            try:
                # Create minimal APIManager for status update
                api = APIManager(url=server_url, room=room)
                api.update_job_status(
                    job_id=job_id,
                    status=JobStatus.FAILED,
                    worker_id=worker_id,
                    error=str(e),
                )
            except Exception as status_error:
                log.error(
                    f"Failed to update job status after setup error: {status_error}"
                )

    finally:
        # Step 7: Clean up resources
        if vis is not None:
            try:
                vis.disconnect()
            except Exception as disconnect_error:
                log.warning(
                    f"Error disconnecting vis for job {job_id}: {disconnect_error}"
                )


def _get_extension_class(category: str, extension: str, public: bool):
    """Get extension class from category and extension name.

    Parameters
    ----------
    category : str
        Extension category
    extension : str
        Extension name
    public : bool
        Whether to look in public (global) extensions

    Returns
    -------
    type
        Extension class

    Raises
    ------
    ValueError
        If category or extension not found
    """
    from zndraw.extensions.analysis import analysis
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
        "analysis": analysis,
    }

    # Validate category
    if category not in category_map:
        raise ValueError(
            f"Unknown category: {category}. Valid categories: {list(category_map.keys())}"
        )

    # Validate extension exists
    if extension not in category_map[category]:
        available = list(category_map[category].keys())
        raise ValueError(
            f"Unknown extension '{extension}' in category '{category}'. "
            f"Available extensions: {available}"
        )

    return category_map[category][extension]
