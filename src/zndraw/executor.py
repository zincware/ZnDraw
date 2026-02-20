"""InternalExecutor implementation for zndraw-joblib.

Runs built-in extensions (modifiers, selections, analysis) in a thread
because the ZnDraw client uses synchronous HTTP and Socket.IO.

Status transitions use the same JobManager HTTP API as external workers.
"""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass
from typing import Any

from pydantic import SecretStr
from zndraw_joblib.models import TaskStatus

logger = logging.getLogger(__name__)


@dataclass
class InternalExtensionExecutor:
    """Executes built-in extensions as the internal worker user.

    Deployment config is captured at creation time and never serialized
    through Redis. Only per-task data is passed at call time.

    Parameters
    ----------
    base_url
        ZnDraw server URL.
    worker_email
        Email of the internal worker user.
    worker_password
        Password for the internal worker user.
    """

    base_url: str
    worker_email: str
    worker_password: SecretStr

    async def __call__(
        self,
        extension_cls: type,
        payload: dict[str, Any],
        room_id: str,
        task_id: str,
    ) -> None:
        """Execute extension via asyncio.to_thread."""
        base_url = self.base_url
        worker_email = self.worker_email
        worker_password = self.worker_password

        def _run() -> None:
            from zndraw_joblib.client import ClaimedTask

            from zndraw.client import ZnDraw

            vis = ZnDraw(
                url=base_url,
                room=room_id,
                user=worker_email,
                password=worker_password,
            )
            task: ClaimedTask | None = None
            try:
                instance = extension_cls(**payload)
                task = ClaimedTask(
                    task_id=task_id,
                    job_name="",
                    room_id=room_id,
                    extension=instance,
                )
                vis.jobs.start(task)
                instance.run(vis)
            except Exception as e:
                try:
                    if task is not None:
                        vis.jobs.fail(task, str(e))
                    else:
                        vis.jobs._update_task(task_id, TaskStatus.FAILED, str(e))
                except Exception:
                    logger.exception("Could not report failure for task %s", task_id)
                raise
            else:
                vis.jobs.complete(task)
            finally:
                vis.disconnect()

        await asyncio.to_thread(_run)
