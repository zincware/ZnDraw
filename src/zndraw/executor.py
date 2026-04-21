"""InternalExecutor implementation for zndraw-joblib.

Runs built-in extensions (modifiers, selections, analysis) in a thread
because the ZnDraw client uses synchronous HTTP and Socket.IO.

Status transitions use the same JobManager HTTP API as external workers.
"""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable

log = logging.getLogger(__name__)


@dataclass
class InternalExtensionExecutor:
    """Executes built-in extensions using a per-task JWT token.

    Only the server base URL is captured at creation time.
    The JWT token is passed per-task at call time.

    Parameters
    ----------
    base_url
        ZnDraw server URL.
    providers_resolver
        Optional zero-arg callable returning a dict mapping
        ``<@internal provider full_name> → <resolved handler>``. Injected
        into ``extension.run(vis, providers=...)`` so server-side modifiers
        (e.g. ``LoadFile``) can reach the same backend handle the
        corresponding @internal provider reads through. ``None`` means no
        injection — matches the pre-@internal-providers behavior.
    """

    base_url: str
    providers_resolver: Callable[[], dict[str, Any]] | None = field(
        default=None, repr=False
    )

    async def __call__(
        self,
        extension_cls: type,
        payload: dict[str, Any],
        room_id: str,
        task_id: str,
        token: str,
    ) -> None:
        """Execute extension via asyncio.to_thread."""
        base_url = self.base_url
        providers_resolver = self.providers_resolver

        def _run() -> None:
            from zndraw.client import ZnDraw
            from zndraw_joblib.client import ClaimedTask

            vis = ZnDraw(
                url=base_url,
                room=room_id,
                token=token,
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
                providers = (
                    providers_resolver() if providers_resolver is not None else {}
                )
                instance.run(vis, providers=providers)
            except Exception as e:
                try:
                    if task is not None:
                        vis.jobs.fail(task, str(e))
                    else:
                        vis.jobs.fail_by_id(task_id, str(e))
                except Exception:
                    log.exception("Could not report failure for task %s", task_id)
                raise
            else:
                vis.jobs.complete(task)
            finally:
                vis.disconnect()

        await asyncio.to_thread(_run)
