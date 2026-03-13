"""Edit lock context manager for the ZnDraw client."""

from __future__ import annotations

import logging
import threading
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from zndraw.client.api import APIManager

log = logging.getLogger(__name__)


@dataclass
class ZnDrawLock:
    """Context manager for room-level edit lock via REST API.

    Acquires the edit lock on enter, releases on exit.
    Automatically refreshes the lock every EDIT_LOCK_REFRESH seconds.

    Parameters
    ----------
    api : APIManager
        The API manager instance.
    msg : str | None
        Optional message describing the lock purpose.
    """

    api: APIManager
    msg: str | None = None
    _stop: threading.Event = field(default_factory=threading.Event, init=False)
    _refresh_thread: threading.Thread | None = field(default=None, init=False)
    _lock_token: str | None = field(default=None, init=False)

    def __enter__(self) -> ZnDrawLock:
        """Acquire the edit lock."""
        result = self.api.edit_lock_acquire(self.msg)
        self._lock_token = result["lock_token"]
        self._stop.clear()
        self._refresh_thread = threading.Thread(target=self._refresh, daemon=True)
        self._refresh_thread.start()
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> bool:
        """Release the edit lock."""
        self._stop.set()
        if self._refresh_thread is not None:
            self._refresh_thread.join(timeout=2.0)
        try:
            self.api.edit_lock_release(self._lock_token)
        except Exception as e:  # noqa: BLE001
            log.warning("Failed to release edit lock: %s", e)
        self._lock_token = None
        return False

    def _refresh(self) -> None:
        """Periodically refresh the edit lock."""
        from zndraw.schemas import EDIT_LOCK_REFRESH

        while not self._stop.wait(EDIT_LOCK_REFRESH):
            if self._lock_token is None:
                break
            try:
                self.api.edit_lock_refresh(self._lock_token, self.msg)
            except Exception as e:  # noqa: BLE001
                log.warning("Failed to refresh edit lock: %s", e)
                break

    @property
    def lock_token(self) -> str | None:
        """The current lock token, or None if not held."""
        return self._lock_token
