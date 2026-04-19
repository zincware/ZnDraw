"""Shared bootstrap for the default @internal filesystem provider.

Used by both the in-process taskiq worker (``zndraw.database.lifespan``)
and the standalone taskiq worker (``zndraw.broker``). Centralizing the
construction keeps both call sites in lockstep.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from zndraw.providers.executor import InternalProviderExecutor

if TYPE_CHECKING:
    from taskiq import AsyncBroker

    from zndraw.config import Settings
    from zndraw_joblib.registry import InternalProviderRegistry


def register_filebrowser_providers(
    broker: AsyncBroker,
    *,
    base_url: str,
    settings: Settings,
) -> InternalProviderRegistry | None:
    """Register the default @internal filesystem provider on ``broker``.

    Returns ``None`` when ``settings.filebrowser_path is None``
    (feature disabled).

    Parameters
    ----------
    broker
        The TaskIQ broker to register provider tasks on.
    base_url
        ZnDraw server URL (e.g. ``http://127.0.0.1:8000``).
    settings
        Application settings supplying ``filebrowser_path`` and
        ``provider_executor_timeout``.
    """
    if settings.filebrowser_path is None:
        return None

    from zndraw.database import _collect_providers
    from zndraw_joblib.registry import register_internal_providers

    executor = InternalProviderExecutor(
        base_url=base_url,
        filebrowser_path=str(Path(settings.filebrowser_path).resolve()),
        timeout_seconds=settings.provider_executor_timeout,
    )
    return register_internal_providers(broker, _collect_providers(), executor)
