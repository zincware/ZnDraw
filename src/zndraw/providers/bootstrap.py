"""Shared bootstrap for the default @internal filesystem provider.

Used by both the in-process taskiq worker (``zndraw.database.lifespan``)
and the standalone taskiq worker (``zndraw.broker``). Centralizing the
construction keeps both call sites in lockstep.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from zndraw.database import _collect_providers
from zndraw.providers.executor import (
    InternalProviderExecutor,
    resolve_internal_provider_handler,
)
from zndraw_joblib.registry import register_internal_providers

if TYPE_CHECKING:
    from collections.abc import Callable

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

    Returns ``None`` when ``settings.filebrowser_enabled`` is False
    (feature disabled).

    Parameters
    ----------
    broker
        The TaskIQ broker to register provider tasks on.
    base_url
        ZnDraw server URL (e.g. ``http://127.0.0.1:8000``).
    settings
        Application settings supplying ``filebrowser_enabled``,
        ``filebrowser_path`` and ``provider_executor_timeout``.
    """
    if not settings.filebrowser_enabled:
        return None

    executor = InternalProviderExecutor(
        base_url=base_url,
        filebrowser_path=str(Path(settings.filebrowser_path).resolve()),
        timeout_seconds=settings.provider_executor_timeout,
    )
    return register_internal_providers(broker, _collect_providers(), executor)


def build_internal_providers_resolver(
    settings: Settings,
) -> Callable[[], dict[str, Any]] | None:
    """Return a zero-arg callable that builds the ``providers`` kwarg for
    @internal extension runs.

    The callable maps each bundled @internal provider's full name to the
    backend handle it reads through — the SAME handle the provider itself
    uses, obtained via the shared ``resolve_internal_provider_handler``.
    Returns ``None`` when ``settings.filebrowser_enabled`` is False.

    Used by ``InternalExtensionExecutor.providers_resolver`` so server-side
    modifiers (e.g. ``LoadFile``) can reach the same filesystem handle.
    """
    if not settings.filebrowser_enabled:
        return None

    filebrowser_path = str(Path(settings.filebrowser_path).resolve())
    provider_classes = _collect_providers()

    def _resolve() -> dict[str, Any]:
        out: dict[str, Any] = {}
        for prov_cls in provider_classes:
            full_name = f"@internal:{prov_cls.category}:{prov_cls.__name__}"
            try:
                out[full_name] = resolve_internal_provider_handler(
                    prov_cls, filebrowser_path=filebrowser_path
                )
            except ValueError:
                # Provider category without a server-side handler — skip.
                continue
        return out

    return _resolve
