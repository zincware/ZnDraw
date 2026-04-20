"""Taskiq-side executor for @internal providers.

Resolves the filesystem handler from configured settings, invokes
``Provider.read(handler)``, and POSTs the result to the server via
the same ``/v1/joblib/providers/{id}/results`` endpoint that remote
providers use.
"""

from __future__ import annotations

import asyncio
import json
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import fsspec
import httpx
from fsspec.implementations.dirfs import DirFileSystem

from zndraw_joblib.exceptions import ProviderExecutionFailed

if TYPE_CHECKING:
    from zndraw_joblib.provider import Provider

log = logging.getLogger(__name__)


@dataclass
class InternalProviderExecutor:
    """Execute an @internal Provider read and POST the result.

    Parameters
    ----------
    base_url
        ZnDraw server URL (e.g. ``http://127.0.0.1:8000``).
    filebrowser_path
        Absolute or relative path rooting the ``filesystem`` handler.
        Relative paths are resolved against the taskiq-worker's cwd
        when the executor is constructed.
    _transport
        Optional httpx transport for testing. Not set in production.
    """

    base_url: str
    filebrowser_path: str
    timeout_seconds: float = 30.0
    _transport: Any = field(default=None, repr=False)

    async def __call__(
        self,
        provider_cls: type[Provider],
        params_json: str,
        provider_id: str,
        request_id: str,
        token: str,
    ) -> None:
        """Execute the provider and POST the result. Raises on upload failure."""
        base_url = self.base_url
        transport = self._transport

        def _run() -> None:
            json_body: dict[str, Any] | None = None
            content: bytes | None = None

            try:
                handler = resolve_internal_provider_handler(
                    provider_cls, filebrowser_path=self.filebrowser_path
                )
                params = json.loads(params_json) if params_json else {}
                instance = provider_cls(**params)
                result = instance.read(handler)
                if provider_cls.content_type == "application/json":
                    content = json.dumps(result).encode()
                else:
                    content = result  # type: ignore[assignment]
                headers = {
                    "Authorization": f"Bearer {token}",
                    "X-Request-Hash": request_id,
                }
            except Exception as err:
                log.exception(
                    "InternalProviderExecutor failed for %s",
                    provider_cls.__name__,
                )
                problem = ProviderExecutionFailed.create(
                    detail=f"{type(err).__name__}: {err}",
                )
                json_body = problem.model_dump(exclude_none=True)
                headers = {
                    "Authorization": f"Bearer {token}",
                    "X-Request-Hash": request_id,
                    "X-Result-Status": "error",
                    "Content-Type": "application/problem+json",
                }

            client_kwargs: dict[str, Any] = {"timeout": self.timeout_seconds}
            if transport is not None:
                client_kwargs["transport"] = transport

            with httpx.Client(**client_kwargs) as client:
                if json_body is not None:
                    resp = client.post(
                        f"{base_url}/v1/joblib/providers/{provider_id}/results",
                        json=json_body,
                        headers=headers,
                    )
                else:
                    resp = client.post(
                        f"{base_url}/v1/joblib/providers/{provider_id}/results",
                        content=content,
                        headers=headers,
                    )
                resp.raise_for_status()

        await asyncio.to_thread(_run)


def resolve_internal_provider_handler(
    provider_cls: type[Provider], *, filebrowser_path: str
) -> Any:
    """Return the fsspec / backend handle a built-in @internal provider reads through.

    Shared by:
    - ``InternalProviderExecutor._run`` (serving the provider's ``read()``).
    - ``InternalExtensionExecutor`` (injecting handlers into extension.run(...,
      providers=) so modifiers like ``LoadFile`` can reach the same handle).

    Raises ``ValueError`` for categories with no server-side handler.
    """
    category = provider_cls.category
    if category == "filesystem":
        return DirFileSystem(path=filebrowser_path, fs=fsspec.filesystem("file"))
    raise ValueError(
        f"No internal handler configured for provider category '{category}'"
    )
