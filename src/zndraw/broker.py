"""Module-level broker for external TaskIQ workers.

Usage::

    taskiq worker zndraw.broker:broker

The broker connects to Redis (via ``ZNDRAW_SERVER_REDIS_URL``) and registers
all built-in extensions. The executor connects back to the FastAPI
server at ``ZNDRAW_SERVER_INTERNAL_URL``.
"""

from taskiq_redis import ListQueueBroker

from zndraw.config import Settings
from zndraw.database import _collect_extensions
from zndraw.executor import InternalExtensionExecutor
from zndraw.providers.bootstrap import (
    build_internal_providers_resolver,
    register_filebrowser_providers,
)
from zndraw_joblib import register_internal_tasks

settings = Settings()

if settings.redis_url is None:
    raise RuntimeError(
        "ZNDRAW_SERVER_REDIS_URL must be set for external TaskIQ workers. "
        "External workers cannot use fakeredis."
    )

broker = ListQueueBroker(settings.redis_url, queue_name=settings.task_queue_name)

server_url = settings.internal_url or f"http://localhost:{settings.port}"

executor = InternalExtensionExecutor(
    base_url=server_url,
    providers_resolver=build_internal_providers_resolver(settings),
)

register_internal_tasks(broker, _collect_extensions(), executor)

register_filebrowser_providers(broker, base_url=server_url, settings=settings)
