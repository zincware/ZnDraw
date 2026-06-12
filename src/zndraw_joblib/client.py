# src/zndraw_joblib/client.py
"""Client SDK for ZnDraw JobLib workers."""

import json
import logging
import random
import signal
import threading
import time
import traceback
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any, ClassVar, Generic, Protocol, TypeVar
from uuid import UUID

import httpx
from pydantic import BaseModel
from zndraw_socketio import SyncClientWrapper

from zndraw_joblib.events import (
    JoinJobRoom,
    JoinProviderRoom,
    LeaveJobRoom,
    LeaveProviderRoom,
    ProviderRequest,
    TaskAvailable,
)
from zndraw_joblib.models import TaskStatus
from zndraw_joblib.provider import Provider
from zndraw_joblib.schemas import (
    JobRegisterRequest,
    ProviderRegisterRequest,
    TaskClaimResponse,
    TaskSubmitRequest,
)

logger = logging.getLogger(__name__)


class Category(StrEnum):
    """Extension category types."""

    MODIFIER = "modifiers"
    SELECTION = "selections"
    ANALYSIS = "analysis"


class Extension(BaseModel, ABC):
    """Base class for all ZnDraw extensions.

    Extensions are Pydantic models that define their parameters as fields.
    The JSON schema is generated from the model for the frontend forms.
    Subclasses must implement run().
    """

    category: ClassVar[Category]

    @abstractmethod
    def run(self, vis: Any, **kwargs: Any) -> Any:
        """Execute the extension logic. Must be overridden in subclasses."""


E = TypeVar("E", bound=Extension)


class ApiManager(Protocol):
    """Protocol for API client that JobManager uses."""

    http: httpx.Client

    def get_headers(self) -> dict[str, str]: ...

    @property
    def base_url(self) -> str: ...

    def raise_for_status(self, response: httpx.Response) -> None:
        """Raise on HTTP error responses.

        Implementations MUST map status codes to these exception types:

        Raises
        ------
        KeyError
            404 responses (resource not found).
        PermissionError
            401/403 responses (authentication/authorization failure).
        ValueError
            409/422 responses (conflict/validation error).
        """
        ...


class ClaimedTask(Generic[E]):
    """A claimed task with its Extension instance and metadata."""

    def __init__(
        self,
        task_id: str,
        job_name: str,
        room_id: str,
        extension: E,
        run_kwargs: dict[str, Any] | None = None,
    ):
        self.task_id = task_id
        self.job_name = job_name
        self.room_id = room_id
        self.extension = extension
        self.run_kwargs: dict[str, Any] = run_kwargs or {}


@dataclass
class _RegisteredProvider:
    """Internal record of a registered provider."""

    id: UUID
    cls: type[Provider]
    handler: Any
    room_id: str


@dataclass
class _OutageState:
    """Tracks server outage for coordinated retry logging across loops."""

    max_unreachable: float
    clock: Callable[[], float] = field(default=time.monotonic)
    _outage_start: float | None = field(default=None, init=False)
    _last_log_time: float = field(default=float("-inf"), init=False)
    _lock: threading.Lock = field(default_factory=threading.Lock, init=False)

    def record_failure(self) -> None:
        """Mark a connection failure; starts the outage clock on first call."""
        with self._lock:
            if self._outage_start is None:
                self._outage_start = self.clock()

    def record_success(self) -> None:
        """Server responded — reset outage state."""
        with self._lock:
            self._outage_start = None

    def elapsed(self) -> float:
        """Seconds since first failure, or 0.0 if not in outage."""
        with self._lock:
            if self._outage_start is None:
                return 0.0
            return self.clock() - self._outage_start

    def should_shutdown(self) -> bool:
        """Return True if outage has exceeded max_unreachable."""
        return self.elapsed() > self.max_unreachable

    def should_log(self, min_interval: float = 5.0) -> bool:
        """Rate-limit log output; returns True at most once per min_interval."""
        with self._lock:
            now = self.clock()
            if now - self._last_log_time >= min_interval:
                self._last_log_time = now
                return True
            return False


class JobManager:
    """Main entry point for workers. Registers jobs and claims tasks.

    When ``execute`` is provided, background threads automatically claim
    and execute tasks after the first ``register()`` or
    ``register_provider()`` call.  Without ``execute``, manual
    ``claim()`` / ``listen()`` still work as before.
    """

    def __init__(
        self,
        api: ApiManager,
        tsio: SyncClientWrapper | None = None,
        *,
        execute: Callable[[ClaimedTask], None] | None = None,
        heartbeat_interval: float = 30.0,
        polling_interval: float = 2.0,
        max_unreachable_seconds: float = 120.0,
        max_startup_retries: int = 5,
    ):
        self.api = api
        self.tsio = tsio
        self._execute = execute
        self._heartbeat_interval = heartbeat_interval
        self._polling_interval = polling_interval
        self._max_unreachable_seconds = max_unreachable_seconds
        self._max_startup_retries = max_startup_retries
        self._registry: dict[str, tuple[type[Extension], dict[str, Any]]] = {}
        self._worker_id: UUID | None = None
        self._providers: dict[str, _RegisteredProvider] = {}

        # Threading state
        self._stop = threading.Event()
        self._task_ready = threading.Event()
        self._threads_started: bool = False
        self._threads: list[threading.Thread] = []

        # Resilience state
        self._outage = _OutageState(max_unreachable=max_unreachable_seconds)
        self._disconnect_lock = threading.Lock()

        # Register SIO handlers up front
        if self.tsio is not None:
            self.tsio.on(ProviderRequest, self._on_provider_request)
            self.tsio.on(TaskAvailable, self._on_task_available)

            # Re-register SIO rooms after transport reconnection
            tsio = self.tsio  # capture for closure (pyright narrowing)

            def _on_reconnect() -> None:
                for job_name in list(self._registry):
                    tsio.emit(
                        JoinJobRoom(
                            job_name=job_name,
                            worker_id=str(self._worker_id),
                        )
                    )
                for full_name in list(self._providers):
                    tsio.emit(
                        JoinProviderRoom(
                            provider_name=full_name,
                            worker_id=str(self._worker_id),
                        )
                    )

            self.tsio.on("connect", _on_reconnect)

    # -- Container protocol --------------------------------------------------

    def __getitem__(self, key: str) -> type[Extension]:
        return self._registry[key][0]

    def __contains__(self, key: str) -> bool:
        return key in self._registry

    def __len__(self) -> int:
        return len(self._registry)

    def __iter__(self) -> Iterator[str]:
        return iter(self._registry)

    def get_run_kwargs(self, key: str) -> dict[str, Any]:
        """Return the run_kwargs stored for *key*."""
        return self._registry[key][1]

    def __enter__(self) -> "JobManager":
        return self

    def __exit__(self, *_exc_info) -> None:
        self.disconnect()

    # -- Lifecycle ------------------------------------------------------------

    def disconnect(self) -> None:
        """Gracefully disconnect the worker.

        Thread-safe — guarded by ``_disconnect_lock`` to prevent concurrent
        execution from signal handler + __exit__.
        """
        with self._disconnect_lock:
            if self._worker_id is None:
                self._stop.set()
                return

            self._stop.set()
            self._task_ready.set()

            for t in self._threads:
                t.join(timeout=10.0)
            self._threads.clear()

            worker_id = self._worker_id
            self._worker_id = None

            if self.tsio is not None:
                for full_name in self._providers:
                    self.tsio.emit(
                        LeaveProviderRoom(
                            provider_name=full_name,
                            worker_id=str(worker_id),
                        )
                    )

            if self.tsio is not None:
                for job_name in self._registry:
                    self.tsio.emit(
                        LeaveJobRoom(job_name=job_name, worker_id=str(worker_id))
                    )

            try:
                resp = self.api.http.delete(
                    f"{self.api.base_url}/v1/joblib/workers/{worker_id}",
                    headers=self.api.get_headers(),
                    timeout=5.0,
                )
                self.api.raise_for_status(resp)
            except Exception:  # noqa: BLE001
                logger.debug(
                    "Failed to delete worker %s (server unreachable?)",
                    worker_id,
                )

            self._providers.clear()
            self._registry.clear()

    def wait(self) -> None:
        """Block until ``disconnect()`` is called or a signal is received.

        When called from the main thread, installs SIGINT/SIGTERM handlers
        that trigger ``disconnect()`` and restores originals on exit.
        From other threads, simply blocks on the stop event.
        """
        is_main = threading.current_thread() is threading.main_thread()

        if is_main:
            original_sigint = signal.getsignal(signal.SIGINT)
            original_sigterm = signal.getsignal(signal.SIGTERM)

            def _shutdown(signum: int, _frame: Any) -> None:
                logger.info("Received signal %s, shutting down...", signum)
                self.disconnect()

            signal.signal(signal.SIGINT, _shutdown)
            signal.signal(signal.SIGTERM, _shutdown)

        try:
            self._stop.wait()
        finally:
            if is_main:
                signal.signal(signal.SIGINT, original_sigint)
                signal.signal(signal.SIGTERM, original_sigterm)

    @property
    def worker_id(self) -> UUID | None:
        """The worker ID for this manager.

        Set after create_worker() or first job registration.
        """
        return self._worker_id

    def _retry_startup(
        self,
        label: str,
        fn: Callable[[], httpx.Response],
    ) -> httpx.Response:
        """Execute *fn* with retry/backoff, re-raising fatal errors immediately.

        Parameters
        ----------
        label
            Human-readable label for log messages (e.g. ``"create_worker"``).
        fn
            Callable that performs the HTTP request and returns the response.
            ``raise_for_status`` is called on the response internally.
        """
        _rng = random.SystemRandom()
        for attempt in range(self._max_startup_retries):
            try:
                resp = fn()
                self.api.raise_for_status(resp)
            except (KeyError, PermissionError):
                raise
            except Exception as e:
                if attempt == self._max_startup_retries - 1:
                    raise
                delay = min(2**attempt, 30) + _rng.uniform(0, 1)
                logger.warning(
                    "%s failed (attempt %d/%d): %s",
                    label,
                    attempt + 1,
                    self._max_startup_retries,
                    e,
                )
                time.sleep(delay)
            else:
                return resp
        raise RuntimeError("Unreachable: retry loop completed without response")

    def create_worker(self) -> UUID:
        """Create a new worker and return its ID.

        The worker is linked to the authenticated user. Use this to explicitly
        create a worker before registering jobs.
        """
        resp = self._retry_startup(
            "create_worker",
            lambda: self.api.http.post(
                f"{self.api.base_url}/v1/joblib/workers",
                headers=self.api.get_headers(),
            ),
        )
        self._worker_id = UUID(resp.json()["id"])
        return self._worker_id

    # -- Job registration -----------------------------------------------------

    def register(
        self,
        extension_class: type[Extension] | None = None,
        *,
        room: str | None = None,
        run_kwargs: dict[str, Any] | None = None,
    ):
        """
        Register an Extension class as a job.

        Parameters
        ----------
        extension_class
            Extension subclass to register.
        room
            Room scope. Defaults to ``"@global"`` when called directly.
            ``ZnDraw.register_job()`` resolves *None* to the client's room.
        run_kwargs
            Extra keyword arguments passed to ``extension.run()`` at
            execution time. Useful for heavy, pre-loaded objects (e.g.
            torch models) that should not be re-created per task.

        Usage:
            @manager.register
            class Rotate(Extension):
                category: ClassVar[Category] = Category.MODIFIER
                angle: float = 0.0

            @manager.register(room="room_123")
            class PrivateRotate(Extension):
                category: ClassVar[Category] = Category.MODIFIER
                angle: float = 0.0
        """

        def decorator(cls: type[Extension]) -> type[Extension]:
            self._register_impl(cls, room, run_kwargs)
            return cls

        if extension_class is None:
            return decorator
        return decorator(extension_class)

    def _register_impl(
        self,
        cls: type[Extension],
        room: str | None,
        run_kwargs: dict[str, Any] | None = None,
    ) -> None:
        """Internal implementation of job registration."""
        room_id = room if room is not None else "@global"
        category = cls.category.value
        name = cls.__name__
        schema = cls.model_json_schema()

        request = JobRegisterRequest(
            category=category,
            name=name,
            schema=schema,
            worker_id=self._worker_id,
        )

        resp = self._retry_startup(
            "register",
            lambda: self.api.http.put(
                f"{self.api.base_url}/v1/joblib/rooms/{room_id}/jobs",
                headers=self.api.get_headers(),
                json=request.model_dump(exclude_none=True, mode="json"),
            ),
        )
        data = resp.json()
        full_name = f"{room_id}:{category}:{name}"
        if data.get("worker_id"):
            self._worker_id = UUID(data["worker_id"])
        if resp.status_code == 200:
            logger.info("Already registered: %s", full_name)
        self._registry[full_name] = (cls, run_kwargs or {})

        if self.tsio is not None:
            self.tsio.emit(
                JoinJobRoom(job_name=full_name, worker_id=str(self._worker_id))
            )
        self._ensure_background_threads()

    # -- Task operations ------------------------------------------------------

    def claim(self) -> ClaimedTask | None:
        """
        Attempt to claim a single task.

        Returns ClaimedTask with the Extension instance, or None if no tasks available.

        Raises:
            ValueError: If worker_id is not set.
                Call create_worker() or register a job first.
        """
        if self._worker_id is None:
            raise ValueError(
                "Worker ID not set. Call create_worker() or register a job first."
            )

        response = self.api.http.post(
            f"{self.api.base_url}/v1/joblib/tasks/claim",
            headers=self.api.get_headers(),
            json={"worker_id": str(self._worker_id)},
        )
        self.api.raise_for_status(response)
        claim_response = TaskClaimResponse.model_validate(response.json())

        if claim_response.task is None:
            return None

        task = claim_response.task
        job_name = task.job_name

        # Look up the Extension class from registry
        if job_name not in self._registry:
            raise KeyError(f"Job '{job_name}' not registered with this JobManager")

        extension_cls, run_kwargs = self._registry[job_name]
        extension = extension_cls.model_validate(task.payload)

        return ClaimedTask(
            task_id=str(task.id),
            job_name=job_name,
            room_id=task.room_id,
            extension=extension,
            run_kwargs=run_kwargs,
        )

    def listen(
        self,
        polling_interval: float = 2.0,
        stop_event: threading.Event | None = None,
    ) -> Iterator[ClaimedTask]:
        """
        Generator that yields claimed tasks indefinitely.

        Polls the server for tasks at the specified interval.
        Pass a threading.Event as stop_event for graceful shutdown.
        """
        while not (stop_event and stop_event.is_set()):
            claimed = self.claim()
            if claimed is not None:
                yield claimed
            else:
                time.sleep(polling_interval)

    def _update_task(
        self, task_id: str, status: TaskStatus, error: str | None = None
    ) -> None:
        """Update a task's status via the server."""
        body: dict[str, str] = {"status": status.value}
        if error is not None:
            body["error"] = error
        response = self.api.http.patch(
            f"{self.api.base_url}/v1/joblib/tasks/{task_id}",
            headers=self.api.get_headers(),
            json=body,
        )
        self.api.raise_for_status(response)

    def start(self, task: ClaimedTask) -> None:
        """Transition a claimed task to RUNNING."""
        self._update_task(task.task_id, TaskStatus.RUNNING)

    def complete(self, task: ClaimedTask) -> None:
        """Transition a running task to COMPLETED."""
        self._update_task(task.task_id, TaskStatus.COMPLETED)

    def fail(self, task: ClaimedTask, error: str) -> None:
        """Transition a running task to FAILED with an error message."""
        self._update_task(task.task_id, TaskStatus.FAILED, error)

    def fail_by_id(self, task_id: str, error: str) -> None:
        """Mark a task as failed by ID, without requiring a ClaimedTask."""
        self._update_task(task_id, TaskStatus.FAILED, error)

    def cancel(self, task: ClaimedTask) -> None:
        """Transition a claimed or running task to CANCELLED."""
        self._update_task(task.task_id, TaskStatus.CANCELLED)

    def heartbeat(self) -> None:
        """Send a heartbeat to keep the worker alive."""
        if self._worker_id is None:
            raise ValueError(
                "Worker ID not set. Call create_worker() or register a job first."
            )
        response = self.api.http.patch(
            f"{self.api.base_url}/v1/joblib/workers/{self._worker_id}",
            headers=self.api.get_headers(),
        )
        self.api.raise_for_status(response)

    def submit(
        self, extension: Extension, room: str, *, job_room: str = "@global"
    ) -> str:
        """
        Submit a task for processing.

        Args:
            extension: The Extension instance with parameters
            room: The room to submit the task to
            job_room: The room where the job is registered (default: @global)

        Returns:
            The task ID
        """
        category = extension.category.value
        name = extension.__class__.__name__
        job_name = f"{job_room}:{category}:{name}"

        # Build request using Pydantic model
        request = TaskSubmitRequest(payload=extension.model_dump())

        response = self.api.http.post(
            f"{self.api.base_url}/v1/joblib/rooms/{room}/tasks/{job_name}",
            headers=self.api.get_headers(),
            json=request.model_dump(),
        )
        self.api.raise_for_status(response)
        return response.json()["id"]

    # -- Provider operations --------------------------------------------------

    def register_provider(
        self,
        provider_cls: type[Provider],
        *,
        name: str,
        handler: Any,
        room: str = "@global",
    ) -> UUID:
        """Register a provider for serving read requests.

        Parameters
        ----------
        provider_cls
            The Provider subclass (defines category and read params schema).
        name
            Unique name for this provider instance (e.g., "local", "s3-bucket").
        handler
            The handler object passed to ``provider.read(handler)`` on dispatch.
        room
            ``"@global"`` or a room_id.

        Returns
        -------
        UUID
            The provider ID assigned by the server.
        """
        schema = provider_cls.model_json_schema()

        request = ProviderRegisterRequest(
            category=provider_cls.category,
            name=name,
            schema=schema,
            content_type=provider_cls.content_type,
            worker_id=self._worker_id,
        )

        resp = self._retry_startup(
            "register_provider",
            lambda: self.api.http.put(
                f"{self.api.base_url}/v1/joblib/rooms/{room}/providers",
                headers=self.api.get_headers(),
                json=request.model_dump(exclude_none=True, mode="json"),
            ),
        )
        data = resp.json()
        provider_id = UUID(data["id"])
        full_name = f"{room}:{provider_cls.category}:{name}"

        if data.get("worker_id"):
            self._worker_id = UUID(data["worker_id"])

        self._providers[full_name] = _RegisteredProvider(
            id=provider_id,
            cls=provider_cls,
            handler=handler,
            room_id=room,
        )

        if self.tsio is not None:
            self.tsio.emit(
                JoinProviderRoom(
                    provider_name=full_name,
                    worker_id=str(self._worker_id),
                )
            )

        self._ensure_background_threads()
        return provider_id

    def unregister_provider(self, full_name: str) -> None:
        """Unregister a provider by full_name (room_id:category:name)."""
        reg = self._providers.pop(full_name, None)
        if reg is None:
            return

        resp = self.api.http.delete(
            f"{self.api.base_url}/v1/joblib/providers/{reg.id}",
            headers=self.api.get_headers(),
        )
        self.api.raise_for_status(resp)

        if self.tsio is not None:
            self.tsio.emit(
                LeaveProviderRoom(
                    provider_name=full_name,
                    worker_id=str(self._worker_id),
                )
            )

    @property
    def handlers(self) -> dict[str, Any]:
        """All registered provider handlers keyed by full_name.

        Used by job execution — the full dict is passed to Extension.run()
        so the Extension can look up the right handler by provider full_name.
        """
        return {full_name: reg.handler for full_name, reg in self._providers.items()}

    # -- Background threads ---------------------------------------------------

    def _ensure_background_threads(self) -> None:
        """Start heartbeat + claim threads on first registration."""
        if self._threads_started:
            return
        self._threads_started = True

        t = threading.Thread(target=self._heartbeat_loop, daemon=True)
        t.start()
        self._threads.append(t)
        if self._execute is not None:
            t = threading.Thread(target=self._claim_loop, daemon=True)
            t.start()
            self._threads.append(t)

    def _heartbeat_loop(self) -> None:
        """Send periodic heartbeats until stopped."""
        while not self._stop.wait(self._heartbeat_interval):
            try:
                self.heartbeat()
                self._outage.record_success()
            except (KeyError, PermissionError):
                logger.exception("Registration lost, shutting down")
                self._stop.set()
                return
            except Exception as e:  # noqa: BLE001
                self._outage.record_failure()
                if self._outage.should_shutdown():
                    logger.error(  # noqa: TRY400 — intentionally no traceback
                        "Server unreachable for >%ss, shutting down. Last error: %s",
                        self._max_unreachable_seconds,
                        e,
                    )
                    self._stop.set()
                    return
                if self._outage.should_log():
                    logger.warning(
                        "Server unreachable — retrying (%.0fs/%.0fs elapsed)",
                        self._outage.elapsed(),
                        self._max_unreachable_seconds,
                    )

    def _claim_loop(self) -> None:
        """Claim and execute tasks until stopped."""
        backoff_attempt = 0
        while not self._stop.is_set():
            self._task_ready.clear()
            try:
                claimed = self.claim()
                self._outage.record_success()
                backoff_attempt = 0
            except (KeyError, PermissionError):
                logger.exception("Registration lost, shutting down")
                self._stop.set()
                return
            except Exception as e:  # noqa: BLE001
                self._outage.record_failure()
                if self._outage.should_shutdown():
                    logger.error(  # noqa: TRY400 — intentionally no traceback
                        "Server unreachable for >%ss, shutting down. Last error: %s",
                        self._max_unreachable_seconds,
                        e,
                    )
                    self._stop.set()
                    return
                if self._outage.should_log():
                    logger.warning(
                        "Server unreachable — retrying (%.0fs/%.0fs elapsed)",
                        self._outage.elapsed(),
                        self._max_unreachable_seconds,
                    )
                wait = min(self._polling_interval * 2**backoff_attempt, 10.0)
                backoff_attempt += 1
                self._task_ready.wait(timeout=wait)
                continue
            if claimed is not None:
                try:
                    self.start(claimed)
                    self._outage.record_success()
                except (KeyError, PermissionError):
                    logger.exception(
                        "Registration lost during task start, shutting down"
                    )
                    self._stop.set()
                    return
                except Exception as e:  # noqa: BLE001
                    logger.warning("Failed to start task %s: %s", claimed.task_id, e)
                    continue
                try:
                    self._execute(claimed)
                except Exception as e:  # noqa: BLE001
                    try:
                        self.fail(claimed, str(e))
                        self._outage.record_success()
                    except (KeyError, PermissionError):
                        logger.exception(
                            "Registration lost during task fail, shutting down"
                        )
                        self._stop.set()
                        return
                    except Exception:
                        logger.exception(
                            "Failed to mark task %s as failed", claimed.task_id
                        )
                    else:
                        logger.exception("Task %s failed", claimed.task_id)
                else:
                    try:
                        self.complete(claimed)
                        self._outage.record_success()
                    except (KeyError, PermissionError):
                        logger.exception(
                            "Registration lost during task completion, shutting down"
                        )
                        self._stop.set()
                        return
                    except Exception:
                        logger.exception(
                            "Failed to mark task %s completed", claimed.task_id
                        )
            else:
                self._task_ready.wait(timeout=self._polling_interval)

    # -- Logging helpers -------------------------------------------------------

    def _log_to_room(self, room_id: str, message: str) -> None:
        """Send a chat message to a room (best-effort, never raises)."""
        try:
            self.api.http.post(
                f"{self.api.base_url}/v1/rooms/{room_id}/chat/messages",
                json={"content": message},
                headers=self.api.get_headers(),
            )
        except Exception:  # noqa: BLE001
            logger.debug("Failed to log message to room %s", room_id, exc_info=True)

    # -- SIO event handlers ---------------------------------------------------

    def _on_task_available(self, _event: TaskAvailable) -> None:
        """Wake the claim loop when a new task is available."""
        self._task_ready.set()

    def _on_provider_request(self, event: ProviderRequest) -> None:
        """Handle incoming ProviderRequest dispatched by the server."""
        reg = self._providers.get(event.provider_name)
        if reg is None:
            logger.warning("Unknown provider: %s", event.provider_name)
            return

        try:
            params = json.loads(event.params)
            instance = reg.cls(**params)
            result = instance.read(reg.handler)
        except Exception:
            logger.exception(
                "Provider %s read failed for request %s",
                event.provider_name,
                event.request_id,
            )
            self._log_to_room(
                reg.room_id,
                f"Provider `{event.provider_name}` read failed for "
                f"request `{event.request_id}`:\n```\n{traceback.format_exc()}```",
            )
            return

        try:
            if reg.cls.content_type == "application/json":
                data = json.dumps(result).encode()
            else:
                data = result

            resp = self.api.http.post(
                f"{self.api.base_url}/v1/joblib/providers/{reg.id}/results",
                content=data,
                headers={**self.api.get_headers(), "X-Request-Hash": event.request_id},
            )
            self.api.raise_for_status(resp)
        except (KeyError, PermissionError):
            logger.exception("Registration lost during provider result upload")
            self._stop.set()
        except Exception:  # noqa: BLE001
            logger.warning(
                "Failed to upload provider result for %s", event.provider_name
            )
