# tests/test_resilience.py
"""Red-phase tests for JobManager resilience features.

These tests verify the planned resilience behaviors:
- Fatal errors (404 / PermissionError) trigger immediate exit
- Transient errors (ConnectError) are tolerated up to max_unreachable_seconds
- disconnect() is idempotent, thread-safe, and tolerates server-gone
- Registration retries with backoff on transient errors
- Provider result upload failures are handled gracefully
"""

import asyncio
import threading
import time
from typing import Any, ClassVar
from unittest.mock import MagicMock

import httpx
import pytest
from sqlalchemy import text

from zndraw_joblib.client import (
    Category,
    ClaimedTask,
    Extension,
    JobManager,
)
from zndraw_joblib.events import ProviderRequest
from zndraw_joblib.provider import Provider

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


class _Ext(Extension):
    """Concrete extension for resilience tests."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self, vis: Any, **kwargs: Any) -> None:
        pass


class _SlowExt(Extension):
    """Extension whose run() sleeps to simulate slow work."""

    category: ClassVar[Category] = Category.MODIFIER
    sleep_seconds: float = 2.0

    def run(self, vis: Any, **kwargs: Any) -> None:
        time.sleep(self.sleep_seconds)


class _FsProvider(Provider):
    """Test provider for filesystem reads."""

    category: ClassVar[str] = "filesystem"
    path: str = "/"

    def read(self, handler: Any) -> Any:
        return handler.list_dir(self.path)


class _UnreachableClient:
    """Proxy that raises ConnectError on any method call."""

    def __getattr__(self, name: str) -> Any:
        def raise_connect(*args: Any, **kwargs: Any) -> Any:
            raise httpx.ConnectError("Simulated network outage")

        return raise_connect


class _ToggleableApi:
    """ApiManager that can simulate network outages.

    Implements the ``ApiManager`` protocol. When ``reachable`` is ``False``,
    any HTTP call raises ``httpx.ConnectError`` instead of reaching the
    test server.
    """

    def __init__(self, test_client: Any) -> None:
        self._client = test_client
        self._reachable = True
        self._unreachable = _UnreachableClient()

    @property
    def reachable(self) -> bool:
        return self._reachable

    @reachable.setter
    def reachable(self, value: bool) -> None:
        self._reachable = value

    @property
    def http(self) -> Any:
        if not self._reachable:
            return self._unreachable
        return self._client

    @property
    def base_url(self) -> str:
        return ""

    def get_headers(self) -> dict[str, str]:
        return {}

    def raise_for_status(self, response: Any) -> None:
        """Conform to ApiManager exception contract.

        Raises
        ------
        KeyError
            404 responses.
        PermissionError
            401/403 responses.
        ValueError
            409/422 responses.
        """
        if response.status_code < 400:
            return
        try:
            detail = response.json().get("detail", response.text)
        except Exception:
            detail = response.text
        if response.status_code == 404:
            raise KeyError(str(detail))
        if response.status_code in {401, 403}:
            raise PermissionError(str(detail))
        if response.status_code in {409, 422}:
            raise ValueError(str(detail))
        response.raise_for_status()


# ===================================================================
# Fatal error tests
# ===================================================================


def test_heartbeat_404_triggers_exit(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Heartbeat loop exits when server returns 404 (worker deleted).

    Registers a worker, starts the heartbeat loop, deletes the worker
    directly from the DB, then asserts ``_stop`` is set within two
    heartbeat intervals.
    """
    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=0.5,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class HB404(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    worker_id = manager.worker_id
    assert worker_id is not None

    # Delete the worker directly from DB
    async def _delete() -> None:
        from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker

        factory = async_sessionmaker(
            threadsafe_engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            await session.execute(
                text("DELETE FROM worker WHERE id = :wid"),
                {"wid": worker_id.hex},
            )
            await session.commit()

    loop = asyncio.new_event_loop()
    loop.run_until_complete(_delete())
    loop.close()

    # _stop should be set within 2 heartbeat intervals (1.0s) + tolerance
    assert manager._stop.wait(timeout=3.0), (
        "_stop was not set after worker deletion (heartbeat 404)"
    )


def test_claim_404_triggers_exit(mock_client_api, threadsafe_client, threadsafe_engine):
    """Claim loop exits when server returns 404 (worker deleted).

    Registers a worker with an execute callback so the claim loop runs,
    deletes the worker from the DB, then asserts ``_stop`` is set within
    two polling intervals.
    """
    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        execute=lambda t: None,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class Claim404(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    worker_id = manager.worker_id
    assert worker_id is not None

    # Delete the worker directly from DB
    async def _delete() -> None:
        from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker

        factory = async_sessionmaker(
            threadsafe_engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            await session.execute(
                text("DELETE FROM worker WHERE id = :wid"),
                {"wid": worker_id.hex},
            )
            await session.commit()

    loop = asyncio.new_event_loop()
    loop.run_until_complete(_delete())
    loop.close()

    # _stop should be set within 2 polling intervals (0.2s) + tolerance
    assert manager._stop.wait(timeout=3.0), (
        "_stop was not set after worker deletion (claim 404)"
    )


def test_complete_404_triggers_exit(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Task completion exits gracefully when worker is deleted mid-task.

    Registers a worker with a slow execute callback, submits a task, then
    deletes the worker while the task is being processed. ``_stop`` should
    be set once the task completion attempt fails.
    """
    executed = threading.Event()

    def slow_execute(task: ClaimedTask) -> None:
        executed.set()
        time.sleep(1.0)

    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        execute=slow_execute,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class Complete404(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    worker_id = manager.worker_id
    assert worker_id is not None

    # Submit a task
    threadsafe_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Complete404",
        json={"payload": {}},
    )

    # Wait for task to start executing
    assert executed.wait(timeout=5.0), "Task never started executing"

    # Delete the worker while task is in progress
    async def _delete() -> None:
        from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker

        factory = async_sessionmaker(
            threadsafe_engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            await session.execute(
                text("DELETE FROM worker WHERE id = :wid"),
                {"wid": worker_id.hex},
            )
            await session.commit()

    loop = asyncio.new_event_loop()
    loop.run_until_complete(_delete())
    loop.close()

    # _stop should eventually be set (after the task finishes or fails)
    assert manager._stop.wait(timeout=10.0), (
        "_stop was not set after worker deletion during task execution"
    )


def test_heartbeat_permission_error_triggers_exit(
    mock_client_api, client_factory, threadsafe_engine
):
    """Heartbeat loop exits on 403 (PermissionError).

    Registers a worker as user A, then swaps the api to user B's client.
    The next heartbeat returns 403 because worker belongs to user A.
    """
    # Build two apps using the threadsafe engine
    user_a_client = client_factory("user_a")
    user_b_client = client_factory("user_b")

    api_a = mock_client_api(user_a_client)
    api_b = mock_client_api(user_b_client)

    manager = JobManager(
        api_a,
        heartbeat_interval=0.5,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class PermJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    worker_id = manager.worker_id
    assert worker_id is not None

    # Swap the api to user B — heartbeat will now 403
    manager.api = api_b

    # _stop should be set within 2 heartbeat intervals (1.0s) + tolerance
    assert manager._stop.wait(timeout=3.0), (
        "_stop was not set after 403 PermissionError on heartbeat"
    )


# ===================================================================
# Transient error tests
# ===================================================================


def test_network_blip_does_not_trigger_exit(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Short network outage does not stop the manager.

    Toggles unreachable for 3s (below max_unreachable_seconds=10),
    then toggles back. Asserts ``_stop`` is NOT set after recovery.
    """
    api = _ToggleableApi(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=0.5,
        polling_interval=0.1,
        max_unreachable_seconds=10.0,
    )

    @manager.register
    class BlipJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    # Simulate a 3-second network outage
    api.reachable = False
    time.sleep(3.0)
    api.reachable = True

    # Wait 5s after recovery — _stop should NOT be set
    time.sleep(5.0)
    assert not manager._stop.is_set(), "_stop was set during a short network blip"
    manager.disconnect()


def test_prolonged_unreachability_triggers_exit(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Prolonged network outage beyond max_unreachable_seconds triggers exit.

    Toggles unreachable permanently and asserts ``_stop`` is set within
    approximately max_unreachable_seconds.
    """
    api = _ToggleableApi(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=0.5,
        polling_interval=0.1,
        max_unreachable_seconds=3.0,
    )

    @manager.register
    class UnreachJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    # Toggle permanently unreachable
    api.reachable = False

    # _stop should be set within ~max_unreachable_seconds + tolerance
    assert manager._stop.wait(timeout=6.0), (
        "_stop was not set after prolonged unreachability"
    )


# ===================================================================
# Disconnect tests
# ===================================================================


def test_disconnect_tolerates_server_gone(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """disconnect() does not raise when the server is unreachable.

    Registers a worker, toggles unreachable, calls disconnect().
    No exception should be raised.
    """
    api = _ToggleableApi(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class GoneJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    assert manager.worker_id is not None

    # Server goes away
    api.reachable = False

    # disconnect() should not raise
    manager.disconnect()


def test_disconnect_still_cleans_local_state(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """disconnect() clears local state even when server DELETE fails.

    After disconnect with unreachable server, ``worker_id`` should be None,
    ``_registry`` empty, and ``_providers`` empty.
    """
    api = _ToggleableApi(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class CleanJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    mock_handler = MagicMock()
    mock_handler.list_dir.return_value = []
    manager.register_provider(_FsProvider, name="local", handler=mock_handler)

    assert manager.worker_id is not None
    assert len(manager) > 0
    assert len(manager._providers) > 0

    # Server goes away
    api.reachable = False

    # disconnect() should clean up local state regardless
    manager.disconnect()

    assert manager.worker_id is None
    assert len(manager._registry) == 0
    assert len(manager._providers) == 0


def test_disconnect_concurrent_calls_safe(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Concurrent disconnect() calls complete without exception.

    Calls disconnect() from two threads simultaneously. Both should
    complete without raising.
    """
    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class ConcJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    errors: list[Exception] = []

    def safe_disconnect() -> None:
        try:
            manager.disconnect()
        except Exception as exc:
            errors.append(exc)

    t1 = threading.Thread(target=safe_disconnect)
    t2 = threading.Thread(target=safe_disconnect)
    t1.start()
    t2.start()
    t1.join(timeout=10.0)
    t2.join(timeout=10.0)

    assert not errors, f"Concurrent disconnect raised: {errors}"


# ===================================================================
# Startup retry tests
# ===================================================================


def test_register_retries_on_connect_error(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """register() retries on transient ConnectError and succeeds.

    Starts with unreachable API. A background thread toggles reachable
    after 2s. register() should retry and eventually succeed.
    """
    api = _ToggleableApi(threadsafe_client)
    api.reachable = False

    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
        max_startup_retries=5,
    )

    def toggle_reachable() -> None:
        time.sleep(2.0)
        api.reachable = True

    t = threading.Thread(target=toggle_reachable, daemon=True)
    t.start()

    # register() should retry and eventually succeed
    @manager.register
    class RetryJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    assert "@global:modifiers:RetryJob" in manager
    manager.disconnect()


def test_register_provider_retries_on_connect_error(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """register_provider() retries on transient ConnectError.

    Same pattern: unreachable initially, toggled after 2s.
    """
    api = _ToggleableApi(threadsafe_client)

    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
        max_startup_retries=5,
    )

    # First register a job to get a worker_id (while reachable)
    @manager.register
    class PreJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    # Now toggle unreachable
    api.reachable = False

    def toggle_reachable() -> None:
        time.sleep(2.0)
        api.reachable = True

    t = threading.Thread(target=toggle_reachable, daemon=True)
    t.start()

    mock_handler = MagicMock()
    mock_handler.list_dir.return_value = []

    # register_provider() should retry and succeed
    provider_id = manager.register_provider(
        _FsProvider, name="retry_local", handler=mock_handler
    )
    assert provider_id is not None
    manager.disconnect()


def test_register_fails_fast_on_404(mock_client_api, threadsafe_client):
    """register() raises KeyError immediately on 404 (no retry).

    Passes a non-existent worker_id. The server returns 404, which
    should be raised as KeyError without retrying.
    """
    import uuid

    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
        max_startup_retries=5,
    )

    # Set a non-existent worker_id
    manager._worker_id = uuid.uuid4()

    with pytest.raises(KeyError):

        @manager.register
        class NoSuchWorker(_Ext):
            category: ClassVar[Category] = Category.MODIFIER


def test_register_exhausts_retries(mock_client_api, threadsafe_client):
    """register() raises after exhausting max_startup_retries.

    With permanently unreachable API and max_startup_retries=2,
    register() should raise after 2 attempts.
    """
    api = _ToggleableApi(threadsafe_client)
    api.reachable = False

    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
        max_startup_retries=2,
    )

    with pytest.raises(httpx.ConnectError):

        @manager.register
        class ExhaustedJob(_Ext):
            category: ClassVar[Category] = Category.MODIFIER


# ===================================================================
# Provider result upload tests
# ===================================================================


def test_provider_result_upload_failure_does_not_crash(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Provider result POST failure on transient error does not crash or stop.

    Registers a provider, toggles unreachable, then manually triggers a
    ProviderRequest event. The handler should not raise, and ``_stop``
    should NOT be set.
    """
    api = _ToggleableApi(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    mock_handler = MagicMock()
    mock_handler.list_dir.return_value = [{"name": "file.xyz"}]

    manager.register_provider(_FsProvider, name="local", handler=mock_handler)

    # Toggle unreachable
    api.reachable = False

    # Trigger a provider request manually
    event = ProviderRequest.from_dict_params(
        request_id="hash123",
        provider_name="@global:filesystem:local",
        params={"path": "/data"},
    )
    # Should not raise
    manager._on_provider_request(event)

    # _stop should NOT be set (transient error)
    assert not manager._stop.is_set(), (
        "_stop was set on transient provider upload failure"
    )
    manager.disconnect()


def test_provider_result_upload_404_triggers_exit(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Provider result POST 404 (fatal) triggers _stop.

    Registers a provider, deletes the worker from DB, then triggers a
    ProviderRequest event. The 404 on result upload should set ``_stop``.
    """
    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    mock_handler = MagicMock()
    mock_handler.list_dir.return_value = [{"name": "file.xyz"}]

    manager.register_provider(_FsProvider, name="local", handler=mock_handler)

    worker_id = manager.worker_id
    assert worker_id is not None

    # Delete the provider from DB so the provider result POST will 404
    provider_reg = manager._providers["@global:filesystem:local"]

    async def _delete() -> None:
        from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker

        factory = async_sessionmaker(
            threadsafe_engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            await session.execute(
                text("DELETE FROM provider WHERE id = :pid"),
                {"pid": provider_reg.id.hex},
            )
            await session.commit()

    loop = asyncio.new_event_loop()
    loop.run_until_complete(_delete())
    loop.close()

    # Trigger a provider request manually
    event = ProviderRequest.from_dict_params(
        request_id="hash456",
        provider_name="@global:filesystem:local",
        params={"path": "/data"},
    )
    manager._on_provider_request(event)

    # _stop should be set (fatal 404 error)
    assert manager._stop.is_set(), (
        "_stop was not set after 404 on provider result upload"
    )


# ===================================================================
# Interaction tests
# ===================================================================


def test_stop_unblocks_wait(mock_client_api, threadsafe_client, threadsafe_engine):
    """wait() unblocks when heartbeat detects deleted worker.

    Registers a worker with heartbeat loop, deletes the worker from DB,
    runs wait() in a separate thread. The thread should unblock within 5s.
    """
    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        heartbeat_interval=0.5,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class WaitJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    worker_id = manager.worker_id
    assert worker_id is not None

    unblocked = threading.Event()

    def wait_thread() -> None:
        manager.wait()
        unblocked.set()

    t = threading.Thread(target=wait_thread, daemon=True)
    t.start()

    # Delete the worker from DB
    async def _delete() -> None:
        from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker

        factory = async_sessionmaker(
            threadsafe_engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            await session.execute(
                text("DELETE FROM worker WHERE id = :wid"),
                {"wid": worker_id.hex},
            )
            await session.commit()

    loop = asyncio.new_event_loop()
    loop.run_until_complete(_delete())
    loop.close()

    # wait() should unblock within 5s (heartbeat detects 404 → _stop set)
    assert unblocked.wait(timeout=5.0), "wait() did not unblock after worker deletion"


def test_task_in_progress_then_exit(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Background task in progress when worker is deleted.

    Submits a slow task, deletes the worker while it is executing.
    ``_stop`` should be set. The task should reach completed or failed.
    """
    executed = threading.Event()
    task_ids: list[str] = []

    def slow_execute(task: ClaimedTask) -> None:
        task_ids.append(task.task_id)
        executed.set()
        time.sleep(2.0)

    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        execute=slow_execute,
        heartbeat_interval=0.5,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class SlowJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    worker_id = manager.worker_id
    assert worker_id is not None

    # Submit a task
    threadsafe_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:SlowJob",
        json={"payload": {}},
    )

    # Wait for task to start
    assert executed.wait(timeout=5.0), "Task never started"

    # Delete the worker from DB while task is executing
    async def _delete() -> None:
        from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker

        factory = async_sessionmaker(
            threadsafe_engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            await session.execute(
                text("DELETE FROM worker WHERE id = :wid"),
                {"wid": worker_id.hex},
            )
            await session.commit()

    loop = asyncio.new_event_loop()
    loop.run_until_complete(_delete())
    loop.close()

    # _stop should be set (heartbeat detects 404)
    assert manager._stop.wait(timeout=10.0), (
        "_stop was not set after worker deletion during task execution"
    )


# ===================================================================
# start() fall-through regression
# ===================================================================


def test_start_failure_skips_execute(
    mock_client_api, threadsafe_client, threadsafe_engine
):
    """Transient start() failure must NOT execute the task.

    If start() (CLAIMED -> RUNNING transition) fails with a transient error
    (e.g. ConnectError), the claim loop must skip execution entirely. The task
    stays CLAIMED and can be reclaimed later.

    Regression test for: start() fall-through bug where _execute() ran despite
    start() raising a transient exception.
    """
    executed = threading.Event()

    def track_execute(task: ClaimedTask) -> None:
        executed.set()

    api = mock_client_api(threadsafe_client)
    manager = JobManager(
        api,
        execute=track_execute,
        heartbeat_interval=30.0,
        polling_interval=0.1,
        max_unreachable_seconds=120.0,
    )

    @manager.register
    class StartFailJob(_Ext):
        category: ClassVar[Category] = Category.MODIFIER

    # Monkeypatch start() to always raise a transient error
    def failing_start(task: ClaimedTask) -> None:
        raise httpx.ConnectError("Simulated start() failure")

    manager.start = failing_start  # type: ignore[assignment]

    # Submit a task
    resp = threadsafe_client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:StartFailJob",
        json={"payload": {}},
    )
    assert resp.status_code in {201, 202}

    # Wait enough time for the claim loop to pick up and attempt start()
    time.sleep(2.0)

    # The execute callback must NOT have been called
    assert not executed.is_set(), (
        "execute() was called despite start() raising a transient error"
    )
    manager.disconnect()
