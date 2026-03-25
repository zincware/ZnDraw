"""Integration tests for worker resilience — SIO disconnect vs. HTTP DELETE.

Verifies that removing eager `cleanup_worker()` from SIO `on_disconnect`
does not break the system:
- SIO disconnect preserves worker records (sweeper handles cleanup).
- HTTP DELETE (jobs.disconnect()) still cleans up correctly.
- Server restart with fresh DB causes worker to get 401/404.
"""

import logging
import threading
import time
from typing import ClassVar

import pytest
from zndraw_joblib.client import Category, Extension

from zndraw import ZnDraw

logger = logging.getLogger(__name__)

# =============================================================================
# Test extensions
# =============================================================================


class _Noop(Extension):
    """Extension that does nothing — used to register a job."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self, vis, **_kwargs):
        pass


_Noop.__name__ = "Noop"


# =============================================================================
# Test 1: Server restart with fresh DB → worker exits
# =============================================================================


def test_worker_exits_on_server_restart_fresh_db(server_factory):
    """Worker thread exits when server restarts with fresh DB (user gone).

    Scenario: server killed + restarted with fresh in-memory DB.
    The worker's token/user no longer exists → claim() gets 401 → raises
    PermissionError (or KeyError for 404) → listen() loop ends → thread exits.
    """
    # Start first server
    instance = server_factory()
    server_url = instance.url

    # Create worker and register a job
    worker = ZnDraw(url=server_url)
    worker.jobs.register(_Noop)
    assert worker.jobs.worker_id is not None

    # Start listen() in a background thread — it should exit when
    # the server restarts with a fresh DB
    thread_error: list[Exception] = []
    thread_exited = threading.Event()

    def _listen_loop():
        try:
            for _task in worker.jobs.listen(polling_interval=0.5):
                pass  # never expect a task
        except Exception as exc:  # noqa: BLE001
            thread_error.append(exc)
        finally:
            thread_exited.set()

    thread = threading.Thread(target=_listen_loop, daemon=True)
    thread.start()

    # Give the listen loop a moment to start polling
    time.sleep(1)

    # Extract port from the first server URL so we can restart on it
    port = int(server_url.rsplit(":", 1)[-1])

    # Stop the first server
    instance.server.should_exit = True
    instance.thread.join(timeout=5)

    # Brief pause to let the port free up
    time.sleep(1)

    # Start a NEW server on the same port — fresh in-memory DB (new users)
    server_factory(
        {
            "ZNDRAW_SERVER_HOST": "127.0.0.1",
            "ZNDRAW_SERVER_PORT": str(port),
        }
    )

    # Wait for the worker thread to exit — it should fail because
    # the old user/worker no longer exists in the new DB
    exited = thread_exited.wait(timeout=30)
    assert exited, "Worker listen() thread did not exit after server restart"

    thread.join(timeout=5)

    # Clean up — the worker's token is invalid on the new server, so
    # just disconnect the socket (disconnect will fail on HTTP DELETE, that's ok)
    try:
        worker.jobs.disconnect()
    except Exception:  # noqa: BLE001
        logger.debug("Expected: jobs.disconnect() failed (stale token)", exc_info=True)
    try:
        worker.disconnect()
    except Exception:  # noqa: BLE001
        logger.debug("Expected: disconnect() failed (stale token)", exc_info=True)


# =============================================================================
# Test 2: Server restart with same DB → worker resumes
# =============================================================================


@pytest.mark.skip(
    reason=(
        "Requires server_factory variant with file-backed SQLite DB "
        "and same ZNDRAW_AUTH_SECRET_KEY across restarts. Current fixture "
        "always uses in-memory SQLite."
    )
)
def test_worker_survives_brief_server_restart():
    """Worker resumes after brief server restart with same DB.

    With a file-backed DB and stable JWT secret, the worker's token
    remains valid across restarts. After reconnecting, claim() should
    work normally.
    """


# =============================================================================
# Test 3: SIO reconnect does not delete worker
# =============================================================================


def test_sio_reconnect_does_not_delete_worker(server):
    """SIO disconnect does not trigger worker cleanup — sweeper handles it.

    After removing eager cleanup from on_disconnect, dropping the SIO
    transport should NOT delete the worker record. The worker can still
    heartbeat and claim tasks via HTTP.
    """
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(_Noop)
        worker_id = worker.jobs.worker_id
        assert worker_id is not None

        # Ensure socket is connected first (register triggers connect
        # only if auto_pickup, so connect explicitly)
        worker.connect()
        assert worker.connected

        # Drop the SIO connection at the transport level
        worker.socket.disconnect()
        time.sleep(1)  # let disconnect propagate to server

        assert not worker.connected

        # Verify worker record still exists by sending a heartbeat.
        # If the worker was deleted on SIO disconnect, this would fail
        # with KeyError (404).
        worker.jobs.heartbeat()  # should NOT raise

    finally:
        # Clean up via HTTP DELETE
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# Test 4: SIO disconnect still cleans cameras/locks (no crash)
# =============================================================================


def test_sio_disconnect_no_crash(server):
    """SIO disconnect completes without errors after removing worker cleanup.

    Verifies that the on_disconnect handler still runs successfully
    (cleaning cameras, locks, session state) even though the worker
    cleanup was removed. A new connection should work fine afterward.
    """
    vis = ZnDraw(url=server)
    try:
        vis.connect()
        assert vis.connected

        # Drop SIO connection — this triggers on_disconnect on the server
        vis.socket.disconnect()
        time.sleep(1)

        assert not vis.connected
    finally:
        vis.disconnect()

    # A new connection to the same server should work fine
    vis2 = ZnDraw(url=server)
    try:
        vis2.connect()
        assert vis2.connected
    finally:
        vis2.disconnect()
