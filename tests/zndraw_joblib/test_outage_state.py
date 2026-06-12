"""Unit tests for _OutageState retry tracker."""

import logging
import threading
from unittest.mock import MagicMock

import pytest

from zndraw_joblib.client import JobManager, _OutageState


def _make_clock(start: float = 0.0):
    """Return a (clock, advance) pair for deterministic time control."""
    t = [start]

    def clock() -> float:
        return t[0]

    def advance(dt: float) -> None:
        t[0] += dt

    return clock, advance


class TestOutageStateInitial:
    def test_elapsed_is_zero(self):
        clock, _ = _make_clock()
        state = _OutageState(max_unreachable=120.0, clock=clock)
        assert state.elapsed() == 0.0

    def test_should_shutdown_is_false(self):
        clock, _ = _make_clock()
        state = _OutageState(max_unreachable=120.0, clock=clock)
        assert state.should_shutdown() is False


class TestOutageStateFailure:
    def test_failure_starts_outage_clock(self):
        clock, advance = _make_clock()
        state = _OutageState(max_unreachable=120.0, clock=clock)
        state.record_failure()
        advance(10.0)
        assert state.elapsed() == pytest.approx(10.0)

    def test_multiple_failures_dont_reset_start(self):
        clock, advance = _make_clock()
        state = _OutageState(max_unreachable=120.0, clock=clock)
        state.record_failure()
        advance(5.0)
        state.record_failure()  # should NOT reset
        advance(5.0)
        assert state.elapsed() == pytest.approx(10.0)

    def test_success_resets_outage(self):
        clock, advance = _make_clock()
        state = _OutageState(max_unreachable=120.0, clock=clock)
        state.record_failure()
        advance(10.0)
        state.record_success()
        assert state.elapsed() == 0.0
        assert state.should_shutdown() is False


class TestOutageStateShutdown:
    def test_should_shutdown_after_max(self):
        clock, advance = _make_clock()
        state = _OutageState(max_unreachable=10.0, clock=clock)
        state.record_failure()
        advance(11.0)
        assert state.should_shutdown() is True

    def test_should_not_shutdown_before_max(self):
        clock, advance = _make_clock()
        state = _OutageState(max_unreachable=10.0, clock=clock)
        state.record_failure()
        advance(9.0)
        assert state.should_shutdown() is False


class TestOutageStateLogging:
    def test_should_log_first_call_returns_true(self):
        clock, _ = _make_clock()
        state = _OutageState(max_unreachable=120.0, clock=clock)
        assert state.should_log(min_interval=5.0) is True

    def test_should_log_rate_limits(self):
        clock, advance = _make_clock()
        state = _OutageState(max_unreachable=120.0, clock=clock)
        assert state.should_log(min_interval=5.0) is True
        advance(2.0)
        assert state.should_log(min_interval=5.0) is False
        advance(4.0)  # total 6s since last log
        assert state.should_log(min_interval=5.0) is True


class TestJobManagerOutageIntegration:
    def test_job_manager_has_outage_state(self):
        api = MagicMock()
        manager = JobManager(api, max_unreachable_seconds=60.0)
        assert isinstance(manager._outage, _OutageState)
        assert manager._outage.max_unreachable == 60.0


class TestHeartbeatLoopLogs:
    def test_heartbeat_loop_logs_countdown_format(self, caplog):
        """Heartbeat loop should log countdown and exit cleanly."""
        api = MagicMock()
        manager = JobManager(
            api,
            heartbeat_interval=0.05,
            max_unreachable_seconds=0.2,
        )
        # Make heartbeat always raise
        manager.heartbeat = MagicMock(side_effect=ConnectionError("Connection refused"))
        # Need a worker_id set so we don't hit other errors
        manager._worker_id = MagicMock()

        with caplog.at_level(logging.WARNING, logger="zndraw_joblib.client"):
            manager._stop.clear()
            t = threading.Thread(target=manager._heartbeat_loop, daemon=True)
            t.start()
            t.join(timeout=5.0)

        # Should contain countdown-style warning
        warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
        assert any("retrying" in r.message and "elapsed" in r.message for r in warnings)

        # Final message should be ERROR, not EXCEPTION (no traceback)
        errors = [r for r in caplog.records if r.levelno == logging.ERROR]
        assert any("shutting down" in r.message for r in errors)

        # No traceback lines in any log record
        for record in caplog.records:
            assert "Traceback" not in record.message
            assert record.exc_info is None or record.exc_info == (None, None, None)


class TestClaimLoopLogs:
    def test_claim_loop_clean_shutdown_no_traceback(self, caplog):
        """Claim loop should shut down cleanly with no raw tracebacks."""
        api = MagicMock()
        manager = JobManager(
            api,
            polling_interval=0.05,
            max_unreachable_seconds=0.3,
        )
        manager.claim = MagicMock(side_effect=ConnectionError("Connection refused"))
        manager._worker_id = MagicMock()
        manager._execute = lambda _task: None  # enable claim loop

        with caplog.at_level(logging.WARNING, logger="zndraw_joblib.client"):
            manager._stop.clear()
            t = threading.Thread(target=manager._claim_loop, daemon=True)
            t.start()
            t.join(timeout=5.0)

        # Should have shutdown error
        errors = [r for r in caplog.records if r.levelno == logging.ERROR]
        assert any("shutting down" in r.message for r in errors)

        # No traceback in any record
        for record in caplog.records:
            assert "Traceback" not in record.message
            assert record.exc_info is None or record.exc_info == (None, None, None)


class TestOutageRecovery:
    def test_outage_recovery_resets_retry(self):
        """After recovery, a new outage should start fresh."""
        clock, advance = _make_clock()
        state = _OutageState(max_unreachable=10.0, clock=clock)

        # First outage
        state.record_failure()
        advance(5.0)
        assert state.elapsed() == pytest.approx(5.0)

        # Recovery
        state.record_success()
        assert state.elapsed() == 0.0

        # New outage starts fresh
        advance(2.0)
        state.record_failure()
        advance(3.0)
        assert state.elapsed() == pytest.approx(3.0)
        assert state.should_shutdown() is False
