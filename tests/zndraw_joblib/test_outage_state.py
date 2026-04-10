"""Unit tests for _OutageState retry tracker."""

import threading

import pytest

from zndraw_joblib.client import _OutageState


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
