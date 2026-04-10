# Worker Retry UX Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace noisy, traceback-heavy worker shutdown logs with a clean countdown and backoff during server outages.

**Architecture:** A shared `_OutageState` dataclass in `client.py` centralizes outage tracking, rate-limited logging, and shutdown decisions for both `_heartbeat_loop` and `_claim_loop`. Injectable clock enables deterministic unit testing.

**Tech Stack:** Python stdlib (`threading`, `time`, `dataclasses`), pytest + `caplog`

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `src/zndraw_joblib/client.py` | Modify | Add `_OutageState` dataclass, refactor both loops |
| `tests/zndraw_joblib/test_outage_state.py` | Create | Unit tests for `_OutageState` + integration tests for loop log output |

---

### Task 1: `_OutageState` Unit Tests + Implementation

**Files:**
- Create: `tests/zndraw_joblib/test_outage_state.py`
- Modify: `src/zndraw_joblib/client.py` (add `_OutageState` before `JobManager` class)

- [ ] **Step 1: Write the test file with all `_OutageState` unit tests**

Create `tests/zndraw_joblib/test_outage_state.py`:

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py -v`
Expected: `ImportError: cannot import name '_OutageState' from 'zndraw_joblib.client'`

- [ ] **Step 3: Implement `_OutageState` in `client.py`**

Add the following **above** the `JobManager` class definition (after the imports and helper classes, before line ~95 where `JobManager` begins). In `src/zndraw_joblib/client.py`:

```python
@dataclasses.dataclass
class _OutageState:
    """Tracks server outage for coordinated retry logging across loops."""

    max_unreachable: float
    clock: Callable[[], float] = dataclasses.field(default=time.monotonic)
    _outage_start: float | None = dataclasses.field(default=None, init=False)
    _last_log_time: float = dataclasses.field(default=0.0, init=False)
    _lock: threading.Lock = dataclasses.field(default_factory=threading.Lock, init=False)

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
```

Also add `import dataclasses` to the imports at the top of the file.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py -v`
Expected: all 9 tests PASS

- [ ] **Step 5: Commit**

```bash
git add tests/zndraw_joblib/test_outage_state.py src/zndraw_joblib/client.py
git commit -m "feat: add _OutageState with unit tests (issue #908)"
```

---

### Task 2: Refactor `JobManager.__init__` to Use `_OutageState`

**Files:**
- Modify: `src/zndraw_joblib/client.py:132-163` (`__init__`)

- [ ] **Step 1: Write a test that `JobManager` exposes `_outage`**

Append to `tests/zndraw_joblib/test_outage_state.py`:

```python
from unittest.mock import MagicMock

from zndraw_joblib.client import JobManager


class TestJobManagerOutageIntegration:
    def test_job_manager_has_outage_state(self):
        api = MagicMock()
        manager = JobManager(api, max_unreachable_seconds=60.0)
        assert isinstance(manager._outage, _OutageState)
        assert manager._outage.max_unreachable == 60.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py::TestJobManagerOutageIntegration::test_job_manager_has_outage_state -v`
Expected: `AttributeError: 'JobManager' object has no attribute '_outage'`

- [ ] **Step 3: Modify `JobManager.__init__`**

In `src/zndraw_joblib/client.py`, replace lines 160-162:

```python
        # Resilience state
        self._last_server_contact: float = time.monotonic()
        self._contact_lock = threading.Lock()
```

with:

```python
        # Resilience state
        self._outage = _OutageState(max_unreachable=max_unreachable_seconds)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py::TestJobManagerOutageIntegration -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/zndraw_joblib/client.py tests/zndraw_joblib/test_outage_state.py
git commit -m "refactor: wire _OutageState into JobManager.__init__ (issue #908)"
```

---

### Task 3: Refactor `_heartbeat_loop` to Use `_OutageState`

**Files:**
- Modify: `src/zndraw_joblib/client.py:698-728` (`_update_last_contact`, `_is_unreachable`, `_heartbeat_loop`)

- [ ] **Step 1: Write integration test for heartbeat loop log format**

Append to `tests/zndraw_joblib/test_outage_state.py`:

```python
import logging
import threading
import time


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
        manager.heartbeat = MagicMock(
            side_effect=ConnectionError("Connection refused")
        )
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py::TestHeartbeatLoopLogs -v`
Expected: FAIL (old log format doesn't match)

- [ ] **Step 3: Refactor `_heartbeat_loop`**

In `src/zndraw_joblib/client.py`, delete `_update_last_contact` (lines 698-701) and `_is_unreachable` (lines 703-708). Then replace `_heartbeat_loop` (lines 710-728) with:

```python
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
            except Exception as e:
                self._outage.record_failure()
                if self._outage.should_shutdown():
                    logger.error(
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
```

- [ ] **Step 4: Update all `_update_last_contact()` calls outside the loops temporarily**

Search for remaining `self._update_last_contact()` calls in `_claim_loop` and replace each with `self._outage.record_success()`. This is needed so the code compiles — the full `_claim_loop` refactor happens in Task 4. The calls to update are at lines (approximate after deletion):

- `self._update_last_contact()` after `self.claim()` succeeds
- `self._update_last_contact()` after `self.start(claimed)`
- `self._update_last_contact()` after `self.fail(claimed, ...)`
- `self._update_last_contact()` after `self.complete(claimed)`

Replace each `self._update_last_contact()` with `self._outage.record_success()`.

Also update `_is_unreachable()` calls in `_claim_loop` — replace `self._is_unreachable()` with `self._outage.should_shutdown()` temporarily (two occurrences in the claim except block). Keep the old `logger.warning`/`logger.exception` calls in `_claim_loop` for now — Task 4 will clean them up.

- [ ] **Step 5: Run all existing tests to verify nothing is broken**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py tests/zndraw_joblib/test_client.py -v`
Expected: all tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/zndraw_joblib/client.py tests/zndraw_joblib/test_outage_state.py
git commit -m "refactor: heartbeat loop uses _OutageState with countdown logs (issue #908)"
```

---

### Task 4: Refactor `_claim_loop` with Backoff and Clean Logging

**Files:**
- Modify: `src/zndraw_joblib/client.py` (`_claim_loop`)
- Modify: `tests/zndraw_joblib/test_outage_state.py` (add claim loop tests)

- [ ] **Step 1: Write integration test for claim loop clean shutdown**

Append to `tests/zndraw_joblib/test_outage_state.py`:

```python
class TestClaimLoopLogs:
    def test_claim_loop_clean_shutdown_no_traceback(self, caplog):
        """Claim loop should shut down cleanly with no raw tracebacks."""
        api = MagicMock()
        manager = JobManager(
            api,
            polling_interval=0.05,
            max_unreachable_seconds=0.3,
        )
        manager.claim = MagicMock(
            side_effect=ConnectionError("Connection refused")
        )
        manager._worker_id = MagicMock()
        manager._execute = lambda t: None  # enable claim loop

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
```

- [ ] **Step 2: Write test for outage recovery resetting retry state**

Append to `tests/zndraw_joblib/test_outage_state.py`:

```python
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
```

- [ ] **Step 3: Run tests to verify claim loop test fails**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py::TestClaimLoopLogs -v`
Expected: FAIL (old log format has tracebacks)

- [ ] **Step 4: Refactor `_claim_loop`**

Replace the `_claim_loop` method in `src/zndraw_joblib/client.py` with:

```python
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
            except Exception as e:
                self._outage.record_failure()
                if self._outage.should_shutdown():
                    logger.error(
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
```

- [ ] **Step 5: Run all tests**

Run: `uv run pytest tests/zndraw_joblib/test_outage_state.py tests/zndraw_joblib/test_client.py -v`
Expected: all tests PASS

- [ ] **Step 6: Commit**

```bash
git add src/zndraw_joblib/client.py tests/zndraw_joblib/test_outage_state.py
git commit -m "refactor: claim loop uses _OutageState with backoff (issue #908)"
```

---

### Task 5: Cleanup and Final Verification

**Files:**
- Modify: `src/zndraw_joblib/client.py` (remove dead code)

- [ ] **Step 1: Verify no remaining references to removed methods**

Run: `grep -n "_update_last_contact\|_is_unreachable\|_last_server_contact\|_contact_lock" src/zndraw_joblib/client.py`
Expected: no matches

- [ ] **Step 2: Run the full test suite**

Run: `uv run pytest tests/zndraw_joblib/ -v`
Expected: all tests PASS

- [ ] **Step 3: Run pre-commit checks**

Run: `uvx prek --all-files`
Expected: all checks pass (or only unrelated warnings)

- [ ] **Step 4: Commit any formatting fixes**

```bash
git add -u
git commit -m "style: pre-commit fixes (issue #908)"
```

(Skip this step if pre-commit made no changes.)
