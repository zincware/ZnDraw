# Worker Shutdown Retry UX — Design Spec

**Issue:** [zincware/ZnDraw#908](https://github.com/zincware/ZnDraw/issues/908)
**Date:** 2026-04-10

## Problem

When the ZnDraw server goes down while a worker is connected, `_heartbeat_loop` (30s) and `_claim_loop` (2s) independently log failures with no countdown, no backoff, and raw tracebacks on exit. Over 120s this produces 60+ noisy, unactionable log lines.

## Approach: Centralized `_OutageState`

A shared dataclass tracks outage timing and coordinates logging across both loops. Replaces `_last_server_contact`, `_contact_lock`, `_update_last_contact()`, and `_is_unreachable()`.

### `_OutageState` (dataclass, private to `client.py`)

- **Fields**: `max_unreachable`, `clock` (injectable, defaults to `time.monotonic`), `_outage_start` (None when healthy), `_last_log_time`, `_lock`
- **`record_failure()`**: Sets `_outage_start` on first failure; subsequent calls are no-ops (don't move the start forward)
- **`record_success()`**: Resets `_outage_start` to None
- **`elapsed()`**: Seconds since first failure, or 0.0 if not in outage
- **`should_shutdown()`**: `elapsed() > max_unreachable`
- **`should_log(min_interval=5.0)`**: Rate-limits log output; returns True at most once per `min_interval` seconds

### Loop Changes

Both loops' except blocks become:

```python
except Exception as e:
    self._outage.record_failure()
    if self._outage.should_shutdown():
        logger.error(
            "Server unreachable for >%ss, shutting down. Last error: %s",
            self._max_unreachable_seconds, e,
        )
        self._stop.set()
        return
    if self._outage.should_log():
        elapsed = self._outage.elapsed()
        logger.warning(
            "Server unreachable — retrying (%.0fs/%.0fs elapsed)",
            elapsed, self._max_unreachable_seconds,
        )
```

On success, both loops call `self._outage.record_success()`.

**Backoff in `_claim_loop`**: During active outage, wait increases: `min(polling_interval * 2 ** attempt, 10.0)`. Resets to `polling_interval` when outage clears. Heartbeat loop keeps its 30s interval (already slow enough).

**Removed**: `_last_server_contact`, `_contact_lock`, `_update_last_contact()`, `_is_unreachable()`.

**Unchanged**: `KeyError`/`PermissionError` -> "Registration lost" path, loop cadences, `_stop` event coordination.

## Expected Log Output

```
[WARNING] Server unreachable — retrying (5s/120s elapsed)
[WARNING] Server unreachable — retrying (15s/120s elapsed)
[WARNING] Server unreachable — retrying (30s/120s elapsed)
[ERROR]   Server unreachable for >120s, shutting down. Last error: Connection refused
```

## Test Strategy

### Unit Tests (`tests/zndraw_joblib/test_outage_state.py`)

All use an injectable clock — no server, threads, or network needed.

1. **`test_initial_state`** — `elapsed() == 0`, `should_shutdown() == False`
2. **`test_failure_starts_outage_clock`** — `record_failure()` + advance clock -> `elapsed()` correct
3. **`test_success_resets_outage`** — `record_failure()` -> `record_success()` -> `elapsed() == 0`
4. **`test_should_shutdown_after_max`** — clock past max -> `True`
5. **`test_should_shutdown_before_max`** — clock under max -> `False`
6. **`test_should_log_rate_limiting`** — first call `True`, second within interval `False`, after interval `True`
7. **`test_multiple_failures_dont_reset_start`** — second `record_failure()` doesn't move `_outage_start`

### Integration Tests (log output verification)

8. **`test_heartbeat_loop_logs_countdown_format`** — `JobManager` with failing heartbeat endpoint, `max_unreachable_seconds=3`, assert `caplog` matches countdown format
9. **`test_claim_loop_clean_shutdown_no_traceback`** — same setup, assert no `"Traceback"` or `"httpx.ConnectError"` in captured logs
10. **`test_outage_recovery_resets_retry`** — fail -> fail -> succeed -> outage resets, next failure starts fresh

### Not in Scope

- Mocking Redis
- Actual server-down integration tests (covered by acceptance criteria, separate)
- Socketio exit path cleanup (separate concern)

## Acceptance Criteria Mapping

| Criterion | Addressed by |
|---|---|
| Retry logging includes elapsed/max time | `should_log()` + loop except block format |
| Final shutdown uses `logger.error()` not `logger.exception()` | Loop except block change |
| Backoff during sustained outage (cap 10s) | `_claim_loop` exponential backoff |
| No raw `httpx.ConnectError` traceback | `logger.error()` with `%s` formatting |
| Unit test: log output format | Tests 8-9 |
| Integration test: clean exit | Tests 8-9 (caplog-based) |
