# ZnDraw Test Suite Guidelines

## Test Hierarchy

1. **Route tests** -- `ASGITransport` + real Redis + real SQLite.
   For HTTP request/response behavior. `MockSioServer` is acceptable
   for verifying broadcast emissions.

2. **Integration/E2E tests** -- Real uvicorn via `server_factory`.
   For flows crossing process boundaries: Socket.IO, CLI client,
   workers, StateFileSource resolution.

3. **Unit tests** -- Pure logic with no I/O. Mocking allowed only here,
   and only when better design can't avoid it.

## Fixture Rules

- All shared fixtures live in `conftest.py`.
- No per-file session/client/Redis fixture definitions.
- `helpers.py` functions are the ONLY way to create test data:
  - `create_test_user_in_db(session, email)` -> `(User, token)`
  - `create_test_room(session, user, description)` -> `Room`
  - `auth_header(token)` -> `dict`
  - `create_test_token(user)` -> `str`
- Fixture scoping: function-scoped for data, session-scoped only for
  truly stateless infrastructure.

## Mocking Rules

### Banned

- `lambda: None` or `lambda: (None)` for Redis
- `AsyncMock()` as a substitute for real Redis or result backends
- `patch("StateFileSource.__call__")` to hide token resolution
- `patch("httpx.Client")` to avoid real HTTP (use `server_factory`)
- `monkeypatch.setattr` on module-level functions to skip service
  interactions when a real server could be used instead

### Allowed

- `@patch("os.kill")` -- prevents real process termination
- `monkeypatch.setattr("webbrowser.open", ...)` -- prevents browser
- `monkeypatch.setattr("time.sleep", ...)` -- prevents test delays
- `MockSioServer` in route tests -- lightweight fake for broadcasts
- `monkeypatch.setenv` / `monkeypatch.delenv` -- environment isolation

### Last resort (requires inline comment explaining WHY)

- `monkeypatch.setattr("zndraw.cli.StateFile", ...)` -- only if the
  code cannot accept a StateFile parameter
- `patch("_is_pid_alive", ...)` -- only for pure-logic unit tests of
  StateFileSource decision logic, not integration tests

## Writing New Tests

1. Plan what to test: which tier (route / E2E / unit), which fixtures.
2. RED/GREEN/REFACTOR: write failing test first, then implement.
3. Review against these guidelines before committing.
