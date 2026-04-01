# Test Suite Overhaul — Design Spec

**Date:** 2026-03-26
**Status:** Draft
**Trigger:** `uv run zndraw file.h5` broken with 401 — no E2E test covered it. Audit revealed systemic testing issues.

## Problem Statement

The test suite has 500+ tests but missed a critical bug because:

1. **Mocked services hide real failures.** 15 route test files override Redis with `AsyncMock()` objects that silently swallow all calls. The shared `conftest.py` `client` fixture uses `lambda: (None)` for Redis. The `StateFileSource` is patched to `return_value={}` in client tests. These patterns prevent tests from catching real integration failures.
2. **Massive fixture duplication.** 15+ route test files each define identical session/client/Redis fixtures with different prefixes (~3000 lines of copy-paste). Helper functions (`_create_user`, `_create_room`) are redefined in 10+ files despite existing in `helpers.py`.
3. **Missing E2E tests.** No test exercises the main user flow (`uv run zndraw file.h5`). No test verifies ZnDraw client persistence across disconnect/reconnect. No test exercises guest auth → room operations end-to-end.
4. **Excessive monkeypatching.** ~50 `monkeypatch.setattr` / `@patch` calls plus ~30 `AsyncMock()` Redis overrides patch module-level functions to avoid testing real interactions. Many could be eliminated by using real servers.

## Constraints

- **Redis is always available** — required for all tests, no mock fallback.
- **SQLite in-memory is acceptable** for the DB layer — SQLModel handles backend abstraction.
- **ASGITransport for route tests** — fast, with real Redis and real DB.
- **Real uvicorn (`server_factory`) for E2E/integration** — Socket.IO, workers, CLI client, cross-process flows.
- **Mocking only for side effects** — `webbrowser.open`, `time.sleep`, `os.kill`. Never for services under test.
- **Big-bang PR** — all changes land in a single PR. Phases represent commit ordering within the branch.
- **Class-to-function conversion is out of scope** — cosmetic, no quality impact.
- **Never modify tests marked `@pytest.mark.protected`.**

## Design

### Phase 1: Testing Guidelines (`tests/README.md`)

A `tests/README.md` co-located with the test suite. First commit in the PR so all subsequent changes have a reference.

#### Test Hierarchy

1. **Route tests** — `ASGITransport` + real Redis + real SQLite. For HTTP request/response behavior. `MockSioServer` is acceptable here for verifying broadcast emissions.
2. **Integration/E2E tests** — Real uvicorn via `server_factory`. For flows crossing process boundaries: Socket.IO, CLI client, workers, StateFileSource resolution.
3. **Unit tests** — Pure logic with no I/O. Mocking allowed only here, and only when better design can't avoid it.

#### Mocking Rules

**Banned:**
- `lambda: None` or `lambda: (None)` for Redis
- `AsyncMock()` (with or without configured return values) as a substitute for real Redis or result backends
- `patch("StateFileSource.__call__")` to hide token resolution
- `patch("httpx.Client")` to avoid real HTTP (use `server_factory` instead)
- `monkeypatch.setattr` on module-level functions to skip service interactions (e.g., `wait_for_server_ready`, `_acquire_admin_jwt`, `_is_url_healthy` when a real server could be used instead)

**Allowed:**
- `@patch("os.kill")` — prevents real process termination
- `monkeypatch.setattr("webbrowser.open", ...)` — prevents browser opening
- `monkeypatch.setattr("time.sleep", ...)` — prevents test delays
- `MockSioServer` in route tests — lightweight fake for broadcast verification
- `monkeypatch.setenv` / `monkeypatch.delenv` — standard environment isolation

**Questionable (evaluate during implementation):**
- `patch.dict("sys.modules", {"PIL": None})` in `test_gif.py` — tests what happens when PIL is missing, but PIL is a required dependency. Consider removing the test entirely.

**Last resort (requires justifying comment in the test):**
- `monkeypatch.setattr("zndraw.cli.StateFile", lambda: StateFile(directory=tmp_path))` — only if the code under test cannot accept a `StateFile` parameter. StateFile is critical infrastructure and should be tested for real. Redirecting to `tmp_path` is acceptable for filesystem isolation but question whether the test even needs it.

**Guiding principle:** If you're patching something to avoid testing it, you're writing the wrong kind of test.

#### Fixture Rules

- All shared fixtures live in `conftest.py` — no per-file session/client/Redis definitions.
- `helpers.py` functions (`create_test_user_in_db`, `create_test_room`, `auth_header`, `create_test_token`) are the ONLY way to create test data. The per-file `_create_user` copies do both user creation AND token creation — the consolidated path uses `create_test_user_in_db()` (returns `(user, token)` tuple) which already handles both.
- Fixture scoping: function-scoped for data, session-scoped only for truly stateless infrastructure.

#### The Iron Law: When Writing New Tests You MUST Follow This Pattern

1. **`/brainstorming` first** — plan what to test, which test tier (route / E2E / unit), which fixtures to use. Do not write tests blindly.
2. **`/test-driven-development`** — RED/GREEN/REFACTOR cycle. Write the failing test first, watch it fail, then implement.
3. **Review against these guidelines** — does the test use real services? Does it avoid banned mock patterns? Does it use shared fixtures?

This applies to ALL test writing — manual, AI-assisted, or agent-driven. No exceptions.

#### Required E2E Coverage

Critical user flows MUST have E2E tests against real servers:
- CLI upload flow (`_acquire_admin_jwt` → store JWT → ZnDraw client → write frames)
- ZnDraw Python client connect → write → read → disconnect → reconnect → verify
- Guest auth → room operations
- StateFileSource with real running server

### Phase 2: Shared Fixture Consolidation

Replace per-file fixture sets with shared fixtures in `tests/zndraw/conftest.py`:

```text
conftest.py fixtures:
  session         — in-memory SQLite (already exists)
  redis_client    — real Redis, flushed per test (already exists)
  client          — ASGITransport + real session + real Redis + MockSioServer
  test_user       — creates user via helpers.create_test_user_in_db()
  test_room       — creates room via helpers.create_test_room()
```

Key changes to the existing `client` fixture:
- Replace `app.dependency_overrides[get_redis] = lambda: (None)` with real `redis_client`
- Add `FrameStorage` backed by real Redis
- Add `ResultBackend` dependency — extract `InMemoryResultBackend` from `zndraw_joblib/conftest.py` to a shared location (e.g., `tests/shared_helpers.py` or `tests/conftest.py`) so it can be used by `tests/zndraw/` without cross-package conftest imports
- Keep `MockSioServer` for route tests (legitimate fake for broadcast verification)

### Phase 3: Route Test File Conversion

All route test files delete their per-file fixtures and use shared ones:

**Files:** `test_routes_bookmarks.py`, `test_routes_edit_lock.py`, `test_routes_figures.py`, `test_routes_frame_selection.py`, `test_routes_frames.py`, `test_routes_geometries.py`, `test_routes_presets.py`, `test_routes_selection_groups.py`, `test_routes_step.py`, `test_screenshots.py`, `test_chat.py`, `test_isosurface.py`, `test_progress.py`, `test_default_camera.py`, `test_frames_provider_dispatch.py`, `test_trajectory.py`

**Per file:**
- Delete per-file session fixture (`bm_session`, `el_session`, etc.)
- Delete per-file client fixture (`bm_client`, `el_client`, etc.)
- Delete per-file `_create_user()`, `_create_room()`, `_auth()` — use `helpers.py`
- Update test function signatures to use shared fixture names

### Phase 4: Mock Cleanup

| Target | Files | Replacement |
|--------|-------|-------------|
| `AsyncMock()` for result backend | `test_isosurface.py`, `test_routes_frames.py` | `InMemoryResultBackend` (extracted to shared location in Phase 2) |
| `patch("StateFileSource.__call__", return_value={})` | `test_client_settings.py`, `test_client_api.py` | Real `StateFile(directory=tmp_path)` via constructor injection |
| `patch("zndraw.auth_utils.httpx.Client")` | `test_resolve_token.py` | Test against real server with `server_factory` |
| `patch("zndraw.cli_agent.auth.httpx.Client")` | `test_cli_auth.py`, `test_cli_agent/test_auth.py` | Test against real server |
| `patch("_is_url_healthy")` + `_is_pid_alive` | `test_state_file_source.py` (~19 occurrences) | Split: pure logic unit tests (no patches) + integration tests with `server_factory` |

**`test_cli.py` evaluation:** This file has ~24 monkeypatch calls that mock `uvicorn.Server.run`, `wait_for_server_ready`, `upload_file`, `_acquire_admin_jwt`, etc. These are evaluated case-by-case during implementation:
- Patches for call-ordering tests (browser-before-upload) may stay if they test orchestration logic
- Patches that skip real service interactions should be converted to `server_factory` integration tests where feasible

### Phase 5: Missing E2E Tests

New test files using `server_factory`:

1. **Client persistence E2E** — connect → write frames → disconnect → reconnect → verify data
2. **Guest auth E2E** — `POST /v1/auth/guest` → JWT → create room → write frames → read back
3. **StateFileSource integration** — start server → StateFileSource discovers it → resolves token → client connects

(CLI upload E2E already added in PR #896)

### Phase 6: Parametrize Opportunities

Convert obvious candidates:
- `test_routes_geometries.py`: 404-assertion tests (room-not-found + geometry-not-found) → parametrized
- `test_routes_edit_lock.py`: 3 similar 403 authorization tests → parametrized
- `test_auth_endpoints.py`: 8 login/register tests → 2 parametrized tests

## Out of Scope

- Class-to-function test conversion (cosmetic)
- `zndraw_auth/conftest.py` refactoring (separate package)
- `zndraw_joblib/conftest.py` refactoring (well-designed fakes, except extracting `InMemoryResultBackend` to shared location)
- Subprocess-based worker tests (tracked in #898)

## Success Criteria

- No Redis dependency override uses `AsyncMock()`, `lambda: None`, or `lambda: (None)` anywhere in `tests/`
- No test file outside `conftest.py` defines any session fixture (including prefixed variants like `bm_session`, `el_session`) or any client fixture
- All 4 critical E2E gaps covered with real-server tests
- `tests/README.md` exists and documents the rules
- Every remaining `monkeypatch.setattr` / `@patch` call has a short inline comment explaining WHY it is justified
- All existing tests still pass
