# Testing Patterns

**Analysis Date:** 2026-03-05

## Backend Test Framework

**Runner:**
- pytest >= 9.0.2 + pytest-asyncio >= 1.3.0
- Config: `pyproject.toml` (`[dependency-groups] dev`)
- Coverage: pytest-cov >= 7.0.0

**Run Commands:**
```bash
uv run pytest tests/              # Run all tests
uv run pytest tests/ -x           # Stop on first failure
uv run pytest tests/ -k "test_bookmarks"  # Run matching tests
uv run pytest tests/ --cov        # With coverage
```

**Important:** Tests can run up to 15 minutes. Be patient (per CLAUDE.md).

## Test File Organization

**Location:** Tests are in a separate `tests/` directory (not co-located):
```
tests/
  conftest.py                    # Shared fixtures (sessions, clients, server factories)
  test_routes_bookmarks.py       # REST API endpoint tests
  test_routes_frames.py
  test_routes_geometries.py
  test_routes_edit_lock.py
  test_routes_figures.py
  test_routes_presets.py
  test_routes_selection_groups.py
  test_routes_frame_selection.py
  test_socketio.py               # Socket.IO integration tests
  test_socketio_scaling.py
  test_redis.py                  # Redis-specific tests
  test_schemas.py                # Schema validation tests
  test_schemas_frames.py
  test_exceptions.py
  test_config.py
  test_auth.py
  test_admin.py
  test_sessions.py
  test_storage_asebytes.py
  test_storage_router.py
  test_analysis.py
  test_connectivity.py
  test_constraints.py
  test_chat.py
  test_screenshots.py
  test_progress.py
  test_trajectory.py
  test_lifespan.py
  test_client_api.py             # Python client library tests
  test_client_source.py
  test_client_token_discovery.py
  test_multi_client.py
  test_cli_agent/                # CLI agent tests (subdirectory)
    conftest.py                  # CLI-specific fixtures
    __init__.py
    test_frames.py
    test_geometries.py
    test_rooms.py
    test_bookmarks.py
    test_auth.py
    test_chat.py
    test_extensions.py
    test_selection.py
    test_step.py
    test_bug_fixes.py
    test_mount.py
    test_room_envvar.py
```

**Naming:**
- `test_<domain>.py` for unit/integration tests
- `test_routes_<resource>.py` for REST API endpoint tests
- `test_cli_agent/test_<command>.py` for CLI tests

## Test Structure

**All tests are functions, not class methods** (per CLAUDE.md rule):

```python
@pytest.mark.asyncio
async def test_list_bookmarks_returns_empty_initially(
    bm_client: AsyncClient,
    bm_session: AsyncSession,
) -> None:
    """Test GET returns empty bookmarks for new room."""
    user, token = await _create_user(bm_session)
    room = await _create_room(bm_session, user)

    response = await bm_client.get(
        f"/v1/rooms/{room.id}/bookmarks",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.json()["items"] == {}
```

**Section organization within test files** using comment separators:
```python
# =============================================================================
# List Bookmarks Tests
# =============================================================================

# =============================================================================
# Set Bookmark Tests
# =============================================================================

# =============================================================================
# Authentication Tests
# =============================================================================
```

**Test naming convention:**
- `test_<action>_<expected_behavior>`: `test_list_bookmarks_returns_empty_initially`
- `test_<action>_<condition>`: `test_set_bookmark_rejects_empty_label`
- `test_<action>_<error_case>`: `test_delete_nonexistent_bookmark_returns_404`

## Fixture Architecture

### Three Tiers of Test Fixtures

**Tier 1: Pure async / unit tests** (fastest, in-memory DB):
- Use `session` fixture (in-memory SQLite via `AsyncSession`)
- Use `client` fixture (httpx `AsyncClient` with ASGI transport)
- Dependency overrides on `app.dependency_overrides`
- Redis mocked out (`get_redis` returns `None` or `AsyncMock`)
- Socket.IO mocked with `MockSioServer` class

```python
# conftest.py
@pytest_asyncio.fixture(name="client")
async def client_fixture(session: AsyncSession) -> AsyncIterator[AsyncClient]:
    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    mock_sio = MockSioServer()
    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_redis] = lambda: None
    app.dependency_overrides[get_tsio] = get_sio_override

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client
    app.dependency_overrides.clear()
```

**Tier 2: Redis integration tests** (requires local Redis):
- Use `redis_client` fixture (connects to `redis://localhost`)
- Flushes DB before/after each test
- Used for edit lock tests, presence tests, etc.

```python
@pytest_asyncio.fixture(name="redis_client")
async def redis_client_fixture():
    from redis.asyncio import Redis
    redis: Redis = Redis.from_url("redis://localhost", decode_responses=True)
    await redis.flushdb()
    yield redis
    await redis.flushdb()
    await redis.aclose()
```

**Tier 3: Full server integration tests** (real uvicorn server):
- Use `server_factory` fixture to spawn real servers
- `server` fixture: basic server
- `server_auth` fixture: server with auth enabled
- Tests connect via httpx/Socket.IO clients to actual HTTP endpoints
- Server spawned in daemon thread, waits for `/v1/health` endpoint

```python
@pytest.fixture(name="server")
def server_fixture(server_factory: ServerFactory) -> str:
    instance = server_factory({})
    return instance.url  # "http://127.0.0.1:<port>"
```

### Environment Configuration for Tests

**Autouse fixture** `test_settings` sets environment variables:
```python
@pytest.fixture(autouse=True)
def test_settings() -> Generator[None, None, None]:
    os.environ["ZNDRAW_DATABASE_URL"] = "sqlite+aiosqlite://"
    os.environ.pop("ZNDRAW_REDIS_URL", None)
    os.environ.pop("ZNDRAW_HOST", None)
    os.environ.pop("ZNDRAW_PORT", None)
    return
```

### Test Helper Functions (conftest.py)

Located in `tests/conftest.py`:
- `create_test_user_model()` — creates User with hashed password
- `create_test_token()` — generates JWT for test user
- `create_test_user_in_db()` — creates user + returns (user, token)
- `create_test_room()` — creates room with owner membership
- `auth_header()` — returns `{"Authorization": "Bearer <token>"}`
- `make_raw_frame()` — converts dict to RawFrame format
- `decode_msgpack_response()` — decodes msgpack HTTP response
- `get_free_port()` — finds available TCP port

### Per-Test-File Fixtures

Route test files often define their own local fixtures to avoid coupling:

```python
# tests/test_routes_bookmarks.py
@pytest_asyncio.fixture(name="bm_session")
async def bm_session_fixture() -> AsyncIterator[AsyncSession]:
    # Own session with own engine

@pytest_asyncio.fixture(name="bm_client")
async def bm_client_fixture(bm_session, mock_sio) -> AsyncIterator[AsyncClient]:
    # Own client with own overrides
```

This keeps tests independent — each test file can have its own DB session.

## Mocking

**Socket.IO Mocking:**
Two approaches depending on test tier:

1. **MockSioServer** (conftest.py) — records emitted events:
```python
class MockSioServer:
    def __init__(self):
        self.emitted: list[dict[str, Any]] = []

    async def emit(self, event_or_model, data=None, *, room=None, ...):
        # Records all emit calls for assertion
```

2. **MagicMock + AsyncMock** — for route-specific tests:
```python
@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MagicMock:
    sio_mock = MagicMock()
    sio_mock.emit = AsyncMock()
    return sio_mock
```

**Redis Mocking:**
- For Tier 1 tests: `AsyncMock()` with specific return values
```python
mock_redis = AsyncMock()
mock_redis.get = AsyncMock(return_value=None)  # No edit lock
```

**Dependency Override Pattern:**
```python
app.dependency_overrides[get_session] = get_session_override
app.dependency_overrides[get_redis] = lambda: mock_redis
app.dependency_overrides[get_tsio] = get_sio_override
```
Always `app.dependency_overrides.clear()` in teardown.

## Assertion Patterns

**Status code + response body:**
```python
response = await client.get(f"/v1/rooms/{room.id}/bookmarks", headers=auth_header(token))
assert response.status_code == 200
assert response.json()["items"] == {}
```

**Validate with Pydantic models:**
```python
StatusResponse.model_validate(response.json())
```

**Problem detail type URI assertion:**
```python
assert response.status_code == 404
assert "bookmark-not-found" in response.json()["type"]
# Or exact match:
assert response.json()["type"] == BookmarkNotFound.type_uri()
```

**Socket.IO broadcast assertions:**
```python
mock_sio.emit.assert_called()
call_args = mock_sio.emit.call_args
model = call_args[0][0]
assert isinstance(model, BookmarksInvalidate)
assert call_args[1]["room"] == f"room:{room.id}"
```

**Database state verification:**
```python
row = await session.get(RoomBookmark, (room.id, 3))
assert row is not None
assert row.label == "New Bookmark"
```

## CLI Agent Tests

Located in `tests/test_cli_agent/`. Use Typer's `CliRunner`:

```python
def test_frames_count(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["frames", "count", "--room", test_room]
    )
    resp = StepResponse.model_validate(data)
    assert resp.total_frames >= 1
```

**Helper function** `invoke_cli()` in `tests/test_cli_agent/conftest.py`:
- Invokes CLI command with `--url` and `--token` injected
- Asserts exit code 0
- Returns parsed JSON from stdout

**Test room fixture** creates a room via CLI before each test:
```python
@pytest.fixture(name="test_room")
def test_room_fixture(server_url, cli_runner, auth_token) -> str:
    result = cli_runner.invoke(app, ["rooms", "create", "--url", server_url, "--token", auth_token])
    return json.loads(result.stdout)["room_id"]
```

## Socket.IO Integration Tests

Located in `tests/test_socketio.py`. Use real server + Socket.IO client:

```python
@pytest.mark.asyncio
async def test_socketio_get_user_authenticated(server_auth: str) -> None:
    async with httpx.AsyncClient() as client:
        # Register user via REST
        register_response = await client.post(f"{server_auth}/v1/auth/register", ...)
        # Login to get token
        login_response = await client.post(f"{server_auth}/v1/auth/jwt/login", ...)
        token = login_response.json()["access_token"]

        # Connect Socket.IO with token
        sio_client = socketio.AsyncClient()
        await sio_client.connect(server_auth, auth={"token": token})

        # Call typed event via wrapper
        tsio_client = wrap(sio_client)
        response = await tsio_client.call(UserGet(), response_model=UserGetResponse)
        assert response.email == "testuser@example.com"

        await sio_client.disconnect()
```

## Test Types Summary

**Unit Tests:**
- Schema validation: `test_schemas.py`, `test_schemas_frames.py`
- Config: `test_config.py`
- Exception types: `test_exceptions.py`
- Pure logic: `test_connectivity.py`, `test_analysis.py`

**Integration Tests (in-memory DB):**
- REST API endpoints: `test_routes_*.py`
- Storage backends: `test_storage_asebytes.py`, `test_storage_router.py`

**Integration Tests (Redis required):**
- Redis operations: `test_redis.py`
- Edit lock flows: `test_edit_lock_integration.py`
- Session management: `test_sessions.py`

**Integration Tests (full server):**
- Socket.IO: `test_socketio.py`, `test_socketio_scaling.py`
- Multi-client: `test_multi_client.py`
- CLI agent: `test_cli_agent/`
- Client library: `test_client_api.py`

## Rules from CLAUDE.md

- **Never use `@pytest.mark.xfail`** — all tests must pass
- Use `pytest.mark.parametrize` where applicable
- Each test must be a function, not a class method
- Tests should be specific and test one thing
- Avoid complex setups
- Redis is available for testing — do not use `isinstance(obj, bytes)` guards

---

## Frontend Tests

### E2E Test Framework

**Runner:**
- Playwright (`@playwright/test` >= 1.58.0)
- Config: `frontend/playwright.config.ts`

**Run Commands:**
```bash
cd frontend && bun run playwright test        # Run all E2E tests
cd frontend && bun run playwright test --ui    # Interactive mode
```

**Configuration:**
- Test directory: `frontend/e2e/`
- Timeout: 90 seconds per test
- Expect timeout: 15 seconds
- Fully parallel with 4 workers
- Browser: Chromium (headless: false)
- Base URL: `ZNDRAW_URL` env var or `http://localhost:8000`
- Screenshots: only on failure
- Trace: retain on failure
- Viewport: 1280x720

### E2E Test File Organization

```
frontend/e2e/
  helpers.ts                    # Shared utilities (CLI, PY, waitForScene)
  registration.spec.ts          # Extension/filesystem registration tests
  chat-features.spec.ts
  socket-sync.spec.ts
  ui-panels-chat.spec.ts
  constraint-visualization.spec.ts
  geometry-drawing.spec.ts
  extensions-analysis.spec.ts
  frame-invalidation.spec.ts
  editing.spec.ts
  camera-session.spec.ts
  docs-screenshots.spec.ts
  frames-navigation.spec.ts
```

### E2E Test Helpers

Located in `frontend/e2e/helpers.ts`:

```typescript
// Run CLI command against test server
export function CLI(cmd: string): string {
    return execSync(`uv run zndraw-cli --url ${BASE_URL} ${cmd}`, { encoding: "utf-8" });
}

// Run Python script
export function PY(code: string): string { ... }

// Spawn background Python process (for extensions that need to stay alive)
export function spawnPY(code: string): ChildProcess { ... }

// Wait for 3D canvas to render
export async function waitForScene(page: Page): Promise<void> {
    await page.waitForTimeout(1000);
    await page.waitForSelector("canvas", { state: "visible", timeout: 15000 });
}

// Create isolated browser context for auth
export async function createAuthContext(browser: Browser): Promise<BrowserContext> { ... }
```

### E2E Test Pattern

```typescript
test.describe("Registration", () => {
    test.describe.configure({ mode: "serial" });

    test("register_job makes extension appear", async ({ page }) => {
        // 1. Setup: create room and data via CLI/Python
        CLI(`rooms create --room-id ${ROOM_EXT}`);
        PY(`...`);

        // 2. Spawn background process (e.g., extension worker)
        const bg = spawnPY(`...`);

        try {
            await waitForBgReady();

            // 3. Navigate and interact
            await page.goto(`${BASE_URL}/rooms/${ROOM_EXT}`);
            await waitForScene(page);

            // 4. Assert UI state
            await expect(option.first()).toBeVisible({ timeout: 15000 });

            // 5. Screenshot for documentation
            await page.screenshot({ path: "e2e/screenshots/..." });
        } finally {
            bg.kill();
        }
    });
});
```

### Unit Tests (Frontend)

**No unit test framework is configured for the frontend.** There is no Jest, Vitest, or similar setup. Frontend testing relies entirely on Playwright E2E tests.

## Coverage

**Requirements:** None enforced (no coverage threshold configured)

**View Coverage:**
```bash
uv run pytest tests/ --cov=zndraw --cov-report=html
```

---

*Testing analysis: 2026-03-05*
