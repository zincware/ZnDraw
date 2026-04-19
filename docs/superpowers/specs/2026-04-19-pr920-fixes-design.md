# PR #920 Fixes — Design Spec

**Date:** 2026-04-19
**Branch:** `spec/dockview-ui-redesign`
**PR:** [zincware/ZnDraw#920](https://github.com/zincware/ZnDraw/pull/920)
**Scope:** Resolve every finding from code review + CodeRabbit (3 Critical, 8 Important, 12 Minor) and fix two active bugs that make the filesystem provider unusable on dev and red on CI.

---

## 1. Problem statement

PR #920 bundles two features on one branch:

1. **Dockview UI redesign** — replaces the legacy sidebar + `react-rnd` window manager with `dockview-react`.
2. **Default `@internal` filesystem provider** — seeds a `ProviderRecord` at `@internal:filesystem:FilesystemRead` on startup, configured via `ZNDRAW_SERVER_FILEBROWSER_PATH`, with a new taskiq dispatch path and `InternalProviderExecutor` that POSTs results back via `/v1/joblib/providers/{id}/results`.

The current state:

- **CI (all 3 pytest jobs, Python 3.11/3.12/3.13) is red.** 13 tests fail across `test_internal.py`, `test_provider.py`, and `test_providers_filesystem.py`.
- **The filesystem provider does not work on `uv run zndraw`.** GET requests to `@internal:filesystem:FilesystemRead` hang indefinitely, and the whole server wedges under subsequent load.
- **Code review raised 23 findings** (3 Critical, 8 Important, 12 Minor) including a silent-exception executor path, a stuck inflight-lock path, dead stub components, a `Plotly.purge` memory leak, a `useLeaveRoom` router desync, and many smaller items.

The goal of this spec is a single design that brings the PR to green-CI and end-to-end-working state, with root-cause fixes only (no xfail, no skip, no try/except suppression), delivered as five logically-grouped commits on the existing branch.

---

## 2. Root cause analysis

Investigation under the systematic-debugging skill identified **two independent bugs** behind the symptoms, plus a set of smaller code-review findings:

### Bug A — Shared Redis causes cross-test state pollution (CI-only)

The `server_factory` fixture at `tests/zndraw/conftest.py:215-298` hardcodes `ZNDRAW_SERVER_REDIS_URL="redis://localhost"` with no per-test namespacing. All test servers share:

- The same `ListQueueBroker` default queue name — tasks enqueued by one test server can be claimed by another test's receiver pointing at a torn-down prior server (→ `httpx.ConnectError` on `receiver.py:282`).
- The same `RedisResultBackend` key space — stale `provider-inflight:…` locks from prior tests block new tests from dispatching; stale `provider-result:…` entries satisfy GETs with wrong data.

This is the sole cause of the 13 CI failures. Main (`bb2198a9`) was green on run `24327808136` before this branch because the PR slightly increases lifespan startup time (added `Path.resolve()` + `register_internal_providers`), which is enough to change the timing window and expose the race.

### Bug B — Global `asyncio.Lock` deadlock on long-polling handlers (dev + CI)

`src/zndraw/database.py:88-102` wraps `app.state.session_maker` with a global non-reentrant `asyncio.Lock` (`_apply_sqlite_locking`). This lock is required because in-memory SQLite and NFS-mounted SQLite don't support WAL — it is the only portable way to serialize SQLite writes across deployments. The lock stays.

The problem is how `zndraw-auth`'s `get_session` dependency (`zndraw_auth/db.py:108-124`) interacts with it:

```python
async def get_session(session_maker: ...) -> AsyncIterator[AsyncSession]:
    async with session_maker() as session:
        yield session
```

This is a **yield-based FastAPI dep**: the `async with` context (and therefore the lock) is held for the **full request lifetime**, not just during DB work.

`SessionDep = Annotated[AsyncSession, Depends(get_session)]` is transitively used by `get_worker_token` at `src/zndraw_joblib/dependencies.py:151-173`, which is exposed as `WorkerTokenDep`.

The new handlers introduced by this PR — `read_provider` (`router.py:1155-1244`) and `upload_provider_result` (`:1281-1322`) — take `SessionMakerDep` for explicit short-lived session work AND `WorkerTokenDep` for token minting. This causes **two independent acquisitions of the same non-reentrant lock within one request**:

1. FastAPI resolves `WorkerTokenDep` → `get_session` opens `async with session_maker()` → acquires lock.
2. Handler body runs → `async with session_maker() as session:` at `router.py:1172` → attempts to acquire the same lock → **deadlock**.

Subagent instrumentation confirmed: the coroutine never progresses past `validate_room_id()`; it blocks at the handler's own `async with session_maker()`. The CPU spike (99-166% at server idle) is baseline fakeredis + socket.io + taskiq polling noise, not caused by the deadlock — it exists before any request hits the server.

### Architectural rule (memorized)

> **Long-polling handlers MUST use factory DI (`SessionMakerDep`), never yield-based session DI (`SessionDep`). Transitively: no sub-dependency of a long-polling handler may hold a `SessionDep` either, because FastAPI keeps yielded deps alive for the whole request.**

This rule is saved as a persistent feedback memory for future review.

### Other review findings (13 Criticals / Importants / Minors)

Enumerated in §4 per commit. Highlights:

- Executor silently swallows exceptions (Critical) — `src/zndraw/providers/executor.py`.
- Inflight lock stuck on dispatch failure (Critical) — `src/zndraw_joblib/router.py:1189-1222`.
- `Plotly.purge` never runs on unmount (Important) — `frontend/src/panels/PlotView.tsx:711-719`.
- `useLeaveRoom` uses `pushState` that doesn't notify react-router (Important).
- `ChatPanel` unread-reset fires when panel is hidden (Important).
- Dead-code: `frontend/src/panels/stubs.tsx` is unused post-migration.
- Duplicated filebrowser bootstrap in `broker.py` + `database.py`.
- `ensure_internal_providers` not concurrent-startup-safe.
- Stale `@internal` `ProviderRecord` rows on disable.
- And 12 minor code-quality items.

---

## 3. Design goals

1. **CI goes from 13 red → all green.**
2. **Filesystem provider works end-to-end on `uv run zndraw`** — browser-verified, not just test-verified.
3. **Every backend fix ships with a regression test** that would fail on pre-fix code.
4. **No lazy fixes** — no `xfail`, no `skip`, no silent `try/except`, no symptom suppression.
5. **Architectural rule enforced** — add a static test that fails if any long-polling handler combines `SessionMakerDep` with a `SessionDep`-holding sub-dep.
6. **All 23 review findings resolved** in 5 logically-grouped commits on the existing branch.

---

## 4. Commit plan

Five commits on `spec/dockview-ui-redesign`, each independently green locally.

### Commit 1 — CI unblocker + deadlock fix

**Goal:** green CI and working `uv run zndraw` filesystem provider.

**Config changes (`src/zndraw/config.py`):**

- Replace the stringly-typed `filebrowser_path: str = "."` sentinel with `filebrowser_path: str | None = "."` — `None` means disabled. Add `env_parse_none_str="none"` to `model_config` so the env-var form `ZNDRAW_SERVER_FILEBROWSER_PATH=none` still parses to `None` cleanly.
- Add `task_queue_name: str = "zndraw:tasks"` and `result_backend_key_prefix: str = "zndraw"` — Pydantic-configurable, documented. Defaults match current implicit behavior so released operators see zero change.

**Bootstrap wiring (`src/zndraw/broker.py`, `src/zndraw/database.py`):**

- All `settings.filebrowser_path.lower() != "none"` → `settings.filebrowser_path is not None`.
- `ListQueueBroker(redis_url, queue_name=settings.task_queue_name)` at both call sites.
- `RedisResultBackend` gets a `key_prefix: str = ""` kwarg that prepends to all keys (including the `notify:` channel derivation). Propagated through `CompositeResultBackend` construction in `database.py:393-396`.

**DI fix — root cause of Bug B (`src/zndraw_joblib/dependencies.py`):**

Rewrite `get_worker_token` from yield-based session DI to factory DI:

```python
async def get_worker_token(
    request: Request,
    session_maker: Annotated[
        async_sessionmaker[AsyncSession], Depends(get_session_maker)
    ],
) -> str:
    settings = request.app.state.settings
    auth_settings = request.app.state.auth_settings
    async with session_maker() as session:
        result = await session.exec(
            select(User).where(User.email == settings.internal_worker_email)
        )
        user = result.one_or_none()
        if user is None:
            raise RuntimeError(
                f"Internal worker user '{settings.internal_worker_email}' not found. "
                "Has the database been initialized?"
            )
    # session + lock released here
    strategy = JWTStrategy(
        secret=auth_settings.secret_key.get_secret_value(),
        lifetime_seconds=auth_settings.token_lifetime_seconds,
    )
    return await strategy.write_token(user)
```

The `expire_on_commit=False` setting on the engine (`database.py:258-260`) ensures the detached `user` object retains its attributes for `write_token`.

**Test fixture changes (`tests/zndraw/conftest.py`):**

Per-server isolation in `server_factory._create_server`:

```python
defaults = {
    "ZNDRAW_SERVER_REDIS_URL": "redis://localhost",
    "ZNDRAW_SERVER_DATABASE_URL": "sqlite+aiosqlite://",
    "ZNDRAW_SERVER_HOST": host,
    "ZNDRAW_SERVER_PORT": str(port),
    "ZNDRAW_SERVER_TASK_QUEUE_NAME": f"zndraw:tasks:{port}",
    "ZNDRAW_SERVER_RESULT_BACKEND_KEY_PREFIX": f"zndraw:{port}",
    "ZNDRAW_SERVER_FILEBROWSER_PATH": "none",
}
```

Tests that exercise the default `@internal` filesystem provider opt back in via `server_factory({"ZNDRAW_SERVER_FILEBROWSER_PATH": "."})` — `test_providers_filesystem.py::test_default_internal_filesystem_*` and `test_filebrowser_path_none_disables_default` already use `server_factory` patterns, update them to opt in where needed.

The 4 failing `test_provider.py` tests (`test_register_fs_creates_provider`, `test_provider_not_visible_in_other_room`, `test_provider_disconnect_removes_provider`, `test_guest_can_register_fs_provider`) pass with zero assertion changes once the default provider is disabled in the `server` fixture.

**New tests:**

- `tests/zndraw_joblib/test_deadlock_regression.py::test_read_provider_completes_without_deadlock` — fire one `read_provider` request against an in-memory SQLite server; assert completes within 2s. Fails on pre-fix `get_worker_token`.
- `tests/zndraw_joblib/test_deadlock_regression.py::test_concurrent_reads_no_deadlock` — 3 parallel reads, all must complete. Catches lock-starvation regressions.
- `tests/zndraw_joblib/test_deadlock_regression.py::test_get_worker_token_releases_lock` — directly call `get_worker_token`, then immediately `async with session_maker()` in the same task; must not deadlock.
- `tests/zndraw_joblib/test_di_architecture.py::test_no_long_polling_handler_holds_full_request_session` — static audit: iterate FastAPI routes, inspect the `Depends` graph, assert no route with `SessionMakerDep` transitively depends on `SessionDep`.
- `tests/zndraw/test_config.py::test_task_queue_name_default` and `test_filebrowser_path_none_parses_to_none` — config parsing correctness.

**Acceptance criteria:**

- `uv run pytest tests/zndraw tests/zndraw_joblib -q` — all green locally.
- CI `pytest (3.11/3.12/3.13, ubuntu-latest)` — all green.
- `uv run zndraw` (no env vars) — default behavior unchanged: server boots with `@internal:filesystem:FilesystemRead` registered at cwd, Redis keys prefixed `zndraw:*`, queue named `zndraw:tasks`.
- **End-to-end browser verification:** open a browser at `http://localhost:8000`, click the filesystem activity-bar icon, see files listed within 1-2 seconds. No 30s timeout.

**Risks:**

- `env_parse_none_str="none"` is model-wide — any string field where `"none"` is a legitimate value would now parse to `None`. Current fields are safe (verified). A comment on `Settings` documents this so whoever adds a new field thinks about it.
- `RedisResultBackend.__init__` signature change may ripple into direct-instantiation in existing tests. Grep-verify before the commit lands.

### Commit 2 — Backend correctness

**Goal:** production-readiness of the `@internal` provider path. Defense-in-depth even though Commit 1 fixes the primary deadlock.

**Extract bootstrap helper (new file `src/zndraw/providers/bootstrap.py`):**

```python
def register_filebrowser_providers(
    broker: AsyncBroker,
    *,
    base_url: str,
    settings: Settings,
) -> InternalProviderRegistry | None:
    """Register the default @internal filesystem provider on a broker.

    Returns None when settings.filebrowser_path is None (disabled).
    Callers: src/zndraw/broker.py (external worker), src/zndraw/database.py (in-process).
    """
    if settings.filebrowser_path is None:
        return None
    executor = InternalProviderExecutor(
        base_url=base_url,
        filebrowser_path=str(Path(settings.filebrowser_path).resolve()),
        timeout_seconds=settings.provider_executor_timeout,
    )
    return register_internal_providers(broker, _collect_providers(), executor)
```

Both `broker.py` (module-level, for `taskiq worker zndraw.broker:broker`) and `database.py` lifespan (in-process worker) collapse to one line each. Resolves CodeRabbit's duplicated-bootstrap finding.

**Inflight-lock release on dispatch failure (`src/zndraw_joblib/router.py:1189-1222`):**

Wrap both branches (the `@internal` kiq dispatch and the remote-provider `emit(...)` branch) in `try/except` that calls `result_backend.release_inflight(inflight_key)` before re-raising. Pattern mirrors `submit_task` at `router.py:642-660`.

**Executor error POST-back (`src/zndraw/providers/executor.py`):**

Restructure `_run` to catch exceptions from handler resolution, param parsing, provider instantiation, `Provider.read`, and serialization; POST an error body with `X-Result-Status: error` header to the same results endpoint. The server-side `upload_provider_result` learns to recognize this header and write the error payload so the long-poll resolves as a 4xx JSON to the client (`{"error": "...", "type": "..."}`) instead of a 504 timeout.

```python
def _run() -> None:
    try:
        handler = self._resolve_handler(provider_cls)
        params = json.loads(params_json) if params_json else {}
        instance = provider_cls(**params)
        result = instance.read(handler)
        content = (
            json.dumps(result).encode()
            if provider_cls.content_type == "application/json"
            else result
        )
        headers = {
            "Authorization": f"Bearer {token}",
            "X-Request-Hash": request_id,
        }
    except Exception as err:
        content = json.dumps({"error": str(err), "type": type(err).__name__}).encode()
        headers = {
            "Authorization": f"Bearer {token}",
            "X-Request-Hash": request_id,
            "X-Result-Status": "error",
        }
    with httpx.Client(timeout=self.timeout_seconds, transport=transport) as client:
        resp = client.post(
            f"{base_url}/v1/joblib/providers/{provider_id}/results",
            content=content,
            headers=headers,
        )
        resp.raise_for_status()
```

**`ensure_internal_providers` concurrent-startup safe (`src/zndraw_joblib/registry.py:162-216`):**

Wrap `session.commit()` in `try/except IntegrityError → re-query + update`. Resolves CodeRabbit's unique-constraint race finding. Matches idempotent-seed pattern elsewhere.

**Stale row cleanup on disable (`src/zndraw/database.py:219`):**

In the `filebrowser_path is None` branch, `DELETE FROM providerrecord WHERE room_id = '@internal'`. Ensures toggling the feature off actually disables it on existing DBs.

**Executor timeout from settings:**

Add `Settings.provider_executor_timeout: float = 30.0`. Threaded through `InternalProviderExecutor.timeout_seconds`. Satisfies CodeRabbit's hardcoded-timeout finding.

**New tests:**

- `tests/zndraw_joblib/test_provider_dispatch.py::test_read_provider_releases_inflight_on_kiq_failure` — monkeypatch `tasks[…].kiq` to raise; assert `inflight_key` absent in Redis afterward.
- `tests/zndraw/test_providers_executor.py::test_internal_filesystem_provider_surfaces_error` — pass an invalid path; assert 4xx with body `{"error": "...", "type": "..."}`.
- `tests/zndraw_joblib/test_registry.py::test_ensure_internal_providers_concurrent_startup_safe` — run two `ensure_internal_providers` calls concurrently against the same session_maker; neither raises, exactly one `ProviderRecord` exists.
- `tests/zndraw/test_providers_filesystem.py::test_filebrowser_path_none_removes_stale_rows` — boot with path set, populate row, shut down, re-boot with `path=None`; assert no `@internal` rows remain.
- `tests/zndraw/test_providers_executor.py::test_executor_timeout_from_settings` — construct executor with `provider_executor_timeout=5`; assert httpx client timeout is 5.

**Acceptance criteria:**

- Commit 1's CI-green state preserved.
- `grep -r "filebrowser_path.lower" src/` → zero results.
- `grep -r "InternalProviderExecutor(" src/` → exactly one construction site (`providers/bootstrap.py`).
- New tests green.

**Out of Commit 2:**

- `Provider.encode_result` hook (CodeRabbit's content-type suggestion) — premature abstraction before a second internal-provider category exists. Stays out.

### Commit 3 — Frontend panel hooks

**Goal:** six behavior-affecting frontend fixes from the review. No dead-code or style polish (reserved for Commit 4).

**Fix 3.1 — `useLeaveRoom` uses react-router navigate (`frontend/src/hooks/useLeaveRoom.ts`):**

Replace `window.history.pushState({}, "", "/")` with `navigate("/")` from `useNavigate()`. Update callers `RoomsPanel.tsx`, `FilesystemPanel.tsx`, `ViewerView.tsx` to thread navigation in. Playwright test: close viewer, assert URL and panel content both reset.

**Fix 3.2 — `useLeaveRoom` api-getter pattern:**

Change `useLeaveRoom` signature to accept `api: DockviewApi | (() => DockviewApi | null)`. Resolve on-call. Update all three call sites to pass `() => getDockviewApi()` instead of `getDockviewApi()`. Vitest test: `getDockviewApi()` returns null initially, populated after onReady; `leaveRoom` invocation sees the populated api.

**Fix 3.3 — `ChatPanel` visibility-gated unread reset (`frontend/src/panels/ChatPanel.tsx:212-217`):**

Use `IDockviewPanelProps.api.isVisible` + subscribe to `api.onDidVisibilityChange`. Reset only when visible. Keep `scrollToBottom` unconditional. Playwright test: hide chat tab, trigger messages, switch back, assert unread count accurate.

**Fix 3.4 — `PlotView` Plotly purge leak (`frontend/src/panels/PlotView.tsx:711-719`):**

Move `plotContainer.current` read inside the cleanup function instead of capturing at mount. Playwright test: open/close plot 10 times, assert `document.querySelectorAll('.plotly').length <= 1`.

**Fix 3.5 — `DockviewLayout` disposable cleanup (`frontend/src/panels/DockviewLayout.tsx:66-84`):**

Collect `onUnhandledDragOverEvent`, `onDidRemovePanel`, `onDidAddPanel` disposables in a ref; dispose + null `sharedApi` in a cleanup `useEffect`. Vitest test: mount, unmount, re-mount; assert disposables invoked.

**Fix 3.6 — `RoomsPanel.handleFiles` drop-to-create (`frontend/src/panels/RoomsPanel.tsx:49-67`):**

Call `await leaveRoom({ skipConfirm: true })` before `navigate`. On upload failure, best-effort `deleteRoom({ room_id: newRoomId })`. Playwright test: in room A with plot open, drop files; assert URL changes, old plot gone, new room loads cleanly.

**Fix 3.7 — `PlotsBrowserPanel` useOpenPlotKeys retries until API ready (`frontend/src/panels/PlotsBrowserPanel.tsx:42`):**

Empty-dep `useEffect` skips registration if `getDockviewApi()` returns `null` at mount (happens on direct URL navigation with `?panel=plots-browser` before `DockviewLayout.onReady` fires). Gate the effect on API readiness using the getter pattern from Fix 3.2 — poll `getDockviewApi()` with a short `setInterval` until it resolves, then register listeners and clear the interval. Vitest test: mount with `getDockviewApi()` returning null, then populated 50ms later; assert listeners registered after.

**Acceptance criteria:**

- `bun --bun tsc --noEmit` → zero errors.
- `bun test` → green.
- Five new Playwright + three new Vitest tests → green.
- Browser smoke: golden path from §6 executes cleanly.

### Commit 4 — Polish + dead code

Cleanup only, no behavior changes.

| Item | File:line | Change |
|---|---|---|
| 4.1 | `frontend/src/panels/stubs.tsx` | **Delete entire file.** Zero imports remain (verified). |
| 4.2 | `frontend/src/panels/roomRowMenu.tsx:97-104` | Flip icons: `room.is_default ? <StarIcon /> : <StarBorderIcon />`. |
| 4.3 | `frontend/src/utils/errors.ts` (new) | `extractDetail(err, fallback): string` — coerce FastAPI validation arrays to strings. Replace inline use in `roomsHeaderActions.tsx:22-26` and `RoomsPanel.handleFiles`. |
| 4.4 | `frontend/src/panels/ChatPanel.tsx:522-527` | `onKeyPress` → `onKeyDown`. |
| 4.5 | `frontend/src/panels/SidebarZone.tsx:14-17, 158-160` | Replace `rgba(25, 118, 210, …)` with `alpha(theme.palette.primary.main, …)`. |
| 4.6 | `frontend/src/hooks/socketHandlers/figureHandlers.ts:52-56` | Drop redundant `getPanel` check; `openPlotTab` already handles exists-case. |
| 4.7 | `frontend/src/panels/groupActions.tsx:12-21` | Remove `onDidLayoutChange` subscription. |
| 4.8 | `frontend/src/stores/slices/activityBarSlice.ts:82-153` | Extract `findSourceBar(state, id)` helper. |
| 4.9 | `frontend/src/panels/ViewerView.tsx:13-14` | Inline `leaveRoomRef.current = leaveRoom` → `useEffect(() => { leaveRoomRef.current = leaveRoom; }, [leaveRoom])`. |
| 4.10 | `frontend/src/panels/DockviewLayout.tsx:38-57` | `onDidDrop` handler: replace manual `addPanel` + `getPanel` check with `plotViewFactory.openPlotTab(api, key, { referenceGroup: event.group, direction: "within" })`. Centralizes "exists → activate, else add" logic. |
| 4.11 | `frontend/src/panels/RoomsPanel.tsx:214-238` | Extract shared `useRoomRowActions(room)` hook covering `onToggleTemplate` + `onToggleLock`. Consume from both `RoomsPanel.tsx` row icons and `roomRowMenu.tsx` menu items. Removes race-condition surface where double-tap fires concurrent API calls. |
| 4.12 | `frontend/src/pages/landingPage.tsx:459-477` | Reset-Layout teardown ordering: close dockview panels first, then call `useAppStore.getState().resetLayout()`, then re-add viewer. Fold the sequence into a single action (`resetLayoutAndTeardown`) consumed by the menu item and `useLeaveRoom` effect. |
| 4.13 | `frontend/src/pages/landingPage.tsx:235-250` | `?panel=` effect: subscribe to the active-bar state via `useAppStore(selector)` instead of imperative `getState()`; remove the redundant `def.default.bar === "editor"` clause since `def.kind !== "tool"` already guards it. |
| 4.14 | `frontend/e2e/dockview-layout.spec.ts:37-46` | "only one panel per bar — switching icons swaps the panel" test: add explicit assertions on `activity-icon-selections` and `activity-icon-modifiers` active-state (e.g. `aria-pressed`) before/after the click, so a regression in the one-panel-per-bar invariant would actually fail. Adding `aria-pressed` / `data-active` to `ActivityBar.tsx` icon buttons may be required to support this — include that. |
| 4.15 | `tests/zndraw_joblib/test_providers.py:145-230, 423-461, 640-704` | Extract shared `_seed_internal_provider` async fixture (User + Worker + `ProviderRecord(room_id="@internal", category="filesystem", name="FilesystemRead")`). Replace the four duplicated `async def seed()` blocks. |

**Acceptance:** `bun --bun tsc --noEmit` clean; `bun test` green; `bun --bun biome check` clean on edited files. Manual browser smoke — no visual regressions.

### Commit 5 — Docs

| Item | File | Change |
|---|---|---|
| 5.1 | `docs/superpowers/specs/2026-04-17-internal-providers-design.md:100` | Table entry: `@internal:filesystem:local` → `@internal:filesystem:FilesystemRead`. |
| 5.2 | `docs/superpowers/plans/2026-04-16-dockview-ui-redesign.md` | Replace `/Users/fzills/tools/zndraw-fastapi` paths with `$(git rev-parse --show-toplevel)`. Add language tags to fences at L193, L217, L238, L421. |
| 5.3 | `frontend/src/panels/dockview-mui.css:6-7` | Update header comment: `<MuiCssVars />` emits `--mui-palette-*`, not `cssVariables: true`. |
| 5.4 | `docs/superpowers/plans/2026-04-17-dockview-ui-fixes.md:746-762` | Merge duplicate `borderStyle` entries into one conditional. |
| 5.5 | `docs/superpowers/plans/2026-04-17-internal-providers.md:1373-1387` | Update the Step 8.3 snippet to match the actual code in `database.py:228`, which already captures `internal_user_id = internal_user.id` before the session closes. CodeRabbit flagged the snippet as having a detached-instance bug; the plan doc is stale, the code is correct. |

**Acceptance:** `uvx prek run --files <edited>` green; RTD build green; `uv run sphinx-build` no new warnings.

---

## 5. Testing strategy

| Commit | New tests | What they catch |
|---|---|---|
| 1 | `test_read_provider_completes_without_deadlock`, `test_concurrent_reads_no_deadlock`, `test_get_worker_token_releases_lock`, `test_no_long_polling_handler_holds_full_request_session`, config-parsing tests | Lock deadlock; factory-DI discipline; fixture isolation; config correctness |
| 2 | `test_read_provider_releases_inflight_on_kiq_failure`, `test_internal_filesystem_provider_surfaces_error`, `test_ensure_internal_providers_concurrent_startup_safe`, `test_filebrowser_path_none_removes_stale_rows`, `test_executor_timeout_from_settings` | Error-path correctness; bootstrap concurrency; feature-disable cleanup |
| 3 | Playwright: `viewer_close_leaves_room_cleanly`, `chat_unread_badge_while_hidden`, `plotview_no_plotly_leak_on_repeated_open_close`, `rooms_panel_drop_leaves_current_room`; Vitest: `useLeaveRoom_api_getter_stale_ref_safe`, `DockviewLayout_disposables_cleanup`, `PlotsBrowserPanel_registers_listeners_when_api_ready` | Panel behavior regressions |
| 4 | None (polish only) | Existing coverage must stay green |
| 5 | None (docs only) | — |

**Constraints (enforced by self-audit before each commit):**

- Zero new `@pytest.mark.xfail` or `@pytest.mark.skip` markers.
- Zero existing tests deleted.
- Zero `# type: ignore` or `noqa` additions unless the reason is written in a comment.
- No `try/except Exception` broader than the specific error the code handles; re-raise or structured error path required.

---

## 6. End-to-end verification (before claiming PR ready)

After all 5 commits land on the branch:

1. `uv sync --reinstall-package zndraw && uv run zndraw` with no env vars.
2. Open `http://localhost:8000` in a browser.
3. **Golden path:**
   - Create a new room via Rooms panel → drop an `.xyz` file → atoms render in viewer.
   - Click filesystem activity-bar icon → file listing of cwd renders within 2 seconds (validates Bug B fix).
   - Click a file → opens in appropriate panel.
   - Open a plot → close → repeat 5 times. Inspect devtools: `document.querySelectorAll('.plotly').length` stays ≤ 1 (validates Fix 3.4).
   - Switch chat panel to an inactive tab → send a message from another session → switch back → unread badge reflected correct count (validates Fix 3.3).
   - Close the viewer → lands on Rooms panel cleanly, URL resets to `/` (validates Fix 3.1).
   - Drop files onto Rooms panel while in an existing room → new room loads, prior room's plots closed (validates Fix 3.6).
4. Browser devtools console — zero errors, zero unexpected warnings.
5. Server logs — only `INFO` startup messages + expected request access logs, no tracebacks.
6. `uv run pytest tests/ -q` — all green.
7. `cd frontend && bun test && bun --bun tsc --noEmit` — all green.
8. `uvx prek run --all-files` — clean.
9. Push branch; `gh run list --limit 1 --workflow ci` — green across 3.11/3.12/3.13.

Only then is the PR marked ready for review.

---

## 7. Out of scope

Explicit punts (not the right scope for this PR):

- **`Provider.encode_result` hook** (CodeRabbit suggestion on `content_type` discipline) — premature abstraction before a second `@internal` provider category exists.
- **Broader SQLite-lock refactor** — the `asyncio.Lock` stays because in-memory SQLite and NFS SQLite don't support WAL. Fix is scoped to the specific DI offender.
- **Pre-existing diagnostic issues** — e.g. pyright warnings on `router.py:194`, `:515`, `:555`, `:644` are pre-existing on `main` and unrelated to this PR's scope.
- **Dropping `fakeredis` from dev mode** — considered during investigation; not needed now that the deadlock is fixed.
- **Splitting Dockview redesign from `@internal` providers** — the PR is already reviewed together; splitting post-hoc is more churn than value.

---

## 8. Rollout

- Branch stays on `spec/dockview-ui-redesign`, new commits pushed to PR #920.
- No force-push; commits append.
- After CI goes green and end-to-end verification passes, re-request review on the PR.
- Spec + plan documents committed alongside the code changes for provenance.

---

## References

- PR: https://github.com/zincware/ZnDraw/pull/920
- Failing CI run: https://github.com/zincware/ZnDraw/actions/runs/24582633191
- Last green main run: https://github.com/zincware/ZnDraw/actions/runs/24327808136
- CodeRabbit review: attached to PR #920 — 17 actionable + 19 nitpicks
- Internal providers design: `docs/superpowers/specs/2026-04-17-internal-providers-design.md`
- Internal providers plan: `docs/superpowers/plans/2026-04-17-internal-providers.md`
