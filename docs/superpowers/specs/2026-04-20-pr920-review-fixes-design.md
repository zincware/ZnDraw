# PR #920 Code-Review Fixes — Design Spec

**Date:** 2026-04-20
**Branch:** `spec/dockview-ui-redesign`
**PR:** [zincware/ZnDraw#920](https://github.com/zincware/ZnDraw/pull/920)
**Scope:** Fixes from the 2026-04-20 code review (2 Opus reviewers — backend + frontend).

---

## 1. Problem

A code review against `main` (base `bb2198a9`, head `56c07969`) surfaced 24 issues across backend and frontend. Four trivial nitpicks have already been applied. This spec covers the remaining 20 issues — 2 Critical, 2 Critical-merged, 10 Important, 6 Minor — plus one review finding that traced to a reviewer misread and will not be fixed.

Severity distribution after the brainstorm:

| | Backend | Frontend |
|---|---|---|
| Critical | 1 | 2 |
| Critical merged into another | 1 → #5 | 1 → #2 |
| Important | 4 | 6 merged → 4 distinct |
| Minor | 3 | 4 merged → 3 distinct |
| Not a defect | 1 | 0 |

Two overarching themes:

1. **Backend**: an `@internal` provider path added route-level deps, content-sniffing, and parallel error shapes that diverge from the codebase's RFC 9457 error convention. These collectively regress remote-provider reads and security invariants.
2. **Frontend**: a `DockviewApi` singleton with cross-tree consumers led to polling workarounds and stale-reference bugs. A redesigned 3-bar drag system shipped with duplicated handler logic and empty-bar UX gaps.

---

## 2. Backend Fixes

### 2.1 Lazy-resolve `worker_token` inside `@internal` branch (Critical)

**Defect:** `src/zndraw_joblib/router.py:1165` added `worker_token: WorkerTokenDep` as a route-level dep on `read_provider`. The token is consumed only inside `if provider.room_id == "@internal":`, but FastAPI resolves it on every call. `get_worker_token` (`src/zndraw_joblib/dependencies.py:158`) raises `RuntimeError` if `app.state.internal_worker_user` is `None`. On a fresh deploy with `init_db_on_startup=False`, the cached user can be `None`; every provider read 500s for the process lifetime.

**Fix:** Drop `worker_token: WorkerTokenDep` from the route signature. Resolve manually inside the `@internal` branch via `request.app.state.internal_worker_user` + the token-minting helper extracted from `get_worker_token`. Remote reads never touch the dep.

**Files:** `src/zndraw_joblib/router.py:1155-1220`, `src/zndraw_joblib/dependencies.py`.

**Test:** regression — a remote provider read succeeds when `app.state.internal_worker_user is None` (e.g., test boots with `init_db_on_startup=False` and skips `ensure_internal_worker`).

---

### 2.2 + 2.5 `ProviderExecutionFailed` ProblemType + preserve `X-Result-Status` header (Critical + Important, merged)

**Defect:**
- `src/zndraw_joblib/router.py:1183-1194, 1255-1265` — cached provider payloads are inspected with `if "error" in parsed and "type" in parsed` and translated to HTTP 400. Any legitimate payload with those two top-level keys (JSON-Schema docs, CloudEvents envelopes) is mis-flagged as an error.
- `src/zndraw/providers/executor.py:72-75` — the error payload is an ad-hoc `{"error": str(err), "type": type(err).__name__}` shape that does not match the codebase's RFC 9457 `ProblemDetail` schema (`src/zndraw_joblib/exceptions.py:17`). No logging on failure.

**Root cause:** `executor.py:79` already sets `X-Result-Status: error` on the upload, but `upload_result` (`router.py:1338`) drops the header and stores only the body. The router then re-derives the signal by sniffing content — and picks a shape incompatible with the rest of the codebase's errors.

**Fix:**

1. Add a new ProblemType in `src/zndraw_joblib/exceptions.py`:

   ```python
   class ProviderExecutionFailed(ProblemType):
       """The provider failed to execute the requested read."""

       title: ClassVar[str] = "Bad Request"
       status: ClassVar[int] = 400
   ```

2. Executor catch clause in `src/zndraw/providers/executor.py`:

   ```python
   except Exception as err:
       log.exception("InternalProviderExecutor failed for %s", provider_cls.__name__)
       problem = ProviderExecutionFailed.create(
           detail=f"{type(err).__name__}: {err}",
       )
       headers = {
           "Authorization": f"Bearer {token}",
           "X-Request-Hash": request_id,
           "X-Result-Status": "error",
           "Content-Type": "application/problem+json",
       }
       client.post(
           url,
           json=problem.model_dump(exclude_none=True),  # httpx serializes dict → JSON
           headers=headers,
       )
   ```

   No manual `.encode()`. `log.exception` uses the already-imported `log` at module top.

3. `upload_result` in `src/zndraw_joblib/router.py:1300-1346` reads the `X-Result-Status` header (`Annotated[str | None, Header()] = None`). When `"error"`, store the status alongside the body. Implementation choice: paired Redis key `{cache_key}:status` with same TTL, cleaner to debug than a prefix byte.

4. `read_provider` in `src/zndraw_joblib/router.py:1180-1268`:
   - Read `{cache_key}:status` alongside `cache_key`.
   - When status == `"error"`: parse body as JSON (`ProblemDetail` shape), extract the `status` field, return `Response(content=body, media_type="application/problem+json", status_code=problem["status"])`. If body is non-JSON (corrupted cache) fall back to 500.
   - Remove both content-sniffing branches.

**Files:** `src/zndraw_joblib/exceptions.py`, `src/zndraw/providers/executor.py`, `src/zndraw_joblib/router.py`.

**Test:** (1) executor failure → client receives `application/problem+json` body matching `ProblemDetail` schema, HTTP 400; (2) provider returning `{"type": "object", "error": null}` as a legitimate JSON payload → client receives it unaltered at `application/json`; (3) `caplog` shows `log.exception` for executor failures.

---

### 2.3 `_resolve_provider` mirrors LIST semantics (Important — security)

**Defect:** `src/zndraw_joblib/router.py:996-1000` — `if room_id not in ("@global",) and provider_room_id not in (...)` short-circuits when caller `room_id == "@global"`. A caller using the `@global` scope can resolve providers from *any* room, including private room providers (e.g., `/v1/joblib/rooms/@global/providers/room-42:filesystem:local`). This contradicts `_room_provider_filter`'s policy (line 977) that `@global` sees only `@global` providers.

**Fix:** Mirror the LIST allowlist:

```python
async def _resolve_provider(session, provider_name, room_id):
    provider_room_id, category, name = provider_name.split(":", 2)  # existing split

    allowed = {"@global"} if room_id == "@global" else {"@global", "@internal", room_id}
    if provider_room_id not in allowed:
        raise ProviderNotFound.exception(
            detail=f"Provider '{provider_name}' not accessible from room '{room_id}'"
        )
    ...
```

**Files:** `src/zndraw_joblib/router.py:984-1005`.

**Test:** regression — a `@global` caller requesting `room-42:filesystem:local` → 404. Two assertions added to `tests/zndraw_joblib/test_providers.py`: `@global` scope cannot resolve `room-42` providers; `@global` scope cannot resolve `@internal` providers.

---

### 2.4 Gate `@internal:filesystem:*` by superuser (Important — security)

**Defect:** `src/zndraw_joblib/router.py:977-981` — `_room_provider_filter` adds `@internal` to every room's visible set. Any authenticated user in any room can call `read_provider` against `@internal:filesystem:FilesystemRead` and read files under `Settings.filebrowser_path`.

**Fix:**

1. New setting in `src/zndraw/config.py`:

   ```python
   filebrowser_require_superuser: bool = True
   """When True, @internal filesystem providers are accessible only
   to superusers. Flip to False to allow all authenticated users."""
   ```

   Secure default. Enforces the memory rule that features must be Pydantic-configurable.

2. Helper in `src/zndraw_joblib/router.py`:

   ```python
   def _require_internal_filesystem_access(
       provider: ProviderRecord, user: User, settings: Settings
   ) -> None:
       if (
           provider.room_id == "@internal"
           and provider.category == "filesystem"
           and settings.filebrowser_require_superuser
           and not user.is_superuser
       ):
           raise Forbidden.exception(
               detail="@internal filesystem access requires superuser"
           )
   ```

3. Three enforcement points:
   - `_list_providers`: post-filter results in Python — cleaner than extending the SQL filter with an additional conditional. Drops gated rows when user is non-superuser and the flag is on.
   - `get_provider_info` after `_resolve_provider` call.
   - `read_provider` after `_resolve_provider`, before dispatch.

**Files:** `src/zndraw/config.py`, `src/zndraw_joblib/router.py` (three handlers).

**Test:** non-superuser → 403 on `get_provider_info` and `read_provider` for `@internal:filesystem:*`; empty list of `@internal:filesystem:*` in list endpoint. Superuser → 200 on all three. Flag `False` → non-superuser sees all three (back to current behavior).

---

### 2.6 `mock_executor` stub adds `token: str` (Important)

**Defect:** `tests/zndraw_joblib/test_registry.py:57-63` — `mock_executor` has 4 params; `InternalExecutor` Protocol declares 5 including `token: str` (`registry.py:31-38`). Pyright diagnostic at 5 call sites.

**Fix:** Add `token: str` to the stub signature. Body remains `pass`. The stub exists deliberately to satisfy the Protocol in package-isolation tests (no concrete executor in `zndraw_joblib`).

**Files:** `tests/zndraw_joblib/test_registry.py:57-63`.

**Test:** pyright passes; existing tests unchanged.

---

### 2.7 `StorageResultBackend` applies `key_prefix` (Important)

**Defect:** `src/zndraw/result_backends.py:151-185` — `RedisResultBackend` prefixes every key via `_k(key)`, but `StorageResultBackend.store/get/delete` pass the raw `key` to `self._storage[key]`. Two zndraw servers sharing a `FrameStorage` backend (Mongo/LMDB) will collide on cached frame results despite isolated `result_backend_key_prefix` settings.

**Fix:** Parity with `RedisResultBackend`:

```python
class StorageResultBackend:
    def __init__(self, storage: FrameStorage, key_prefix: str = "") -> None:
        self._storage = storage
        self._key_prefix = key_prefix

    def _k(self, key: str) -> str:
        return f"{self._key_prefix}:{key}" if self._key_prefix else key

    async def store(self, key: str, data: bytes, ttl: int) -> None:
        io = self._storage[self._k(key)]
        ...
    # apply _k in get/delete too
```

Wire from `Settings.result_backend_key_prefix` at composition time (wherever the backend is built in lifespan).

**Files:** `src/zndraw/result_backends.py`.

**Test:** two `StorageResultBackend` instances share a memory `FrameStorage` with different prefixes; store under the same `key` → no collision.

---

### 2.8 Move `LoadFile` import to top of `modifiers.py` (Minor)

**Defect:** `src/zndraw/extensions/modifiers.py:266` — lazy import with comment "avoid a cycle via zndraw_joblib.client". No cycle exists (verified: `filesystem.py` imports only from `zndraw.extensions.abc`, which imports only pydantic + stdlib). The lazy import + stale comment is vestigial.

**Fix:** Move `from zndraw.extensions.filesystem import LoadFile` to the top with the other imports, delete the comment and `# noqa: E402`. Keep `modifiers[LoadFile.__name__] = LoadFile` at the bottom where the dict is finalized.

**Files:** `src/zndraw/extensions/modifiers.py`.

**Test:** existing tests exercise the module; no new test needed.

---

### 2.9 Split `filebrowser_path` into `filebrowser_enabled: bool` + `filebrowser_path: str` (Minor)

**Defect:** `src/zndraw/config.py:27-31, 96` — `env_parse_none_str="none"` is model-wide. Legitimate `"none"` values in any other string field (e.g., `guest_password`) silently become `None`. Docstring warns but "audit new fields" is fragile.

**Fix:** Eliminate the sentinel entirely:

```python
# src/zndraw/config.py
class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_SERVER_",
        pyproject_toml_table_header=("tool", "zndraw", "server"),
        # env_parse_none_str removed
    )
    ...
    filebrowser_enabled: bool = True
    filebrowser_path: str = "."
```

Drop `env_parse_none_str` from `model_config`. All callsites that check `settings.filebrowser_path is None` → `not settings.filebrowser_enabled`.

**Files:** `src/zndraw/config.py`, any `settings.filebrowser_path is None` callsites (grep).

**Test:** new test — env vars `ZNDRAW_SERVER_FILEBROWSER_ENABLED=false` disables the default provider; `ZNDRAW_SERVER_FILEBROWSER_PATH=/custom` + default enabled → provider uses custom path.

---

### 2.10 Delete weak test; rename remaining two for new flag (Minor)

**Defect:** `tests/zndraw/test_providers_filesystem.py:360-368` — `assert resp.status_code in (200, 401)` tautologically passes; comment "we just want to confirm boot" reveals it adds no coverage beyond the other two tests in the trio.

**Fix:**
- Delete `test_filebrowser_path_none_disables_default_provider` (lines 360-368).
- Rename the remaining two from `test_filebrowser_path_none_*` → `test_filebrowser_disabled_*` and update them to use `ZNDRAW_SERVER_FILEBROWSER_ENABLED=false` per #2.9.

**Files:** `tests/zndraw/test_providers_filesystem.py`.

---

### 2.11 No change — review misread (Minor)

**Finding:** `tests/zndraw_joblib/test_registry.py:188` uses `_DummyProvider.category: ClassVar[str] = "filesystem"`. Reviewer compared against `Extension.category: ClassVar[Category]` (StrEnum). But `Provider.category: ClassVar[str]` (`src/zndraw_joblib/provider.py:29`) is a plain string by design — providers live in the framework layer where open strings let third-party authors add categories without extending an enum.

**Decision:** no change. Test matches its actual base class.

---

## 3. Frontend Fixes

### 3.1 Delete temporary e2e spec files (Critical, downgraded to cleanup)

**Defect:** `frontend/e2e/dockview-layout.spec.ts` and `frontend/e2e/pr920-fixes.spec.ts` reference selectors (`data-sliver-state`, text `Drop to dock right`) that don't exist in the implementation. Frontend e2e tests are not run in CI, so these are misleading-docs masquerading as coverage. Per user: both files were temporary.

**Fix:** Delete both files.

**Files:** `frontend/e2e/dockview-layout.spec.ts`, `frontend/e2e/pr920-fixes.spec.ts`.

**Follow-up (out of scope):** wire e2e into CI in a separate PR.

---

### 3.2 + 3.6 Zustand `useDockviewApi` store replaces `sharedApi` singleton (Critical + Important, merged)

**Defect:**
- `frontend/src/panels/DockviewLayout.tsx:25` — module-level `let sharedApi: DockviewApi | null`. Mutable singleton, no subscription mechanism.
- `frontend/src/panels/PlotsBrowserPanel.tsx:42-46` — `setInterval(..., 50)` polls `getDockviewApi()` until the API binds. Workaround for consumers mounting before `onReady` fires.
- `frontend/src/panels/ViewerView.tsx:31-36` — `onLeaveRoom` cascade reads `getDockviewApi()`; when `DockviewLayout` cleanup nulls `sharedApi` before `ViewerView`'s effect attaches, the listener silently no-ops.

**Fix:** New Zustand store:

```ts
// frontend/src/stores/dockviewApiStore.ts
import { create } from "zustand";
import type { DockviewApi } from "dockview-react";

interface DockviewApiStore {
  api: DockviewApi | null;
  setApi: (api: DockviewApi | null) => void;
}

export const useDockviewApi = create<DockviewApiStore>((set) => ({
  api: null,
  setApi: (api) => set({ api }),
}));
```

- `DockviewLayout.onReady`: `useDockviewApi.getState().setApi(event.api)`.
- `DockviewLayout` cleanup: `useDockviewApi.getState().setApi(null)`.
- Component consumers: `const api = useDockviewApi((s) => s.api);` — re-renders on bind, `useEffect` dep on `api` subscribes cleanly.
- Non-component callers (socket handlers): `useDockviewApi.getState().api` — synchronous read.

Delete `let sharedApi` and the exported `getDockviewApi()` function. Migrate all 5-6 call sites.

**Files:**
- `frontend/src/panels/DockviewLayout.tsx` (remove singleton, wire store).
- `frontend/src/panels/PlotsBrowserPanel.tsx` (remove polling, subscribe to store).
- `frontend/src/panels/ViewerView.tsx` (subscribe to store).
- `frontend/src/hooks/useLeaveRoom.ts` (lazy resolution).
- `frontend/src/hooks/socketHandlers/figureHandlers.ts`, `chatHandlers.ts` (use `.getState()`).
- New file: `frontend/src/stores/dockviewApiStore.ts`.

**Test:** `PlotsBrowserPanel` mounted via direct URL before `onReady` fires — listener binds once API appears, no polling. Unit test on the store (set + subscribe).

---

### 3.3 Drop `@vitejs/plugin-react` dep, regenerate lock (Critical)

**Defect:** `frontend/package.json:77` installs `@vitejs/plugin-react` alongside `@vitejs/plugin-react-swc`. Only the SWC variant is imported (`vite.config.ts:2`). Dead dep pulls ~7 transitive Babel packages.

**Fix:** Remove `"@vitejs/plugin-react": "^4"` from `devDependencies`. Run `bun install` to regenerate `bun.lock`.

**Files:** `frontend/package.json`, `frontend/bun.lock`.

**Test:** `bun run build` + `bun run dev` succeed; `bunx depcheck` (if available) reports no orphaned devDeps.

---

### 3.4 + 3.12 `useDragHover(position)` hook + shared `shimmer` constant (Important + Minor, merged)

**Defect:** `ActivityBar.tsx`, `SidebarZone.tsx`, `BottomZone.tsx` each duplicate: `dragDepth` useRef, `onDragEnter`/`onDragOver`/`onDragLeave` with `setTimeout(0)` re-check, `keyframes shimmer` declaration.

**Fix:**

```ts
// frontend/src/panels/dragStyles.ts
import { keyframes } from "@mui/system";

export const shimmer = keyframes`
  0%   { background-color: rgba(25, 118, 210, 0.14); }
  50%  { background-color: rgba(25, 118, 210, 0.24); }
  100% { background-color: rgba(25, 118, 210, 0.14); }
`;
```

```ts
// frontend/src/panels/useDragHover.ts
export function useDragHover(position: BarPosition) {
  const dragDepth = useRef(0);
  const hoverBar = useAppStore((s) => s.hoverBar);
  const setHoverBar = useAppStore((s) => s.setHoverBar);
  const isHovered = hoverBar === position;

  const onDragEnter = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    if (++dragDepth.current === 1) setHoverBar(position);
  }, [position, setHoverBar]);

  const onDragLeave = useCallback(() => {
    dragDepth.current = Math.max(0, dragDepth.current - 1);
    if (dragDepth.current === 0 && hoverBar === position) {
      setTimeout(() => { if (dragDepth.current === 0) setHoverBar(null); }, 0);
    }
  }, [position, hoverBar, setHoverBar]);

  const onDragOver = useCallback((e: React.DragEvent) => e.preventDefault(), []);

  return { isHovered, dragHandlers: { onDragEnter, onDragLeave, onDragOver } };
}
```

Each bar component imports and uses the hook. `shimmer` imported from `dragStyles.ts`.

**Files:** new `frontend/src/panels/dragStyles.ts`, new `frontend/src/panels/useDragHover.ts`, updates to `ActivityBar.tsx`, `SidebarZone.tsx`, `BottomZone.tsx`.

**Test:** unit test on the hook — three mounted instances share `hoverBar` state correctly; entering position A while hovering position B reassigns cleanly.

---

### 3.5 Always render sliver (Important)

**Defect:** `frontend/src/panels/ActivityBar.tsx:137` — `if (visibleIcons.length === 0 && !isDragActive) return null;`. An empty bar disappears. Dragging an icon away leaves no drop target to drag a new one back.

**Fix:** Drop the `return null`. Always render the Box. When `visibleIcons.length === 0 && !isDragActive`, render at narrow sliver width (4px) with `bgcolor: "background.paper"`. When `isDragActive`, expand to full width with shimmer. When `isHovered`, highlight. CSS transition on width for smooth visual.

Same pattern applies to `SidebarZone.tsx` and `BottomZone.tsx` — consistency.

**Files:** `frontend/src/panels/ActivityBar.tsx`, `SidebarZone.tsx`, `BottomZone.tsx`.

**Test:** empty bar has non-zero bounding box. Component test: render empty ActivityBar, assert it is in the DOM.

---

### 3.7 Extract `resetDockview` + `ensureViewerPanel` helpers (Important)

**Defect:** Three code paths redeclare the viewer panel's shape (`id, component, title`):
- `DockviewLayout.tsx:31-37` — `addViewerPanel` (initial mount).
- `landingPage.tsx:219-230` — useEffect re-adds on `roomId` change.
- `landingPage.tsx:459-473` — Reset-layout MenuItem inline closes panels and re-adds viewer.

**Fix:** In `DockviewLayout.tsx`:

```ts
function addViewerPanel(api: DockviewApi) { /* existing private helper */ }

export function resetDockview(api: DockviewApi): void {
  for (const p of api.panels) p.api.close();
  addViewerPanel(api);
}

export function ensureViewerPanel(api: DockviewApi): void {
  if (!api.getPanel("viewer")) addViewerPanel(api);
}
```

- `landingPage.tsx:459-473` Reset-layout MenuItem → `resetDockview(api); useAppStore.getState().resetLayout();`
- `landingPage.tsx:219-230` useEffect → `ensureViewerPanel(api);`
- `onReady` keeps calling `addViewerPanel(api)` internally.

Viewer config (`id: "viewer"`, `component: "viewer"`) lives in one file.

**Files:** `frontend/src/panels/DockviewLayout.tsx`, `frontend/src/pages/landingPage.tsx`.

**Test:** Reset-layout closes all panels and the viewer reappears. Direct-URL navigation with viewer closed → `ensureViewerPanel` recreates it.

---

### 3.8 Shared `useFilesystemProviders()` hook (Important)

**Defect:** Two consumers query `["filesystemProviders", roomId]`:
- `frontend/src/hooks/useHasFilesystemProviders.ts:22` → `staleTime: 5_000`.
- `frontend/src/panels/FilesystemPanel.tsx:78` → no `staleTime` (default 0).

Per-observer `staleTime` semantics mean the Panel's `0` defeats the Hook's `5_000` dedup whenever the Panel is open.

**Fix:** Extract a single hook:

```ts
// frontend/src/hooks/useFilesystemProviders.ts
export function useFilesystemProviders() {
  const roomId = useAppStore((s) => s.roomId);
  const queryClient = useQueryClient();
  const result = useQuery({
    queryKey: ["filesystemProviders", roomId],
    queryFn: () => listProviders(roomId!, "filesystem"),
    enabled: !!roomId,
    retry: false,
    staleTime: 5_000,
  });
  useEffect(() => {
    if (!roomId) return;
    const handle = () => {
      queryClient.invalidateQueries({ queryKey: ["filesystemProviders", roomId] });
    };
    socket.on("providers_invalidate", handle);
    return () => { socket.off("providers_invalidate", handle); };
  }, [roomId, queryClient]);
  return result;
}
```

- `useHasFilesystemProviders` becomes a one-line wrapper.
- `FilesystemPanel.tsx` imports `useFilesystemProviders` instead of its inline `useQuery`.

**Files:** new `frontend/src/hooks/useFilesystemProviders.ts`, rewrite `frontend/src/hooks/useHasFilesystemProviders.ts`, update `frontend/src/panels/FilesystemPanel.tsx`.

**Test:** mount both consumers simultaneously → single network call; socket `providers_invalidate` event triggers refetch for both observers.

---

### 3.9 Biome-clean `src/panels` (Important)

**Defect:** `bunx biome check src/panels` reports 10 errors + 23 warnings (verified 2026-04-20).

**Fix (two passes):**

1. `bunx biome check --write src/panels` — auto-applies safe fixes (import order, `useOptionalChain`, `useTemplate`).
2. Manual fixes for remaining:
   - `PlotsBrowserPanel.tsx:62` unused biome-ignore — becomes moot after 3.2 (the polling hook is deleted).
   - `ChatPanel.tsx:207/221` exhaustive-deps — add missing deps, verify no re-render loop.
   - `FilesystemPanel.tsx:79/116` `noNonNullAssertion` — `roomId!` → early return `if (!roomId) return null;`. The query already gates on `enabled: !!roomId`.

**Files:** all files under `frontend/src/panels/`.

**Out of scope:** `PlotView.tsx` 14 `any` types — pre-existing from legacy `FigureWindow`, not introduced by PR #920. Follow-up.

**Test:** `bunx biome check src/panels` returns zero errors, zero warnings (excluding `PlotView.tsx` `any` warnings noted above).

---

### 3.10 Add `activityBarSlice` unit tests (Minor)

**Defect:** `bun run test` runs 2 tests total. `frontend/src/stores/slices/activityBarSlice.ts:162` has `moveIconToBar`, `dropIconOnPanel`, `resetLayout`, and the active-bar state machine, all untested.

**Fix:** ~5-8 vitest tests covering:
- `moveIconToBar(id, position)` — icon moves from one bar to another.
- `moveIconToBar(id, position, overIdx)` — icon inserts at index.
- `dropIconOnPanel(id)` — removes from bar.
- `resetLayout` — restores `initialState()`.
- `toggleActive(bar, panel)` — opens, closes, switches.

**Files:** new `frontend/src/stores/slices/__tests__/activityBarSlice.test.ts`.

---

### 3.11 Selector refactor in `landingPage.tsx:240-249` (Minor)

**Defect:** The `?panel=` deep-link effect reads store state via `useAppStore.getState()` inside an effect. Non-idiomatic; effect doesn't re-run when the store changes.

**Fix:** Replace the `getState()` block with selectors:

```ts
const activeLeft = useAppStore((s) => s.activeLeft);
const activeRight = useAppStore((s) => s.activeRight);
const activeBottom = useAppStore((s) => s.activeBottom);
const toggleActive = useAppStore((s) => s.toggleActive);

useEffect(() => {
  const panel = searchParams.get("panel") as PanelId | null;
  if (!panel || !(panel in PANELS)) return;
  const def = PANELS[panel];
  if (def.kind !== "tool" || def.default.bar === "editor") return;
  const activeKey = def.default.bar === "left" ? activeLeft
    : def.default.bar === "right" ? activeRight : activeBottom;
  if (activeKey !== panel) toggleActive(def.default.bar, panel);
}, [searchParams, activeLeft, activeRight, activeBottom, toggleActive]);
```

**Files:** `frontend/src/pages/landingPage.tsx:240-249`.

---

### 3.13 — 3.14 out of scope

- `PlotView.tsx` 14 `any` types — pre-existing legacy debt. Not blocking. Follow-up PR.
- Wire frontend e2e into CI — infrastructure work, separate PR.

---

## 4. Acceptance Criteria

1. Backend: `uv run pytest tests/zndraw tests/zndraw_joblib -q` — all green. New regression tests for 2.1, 2.2+2.5, 2.3, 2.4, 2.7 pass.
2. Frontend: `bun run test` — all green. New slice tests from 3.10 pass.
3. `bunx @biomejs/biome check src/panels` — zero errors, zero warnings (excluding `PlotView.tsx` `any`).
4. `uvx prek --all-files` — green.
5. Manual smoke (per CLAUDE.md, playwright on dev server):
   - Empty ActivityBar is visible as a sliver.
   - Reset-layout recreates viewer and closes all other panels.
   - Filesystem icon appears/disappears per server provider state.
   - Direct-URL navigation (`?panel=filesystem`) opens the panel once the API binds.
6. No regressions: existing tests pass, `uv run zndraw` boots clean in dev mode, frontend dev server (`bun run dev`) renders.

---

## 5. Non-Goals

- **E2E in CI**: flagged as follow-up.
- **`PlotView.tsx` `any` types**: pre-existing debt.
- **Alembic migrations**: `ProviderRecord.worker_id` nullability from PR #920 already shipped without a migration; nothing in this spec changes the schema.
- **Unshipped compat**: PR #920 has not merged to a stable tag; per the "no unshipped legacy" rule, free to break in-flight experiments.

---

## 6. Commit Plan

Single branch, multiple commits (one per numbered fix or merged pair) so the follow-up review can walk diffs:

```text
fix(providers): lazy-resolve worker_token for @internal branch only (#2.1)
fix(providers): ProviderExecutionFailed + preserve X-Result-Status header (#2.2+2.5)
fix(providers): _resolve_provider mirrors LIST semantics (#2.3 — security)
feat(providers): filebrowser_require_superuser gate (#2.4 — security)
test(providers): add token param to mock_executor stub (#2.6)
fix(backends): StorageResultBackend applies key_prefix (#2.7)
chore(extensions): remove vestigial lazy LoadFile import (#2.8)
refactor(config): split filebrowser_enabled / filebrowser_path (#2.9 + 2.10)
fix(frontend): remove temporary e2e spec files (#3.1)
refactor(frontend): useDockviewApi zustand store (#3.2+3.6)
chore(frontend): drop unused @vitejs/plugin-react dep (#3.3)
refactor(panels): useDragHover hook + shared shimmer (#3.4+3.12)
fix(panels): always render empty bar as sliver (#3.5)
refactor(panels): resetDockview + ensureViewerPanel helpers (#3.7)
refactor(hooks): shared useFilesystemProviders hook (#3.8)
chore(panels): biome-clean (#3.9)
test(store): activityBarSlice unit tests (#3.10)
refactor(landing): selector-based panel deep-link effect (#3.11)
```
