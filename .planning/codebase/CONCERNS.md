# Codebase Concerns

**Analysis Date:** 2026-03-05

## Tech Debt

**Client module is excessively large (2297 lines):**
- Issue: `src/zndraw/client.py` is 2297 lines, combining API wrapper, socket manager, lock context manager, and 10+ accessor classes in a single file.
- Files: `src/zndraw/client.py`
- Impact: Difficult to maintain, test in isolation, or extend. Contributors must navigate a monolithic file. Contains 149 function definitions.
- Fix approach: Extract `APIManager`, `SocketManager`, `ZnDrawLock`, and accessor classes (`Selections`, `Bookmarks`, `Geometries`, `Figures`, etc.) into separate modules under `src/zndraw/client/`.

**`useSocketManager` hook is a 949-line monolith:**
- Issue: `frontend/src/hooks/useSocketManager.ts` registers 25+ socket event handlers, manages auth/reconnect logic, room join/leave lifecycle, and error recovery in a single `useEffect`. This makes it fragile to modify.
- Files: `frontend/src/hooks/useSocketManager.ts`
- Impact: Any change to one handler risks breaking others. The massive dependency array (30+ items) makes it hard to reason about re-render triggers.
- Fix approach: Extract event handlers into separate hooks or a socket event registry pattern. Group related handlers (e.g., chat events, geometry events, frame events) into dedicated modules.

**`sceneSlice` carries too much responsibility (596 lines):**
- Issue: `frontend/src/stores/slices/sceneSlice.ts` manages geometries, selections, selection groups, drawing mode, editing mode, transform controls, camera attachment, curve refs, frame edits, and dynamic positions in one slice.
- Files: `frontend/src/stores/slices/sceneSlice.ts`
- Impact: State changes in one concern (e.g., editing mode) trigger re-renders in components that only care about another concern (e.g., selections). The slice interface has 50+ methods/properties.
- Fix approach: Split into dedicated slices: `geometrySlice`, `selectionSlice`, `editingSlice`, `drawingSlice`. Zustand supports this well with the existing `StateCreator` pattern.

**Broad `except Exception` catches in backend:**
- Issue: Nine instances of `except Exception` (bare catch-all) across the backend code.
- Files: `src/zndraw/executor.py` (lines 77, 83), `src/zndraw/routes/trajectory.py` (line 285), `src/zndraw/cli.py` (line 188), `src/zndraw/client.py` (lines 224, 238), `src/zndraw/accessors.py` (lines 544, 614), `src/zndraw/cli_agent/gif.py` (line 368)
- Impact: Swallows unexpected exceptions, making debugging harder. Can mask real bugs.
- Fix approach: Replace with specific exception types. For executor/client code, catch `httpx.HTTPError`, `socketio.exceptions.ConnectionError`, etc.

**74 `# type: ignore` comments across backend:**
- Issue: 74 `# type: ignore` comments spread across 20 files, primarily in Redis calls and SQLModel queries.
- Files: `src/zndraw/dependencies.py` (9), `src/zndraw/storage/asebytes_backend.py` (8), `src/zndraw/routes/geometries.py` (8), `src/zndraw/routes/screenshots.py` (8), `src/zndraw/socketio.py` (7)
- Impact: Pyright cannot catch real type errors in these locations. The Redis async stub issues are known but mask other problems.
- Fix approach: For Redis calls, create typed wrapper functions that handle the `decode_responses=True` return type properly. For SQLModel, use proper type annotations on session.exec results.

**Editable path dependencies in `pyproject.toml`:**
- Issue: `[tool.uv.sources]` contains `asebytes = { path = "../asebytes", editable = true }` and `zndraw-joblib = { path = "../zndraw-joblib", editable = true }`.
- Files: `pyproject.toml` (lines 209-211)
- Impact: These relative paths break in CI, worktrees, and when other developers clone the repo. Only works on the maintainer's local machine.
- Fix approach: Move editable overrides to `uv.toml` (local-only, gitignored) or use `--override` at install time.

## Known Bugs

**Missing FK CASCADE on most Room child tables:**
- Symptoms: Orphaned rows left behind when a room is deleted manually from DB.
- Files: `src/zndraw/models.py` (lines 58, 69, 80, 90, 98, 106, 137)
- Trigger: Only `Screenshot` has `ondelete="CASCADE"` (line 117). `Message`, `RoomMembership`, `RoomGeometry`, `RoomBookmark`, `SelectionGroup`, `RoomFigure`, and `RoomPreset` use plain `foreign_key="room.id"` without cascade.
- Workaround: No room deletion endpoint exists currently, but if one is added or manual DB cleanup is done, orphans will accumulate.

**No room deletion endpoint:**
- Symptoms: Rooms accumulate indefinitely. No way to clean up old rooms except direct DB access.
- Files: `src/zndraw/routes/rooms.py` (no DELETE endpoint), `src/zndraw/routes/admin.py`
- Trigger: Over time, the room table grows without bound.
- Workaround: Manual database cleanup.

## Security Considerations

**Socket.IO CORS set to wildcard:**
- Risk: `cors_allowed_origins="*"` allows any domain to establish WebSocket connections, enabling CSRF-like attacks on authenticated sessions.
- Files: `src/zndraw/socketio.py` (line 50)
- Current mitigation: JWT token required on connect, so connections still need valid auth.
- Recommendations: Set `cors_allowed_origins` to the actual frontend origin(s) via settings. Even with JWT, wildcard CORS enables credential theft from XSS on third-party sites.

**Default passwords in settings:**
- Risk: `guest_password` defaults to `"zndraw"` and `worker_password` defaults to `"zndraw-worker"`. If deployed without overriding, anyone can authenticate as the internal worker (superuser).
- Files: `src/zndraw/config.py` (lines 28-29)
- Current mitigation: These are SecretStr fields and the worker user has a fixed email (`worker@internal.user`).
- Recommendations: Log a warning at startup if defaults are in use. Consider requiring explicit password configuration in production mode.

**JWT token stored in localStorage:**
- Risk: localStorage is accessible to any JavaScript running on the same origin, making it vulnerable to XSS attacks.
- Files: `frontend/src/utils/auth.ts` (lines 12, 72, 96, 182, 189)
- Current mitigation: Auto-acquireToken coalesces concurrent calls and validates server-side.
- Recommendations: Consider httpOnly cookies for token storage, or accept the trade-off with strong CSP headers.

**No rate limiting on API endpoints:**
- Risk: No rate limiting on auth endpoints (`/v1/auth/guest`, `/v1/auth/jwt/login`) or data-heavy endpoints (frame uploads, trajectory uploads).
- Files: All route files under `src/zndraw/routes/`
- Current mitigation: None.
- Recommendations: Add rate limiting middleware (e.g., `slowapi` or custom Redis-based limiter). Priority: auth endpoints first to prevent brute-force attacks.

**No CSRF protection on REST API:**
- Risk: JWT Bearer auth in `Authorization` header prevents cookie-based CSRF, but the guest login endpoint (`POST /v1/auth/guest`) requires no credentials and creates user accounts.
- Files: `src/zndraw/routes/auth.py`
- Current mitigation: Guest accounts have limited privileges.
- Recommendations: Consider rate-limiting guest account creation.

## Performance Bottlenecks

**SmilesEditDialog bundle is 22MB (uncompressed):**
- Problem: The Ketcher molecule editor produces a 22MB JavaScript chunk. Even with lazy loading, this is extremely large.
- Files: `frontend/src/components/jsonforms-renderers/SmilesEditDialog.tsx`, built as `src/zndraw/static/assets/SmilesEditDialog-BReLlaVd.js`
- Cause: Ketcher bundles its entire chemical drawing engine including SVG rendering, bond calculations, and chemical validation.
- Improvement path: Verify the chunk is truly lazy-loaded (only fetched when user opens SMILES editor). Consider loading Ketcher from a CDN. Total JS payload for a room view: ~10MB uncompressed (index 2MB + three 1.2MB + plotly 4.7MB + mui 817K).

**Plotly chunk is 4.7MB uncompressed:**
- Problem: Plotly.js is a heavyweight charting library loaded as a separate chunk.
- Files: `frontend/src/components/FigureWindow.tsx`, built as `src/zndraw/static/assets/plotly-CEtMjGGD.js`
- Cause: Full Plotly distribution bundled despite using only basic chart types.
- Improvement path: Use `plotly.js-basic-dist-min` instead of `plotly.js-dist-min` to reduce size by ~60%. Or use `plotly.js-cartesian-dist-min` if only 2D charts are needed.

**Full geometry list fetched on every socket invalidation:**
- Problem: When `onSelectionsInvalidate` fires, it calls `listGeometries(roomId)` which fetches ALL geometries including their full config data, just to extract selection arrays.
- Files: `frontend/src/hooks/useSocketManager.ts` (lines 698-710)
- Cause: No dedicated "list selections only" endpoint exists.
- Improvement path: Add a `GET /v1/rooms/{room_id}/selections` endpoint that returns only `{key: selection[]}` without geometry configs. Or include only the changed selection in the socket event payload.

**SQLite locking serializes all DB writes:**
- Problem: In dev mode (default SQLite), `_apply_sqlite_locking` wraps all sessions with a single `asyncio.Lock`, serializing concurrent requests.
- Files: `src/zndraw/database.py` (lines 77-91)
- Cause: SQLite does not support concurrent writes.
- Improvement path: This is acceptable for development. Document that PostgreSQL should be used for production multi-user deployments. Add a startup log message when SQLite locking is active.

**`fakeredis` is a production dependency:**
- Problem: `fakeredis>=2.33.0` is in main dependencies (not dev-only), adding unnecessary weight to production Docker images that always have real Redis.
- Files: `pyproject.toml` (line 17), `src/zndraw/database.py` (line 215)
- Cause: Used as fallback when `REDIS_URL` is not configured (dev/single-user mode).
- Improvement path: Move to optional dependency group and make it conditional import. In Docker, always require `REDIS_URL`.

## Fragile Areas

**`useSocketManager` cleanup logic during room switching:**
- Files: `frontend/src/hooks/useSocketManager.ts` (lines 862-917)
- Why fragile: The cleanup function must distinguish between room switches, navigation to overview, and actual unmount. The `cancelled` flag, `effectRoomId`/`effectIsOverview` capture, and "do NOT emit room_leave during room switching" comment reveal past race conditions.
- Safe modification: Add tests covering room switch, overview navigation, and tab close scenarios. The cleanup relies on comparing captured values to current state which is inherently racy.
- Test coverage: E2E tests exist (`e2e/socket-sync.spec.ts`) but unit tests for cleanup logic are missing.

**Lock renewal interval timer lifecycle:**
- Files: `frontend/src/stores/slices/lockSlice.ts` (lines 44-85)
- Why fragile: `setInterval` for lock renewal stores its ID in Zustand state. If the component unmounts without proper cleanup, the interval keeps running and sending requests. The `stopLockRenewal` must be called explicitly.
- Safe modification: Ensure all mode exits (`exitDrawingMode`, `exitEditingMode`) call `releaseLock` which calls `stopLockRenewal`. Add a global cleanup on socket disconnect.
- Test coverage: `tests/test_edit_lock_integration.py` covers server-side. No frontend unit tests for timer lifecycle.

**Geometry invalidation handler has three code paths:**
- Files: `frontend/src/hooks/useSocketManager.ts` (lines 527-628)
- Why fragile: The `onGeometriesInvalidate` handler has separate paths for: delete-with-key, set-with-key (fetch single), and set-without-key (fetch all). The "fetch all" fallback fetches each geometry individually then rebuilds the store, creating a window where partial state is visible.
- Safe modification: Always include `key` in geometry invalidation events. Remove the "no key" fallback path.
- Test coverage: E2E tests exist but do not cover the race condition during bulk re-fetch.

**App state accessed via `useAppStore.getState()` inside event handlers:**
- Files: `frontend/src/hooks/useSocketManager.ts` (lines 131, 174, 262, 266, 353, 478, 504, 748-776)
- Why fragile: Direct store access via `getState()` bypasses React's render cycle. State reads inside event callbacks may see stale values if multiple events fire rapidly. This pattern is intentional (avoids closure stale captures) but makes reasoning about state consistency difficult.
- Safe modification: Accept this as a Zustand-idiomatic pattern. Document which state reads are "point-in-time snapshots" vs. "reactive."
- Test coverage: No unit tests for concurrent event handling.

## Scaling Limits

**In-memory SQLite as default database:**
- Current capacity: Single-process, single-user.
- Limit: No persistence across restarts. All rooms, users, messages lost on restart. No concurrent writes.
- Scaling path: Set `ZNDRAW_DATABASE_URL` to PostgreSQL (`postgresql+asyncpg://...`). Already supported by the codebase.

**Redis for all ephemeral state:**
- Current capacity: Single Redis instance handles cameras, locks, progress, provider caching.
- Limit: No Redis Cluster support configured. Provider result caching can fill memory if many providers are active.
- Scaling path: Redis key patterns in `src/zndraw/redis.py` use room-scoped keys, which would support sharding by room_id.

**Frame storage in-memory or LMDB:**
- Current capacity: `memory://` (default) or local LMDB file.
- Limit: Memory storage lost on restart. LMDB is local to one node.
- Scaling path: MongoDB storage backend exists via `asebytes[mongodb]`. Configure `ZNDRAW_STORAGE=mongodb://host/db`.

## Dependencies at Risk

**Ketcher at release candidate version:**
- Risk: `ketcher-core`, `ketcher-react`, `ketcher-standalone` are pinned to `3.9.0-rc.3` with resolution overrides in `package.json` (lines 9-12).
- Impact: RC versions may have breaking changes. The resolution override forces all transitive deps to use this version.
- Migration plan: Upgrade to stable Ketcher 3.9.x when released. Remove resolution overrides.

**Multiple overlapping charting libraries:**
- Risk: Both `plotly.js` + `plotly.js-dist-min` AND `recharts` are listed as dependencies.
- Impact: Bundle bloat. Two charting libraries that serve overlapping purposes.
- Files: `frontend/package.json` (lines 48-49, 58)
- Migration plan: Audit usage. If `recharts` is only used for the progress bar or simple charts, consolidate to one library.

**Both `@vitejs/plugin-react` and `@vitejs/plugin-react-swc` in devDependencies:**
- Risk: Two competing React transform plugins. Only `plugin-react-swc` is used in `vite.config.ts`.
- Files: `frontend/package.json` (lines 74-75), `frontend/vite.config.ts` (line 2)
- Impact: Extra install time. Potential confusion about which is active.
- Migration plan: Remove `@vitejs/plugin-react` from devDependencies.

## Missing Critical Features

**No room deletion API:**
- Problem: Rooms cannot be deleted through the API. Only admin user deletion exists.
- Blocks: Room cleanup, data management, GDPR compliance. Rooms and their associated data (frames, geometries, messages, screenshots) accumulate indefinitely.

**No frontend unit tests:**
- Problem: Zero unit tests for frontend components, hooks, stores, or utilities. Only E2E tests (13 Playwright specs) exist.
- Blocks: Safe refactoring of the 35K-line frontend codebase. State management logic in slices is complex (lock timers, selection equality, frame edit batching) and untested.

## Test Coverage Gaps

**Frontend state management untested:**
- What's not tested: All Zustand store slices (`sceneSlice`, `lockSlice`, `playbackSlice`, `uiSlice`, `connectionSlice`), the form store, and window manager store have zero unit tests.
- Files: `frontend/src/stores/slices/*.ts`, `frontend/src/formStore.ts`, `frontend/src/stores/windowManagerStore.ts`
- Risk: Lock renewal timer logic, selection comparison, frame edit batching, and mode transitions could break silently.
- Priority: High

**Frontend hooks and utilities untested:**
- What's not tested: `useSocketManager`, `useFrameBatch`, `useStepControl`, `useGeometryEditing`, `useFrameEditing`, `useCameraControls`, and all utilities (`msgpack-numpy.ts`, `transformProcessor.ts`, `geometryEditing.ts`, `colorUtils.ts`).
- Files: `frontend/src/hooks/*.ts`, `frontend/src/utils/*.ts`
- Risk: The `msgpack-numpy.ts` utility (485 lines) handles binary encoding/decoding critical to frame data transfer. Any regression here corrupts all frame data silently.
- Priority: High (especially `msgpack-numpy.ts`)

**API client error handling untested:**
- What's not tested: The axios interceptor retry logic in `frontend/src/myapi/client.ts` (401 retry, token refresh, lock token injection). The `WeakSet`-based retry guard and `acquireToken()` coalescing in `frontend/src/utils/auth.ts`.
- Files: `frontend/src/myapi/client.ts` (lines 91-141), `frontend/src/utils/auth.ts` (lines 155-176)
- Risk: Auth flow breakage would silently downgrade all users to disconnected state.
- Priority: Medium

**Socket event handler coverage:**
- What's not tested: The 25+ socket event handlers registered in `useSocketManager` have no isolated tests. E2E tests cover happy paths but not error conditions, race conditions, or partial failure scenarios.
- Files: `frontend/src/hooks/useSocketManager.ts`
- Risk: Concurrent events (e.g., geometry invalidation during room switch) may produce inconsistent state.
- Priority: Medium

**Backend accessibility:**
- What's not tested: The frontend has minimal ARIA attributes (33 occurrences across 15 files, mostly from MUI defaults). No accessibility testing suite exists.
- Files: All components under `frontend/src/components/`, `frontend/src/pages/`
- Risk: The 3D visualization (Three.js canvas) is inherently inaccessible. Non-3D UI elements (panels, dialogs, forms) rely on MUI defaults but have no explicit accessibility testing.
- Priority: Low (scientific visualization tool, not public-facing)

**`any` type proliferation in frontend:**
- What's not tested: 243 occurrences of `any` type across 58 frontend files. Particularly concentrated in `myapi/client.ts` (13), `useSocketManager.ts` (26), and `sceneSlice.ts` (12).
- Files: Concentrated in `frontend/src/hooks/useSocketManager.ts`, `frontend/src/myapi/client.ts`, `frontend/src/stores/slices/sceneSlice.ts`
- Risk: TypeScript's type safety is bypassed in critical data flow paths. Socket event payloads, geometry data, and frame data are all typed as `any`.
- Priority: Medium

---

*Concerns audit: 2026-03-05*
