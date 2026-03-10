# Domain Pitfalls

**Domain:** Splitting monolithic Python modules and TypeScript React/Zustand files
**Researched:** 2026-03-05
**Overall confidence:** HIGH (based on direct codebase analysis)

## Critical Pitfalls

Mistakes that cause test suite failures, runtime errors, or force reverts.

### Pitfall 1: Breaking the Public Import Contract

**What goes wrong:** Tests and internal code import specific names from `zndraw.client` (not just `ZnDraw`). When `client.py` becomes `client/__init__.py`, every import path that references a name no longer exported from `__init__.py` breaks silently at import time or at first use.

**Why it happens:** The refactor moves classes/functions into submodules but the `__init__.py` re-exports are incomplete. Easy to miss private helpers like `_estimate_frame_size` or exception classes like `RoomLockedError` and `ZnDrawError`.

**Consequences:** 499 backend tests fail at import time. Production code in `src/zndraw/exceptions.py`, `src/zndraw/cli_agent/connection.py`, and `src/zndraw/accessors.py` all import directly from `zndraw.client`.

**Specific imports that must survive the split:**
- `from zndraw.client import ZnDraw` (13+ files)
- `from zndraw.client import RoomLockedError` (3 files: `exceptions.py`, `test_geometry_ownership.py`, `cli_agent/connection.py`)
- `from zndraw.client import ZnDrawError` (2 files: `test_client.py`, `cli_agent/connection.py`)
- `from zndraw.client import Sessions` (1 file: `test_sessions.py`)
- `from zndraw.client import _estimate_frame_size` (3 test uses in `test_client.py`)
- `from zndraw.client import atoms_to_json_dict` (1 file: `test_routes_frames.py`)
- `from zndraw.client import APIManager` (2 files: `accessors.py`, `cli_agent/rooms.py` -- both under `TYPE_CHECKING`)

**Prevention:**
1. Before touching any code, run `grep -rn "from zndraw.client import" src/ tests/` and record every unique import.
2. Create `client/__init__.py` that re-exports every name that was previously importable.
3. Run `uv run pytest tests/ -x` immediately after the structural change, before any further refactoring.
4. Use `__all__` in `client/__init__.py` to document the public surface explicitly.

**Detection:** Tests fail at collection time with `ImportError`. This is loud and immediate -- the danger is doing the split in a branch and not running tests until multiple files have changed.

**Phase:** Must be addressed in the first commit of the `client.py` split.

---

### Pitfall 2: Zustand Cross-Slice State Access Breaking After Split

**What goes wrong:** `sceneSlice` calls `get().roomId`, `get().acquireLock()`, `get().releaseLock()`, and `get().showSnackbar()` -- all from other slices (ConnectionSlice, LockSlice, UISlice). If the scene slice is split into `geometrySlice`, `selectionSlice`, `editingSlice`, and `drawingSlice`, every sub-slice must still receive the full `AppState` via `StateCreator<AppState, [], [], SubSlice>`. If the type parameter is narrowed to the sub-slice type, `get()` will not expose cross-slice methods.

**Why it happens:** Developers instinctively type `StateCreator<DrawingSlice>` for a drawing slice, not realizing that Zustand's `get()` type is derived from the first type parameter. The `AppState` intersection type is what makes cross-slice access work.

**Consequences:** TypeScript compilation errors (loud), or worse: using `(get() as any).acquireLock()` to suppress the error, hiding real bugs. Drawing mode and editing mode both acquire distributed locks -- if `acquireLock` silently returns `undefined` instead of a boolean, users can enter editing mode without a lock, causing data corruption.

**Prevention:**
1. Every new sub-slice must use `StateCreator<AppState, [], [], SubSlice>` as its type.
2. The `AppState` type must be updated to include all new sub-slices: `AppState = ConnectionSlice & PlaybackSlice & GeometrySlice & SelectionSlice & EditingSlice & DrawingSlice & LockSlice & UISlice`.
3. After splitting, verify that `enterDrawingMode`, `exitDrawingMode`, `enterEditingMode`, `exitEditingMode` still call `acquireLock`/`releaseLock` successfully by running E2E tests.

**Detection:** TypeScript compiler errors if types are correct. If types are suppressed with `any`, it fails silently at runtime -- user enters drawing mode without acquiring lock.

**Phase:** Must be addressed during `sceneSlice` split.

---

### Pitfall 3: Socket Event Handler Registration/Cleanup Mismatch

**What goes wrong:** `useSocketManager` registers 25 socket event handlers with `socket.on()` and must unregister all 25 with `socket.off()` in the cleanup function. If handlers are extracted into separate modules/hooks but the cleanup does not exactly mirror the registration, stale handlers accumulate on reconnect, causing duplicate event processing (double state updates, double API calls).

**Why it happens:** When splitting handlers into groups (e.g., `useChatHandlers`, `useGeometryHandlers`), each group registers its own handlers. If a child hook unmounts at a different time than the parent, or if one group forgets cleanup, stale handlers pile up. React StrictMode double-mount makes this worse -- the existing `cancelled` flag pattern must be preserved in each sub-hook.

**Consequences:** Duplicate geometry fetches, duplicate snackbar messages, race conditions where an old handler reads stale closure state and overwrites state from a newer handler. These are intermittent and hard to reproduce.

**Prevention:**
1. If extracting handlers into separate functions/modules (not separate hooks), keep the single `useEffect` with single registration/cleanup. The handlers can be defined in external files and imported, but registration must stay centralized.
2. If extracting into separate hooks, each hook must manage its own `socket.on()`/`socket.off()` lifecycle symmetrically. Write a test that counts registered listeners for each event and asserts exactly 1.
3. Do not split into separate `useEffect` calls within the same component -- the dependency arrays would diverge and cause independent re-registration cycles.

**Detection:** In development, React StrictMode double-mount reveals handler duplication. In production, watch for duplicate console logs on socket events. Add a dev-mode assertion: `if (socket.listeners(eventName).length > 1) console.error(...)`.

**Phase:** Must be addressed during `useSocketManager` split.

---

### Pitfall 4: Circular Import Between Split Modules (Python)

**What goes wrong:** `client.py` currently has `ZnDraw`, `APIManager`, `SocketManager`, `ZnDrawLock`, and helper functions in one file. `SocketManager` references `ZnDraw` (via `self.zndraw: ZnDraw`). `ZnDraw` creates `APIManager` and `SocketManager`. If these become separate files, `socket_manager.py` imports `ZnDraw` from `zndraw_class.py` and `zndraw_class.py` imports `SocketManager` from `socket_manager.py` -- circular import.

**Why it happens:** The classes were designed to live in one file. Bidirectional references between classes that are colocated become circular imports when separated.

**Consequences:** `ImportError: cannot import name 'X' from partially initialized module` at startup. All 499 tests fail.

**Prevention:**
1. Use `TYPE_CHECKING` guards for type-only references. `SocketManager.zndraw` is typed as `ZnDraw` but only needs the type at static analysis time -- at runtime it receives the instance via `__init__`.
2. Move the `ZnDraw` type annotation in `SocketManager` behind `if TYPE_CHECKING:` and use a string annotation `"ZnDraw"` or the existing `from __future__ import annotations`.
3. The file already uses `from __future__ import annotations` (line 12), which defers all annotation evaluation. Leverage this: types referenced only in annotations will not cause circular imports.
4. If runtime access is needed (e.g., `isinstance(x, ZnDraw)`), restructure to use a protocol or move the check to the module that owns the class.

**Detection:** Immediate `ImportError` on any import of the package. Loud failure.

**Phase:** Must be addressed during `client.py` split, specifically when separating `SocketManager` from `ZnDraw`.

---

### Pitfall 5: Re-export Chain Breakage (TypeScript)

**What goes wrong:** `store.tsx` re-exports `getActiveCurves` and `selectPreferredCurve` from `sceneSlice.ts`. Components import these from `store.tsx` (e.g., `Canvas.tsx` does `import { getActiveCurves, selectPreferredCurve, useAppStore } from "../store"`). If these helpers move to a new file (e.g., `geometrySlice.ts` or `utils/curves.ts`) but the re-export in `store.tsx` is not updated, the build breaks.

**Why it happens:** The re-export in `store.tsx` (lines 14-17) explicitly names `getActiveCurves` and `selectPreferredCurve` from `"./stores/slices/sceneSlice"`. If the source file changes, this import breaks.

**Consequences:** TypeScript/Vite build fails. E2E tests cannot run.

**Prevention:**
1. Update re-exports in `store.tsx` when moving helpers to new files.
2. Alternatively, have components import helpers directly from the new module instead of going through `store.tsx`. This removes the re-export indirection.
3. Run `bunx tsc --noEmit` after each file move to catch broken imports before running the full build.

**Detection:** TypeScript compilation error. Loud and immediate.

**Phase:** Must be addressed during `sceneSlice` split.

## Moderate Pitfalls

### Pitfall 6: Zustand Selector Stability After Slice Reorganization

**What goes wrong:** Components use fine-grained selectors like `useAppStore((state) => state.setGeometries)`. After splitting `SceneSlice` into sub-slices, if the property name or location changes (e.g., `setGeometries` moves from `SceneSlice` to `GeometrySlice`), the selector still works at the type level (because `AppState` is the intersection of all slices) but the property may accidentally be defined in two places if the old slice is not fully removed.

**Why it happens:** Zustand's `create<AppState>((...a) => ({ ...createSliceA(...a), ...createSliceB(...a) }))` merges all slices with object spread. If two slices define the same property, the last one wins silently. During an incremental split, a property may temporarily exist in both the old `SceneSlice` (partially cleaned) and the new `GeometrySlice`.

**Consequences:** The wrong implementation is used at runtime. For example, if `removeGeometry` exists in both old and new slices with slightly different cleanup logic (one removes selection, the other does not), the spread order determines which version runs. This produces intermittent bugs visible only when specific geometry operations are performed.

**Prevention:**
1. Split atomically: remove from `SceneSlice` and add to `GeometrySlice` in the same commit.
2. After splitting, add a dev-mode check: iterate `Object.keys` of each slice's initial state and assert no duplicates across slices.
3. Use TypeScript to enforce: if `GeometrySlice` and `SceneSlice` both define `setGeometries`, the intersection type `GeometrySlice & SceneSlice` will error if the signatures differ. But if signatures are identical, no error -- so rely on the key uniqueness check above.

**Detection:** Silent at compile time if signatures match. Behavioral regression in E2E tests if cleanup logic differs.

**Phase:** Must be addressed during `sceneSlice` split.

---

### Pitfall 7: useEffect Dependency Array Divergence

**What goes wrong:** The current `useSocketManager` has a dependency array with 28 items (lines 918-948). When splitting into sub-hooks, each sub-hook needs its own dependency array. If a dependency is forgotten, the hook uses a stale closure value. If a dependency is added that was not previously there, the hook re-runs more often than before, causing unnecessary socket re-registration.

**Why it happens:** The monolithic hook captures all Zustand setters in the dependency array. These setters are stable references (Zustand guarantees this), so the dependency array effectively only reacts to `roomId` and `isOverview` changes. But when splitting, a developer might add a reactive value (not a setter) to the dependency array of a sub-hook, causing it to re-register handlers on every state change.

**Consequences:**
- Missing dependency: stale handler reads old `roomId`, sends events to wrong room.
- Extra dependency: handler constantly re-registers, causing brief windows where no handler is active (between cleanup and re-registration), dropping socket events.

**Prevention:**
1. The primary approach in the existing code (registering handlers inside a single `useEffect` that depends on `roomId` and `isOverview`, plus stable setters) should be preserved even after extraction. Move handler *definitions* out of the hook, not the *registration*.
2. If separate hooks are needed, keep handler registration inside the hook and make each sub-hook depend only on `roomId`/`isOverview` plus its specific stable setters.
3. Use `eslint-plugin-react-hooks` (exhaustive-deps rule) to catch missing dependencies.
4. Verify with React DevTools Profiler that the hook does not re-run more often after the split.

**Detection:** Intermittent event drops or duplicate processing. Very hard to reproduce -- requires specific timing of room switches + socket events.

**Phase:** Must be addressed during `useSocketManager` split.

---

### Pitfall 8: Python Package `__init__.py` Import Order Causing Side Effects

**What goes wrong:** When converting `client.py` to `client/__init__.py` + submodules, the import order in `__init__.py` matters if submodules have module-level side effects (logging configuration, module-level exception registration, etc.). The current `client.py` defines `log = logging.getLogger(__name__)` at module level and registers exception classes. If submodules are imported in the wrong order, `__name__` resolves differently and logging configuration breaks.

**Why it happens:** In `client.py`, `__name__` is `zndraw.client`. In `client/socket_manager.py`, `__name__` is `zndraw.client.socket_manager`. Log messages shift to a different logger hierarchy, which may not be captured by existing log filters.

**Consequences:** Log messages from the client disappear or appear under unexpected logger names. Not a functional failure, but hampers debugging in production. Exception classes registered with the wrong module path may confuse error tracking tools.

**Prevention:**
1. Each submodule should use `logging.getLogger(__name__)` (standard practice). Accept that logger names will change.
2. If log filtering by module name is used anywhere (tests, production config), update the filter patterns.
3. Check if any code does `logging.getLogger("zndraw.client")` explicitly -- this would need updating to `"zndraw.client.socket_manager"` etc.

**Detection:** Missing log output in test runs or production. Grep for `getLogger("zndraw.client")` to find hardcoded references.

**Phase:** Should be verified after `client.py` split, but is not blocking.

---

### Pitfall 9: SocketManager's Tight Coupling to ZnDraw Instance

**What goes wrong:** `SocketManager` directly accesses `self.zndraw.url`, `self.zndraw.api.token`, `self.zndraw.room`, `self.zndraw.copy_from`, `self.zndraw.api.session_id`, `self.zndraw._cached_length`, `self.zndraw._mount`, `self.zndraw._mount_name`, and `self.zndraw.api.create_room()`. This is 10+ attribute accesses on the `ZnDraw` instance, including private attributes (`_cached_length`, `_mount`, `_mount_name`). Simply putting `SocketManager` in a separate file does not reduce coupling -- it just spreads the coupling across files.

**Why it happens:** `SocketManager` was designed as an internal component of `ZnDraw`, not as an independent module. The split separates files but not concerns.

**Consequences:** The split creates two tightly coupled files that must change together. Any refactoring of `ZnDraw`'s internal state (`_cached_length`, `_mount`) requires changing `SocketManager` too. This is worse than the monolith because now the coupling is hidden across files.

**Prevention:**
1. Accept tight coupling for now. The goal is file organization, not architectural decoupling. Document that `SocketManager` is an internal implementation detail of `ZnDraw`.
2. Do not try to introduce abstractions (protocols, interfaces) to decouple them during this refactor -- that is scope creep that changes behavior.
3. Keep both in the same `client/` package to signal they are related.

**Detection:** No runtime issue -- this is a design quality concern. If a future refactoring of `ZnDraw` breaks `SocketManager`, the test suite will catch it.

**Phase:** Acknowledge during `client.py` split. Do not try to fix coupling in this milestone.

## Minor Pitfalls

### Pitfall 10: Relative Import Path Changes in Test Files

**What goes wrong:** Test files like `test_client.py` use local imports inside test functions: `from zndraw.client import _estimate_frame_size`. After the split, this helper may live in `zndraw.client.serialization` or similar. The tests must be updated.

**Why it happens:** Tests import private helpers directly rather than through the public API. Private helpers are the most likely to move during a split.

**Prevention:**
1. Update test imports in the same commit as the module restructure.
2. Consider re-exporting `_estimate_frame_size` from `client/__init__.py` for backward compatibility, since tests depend on it.
3. Run `uv run pytest tests/ -x --co` (collect-only) to verify all imports resolve before running the full suite.

**Detection:** `ImportError` during test collection. Loud and immediate.

**Phase:** Must be addressed alongside `client.py` split.

---

### Pitfall 11: Git Diff Noise from File Moves

**What goes wrong:** When a 2297-line file is deleted and its contents are split into 4-5 new files, `git diff` and `git blame` lose history. Code review becomes difficult because the PR shows the entire old file as deleted and all new files as added, even though most code is unchanged.

**Why it happens:** Git detects renames/moves only if the similarity threshold is met (default 50%). If a 2297-line file is split into 5 files of ~400-500 lines each, no single file has >50% similarity to the original.

**Prevention:**
1. Split in multiple commits: first commit creates `client/__init__.py` that imports from `client.py` (no code change, just the structural setup). Second commit moves one class at a time. This gives Git better rename detection per commit.
2. Use `git log --follow` for blame history after the split.
3. In the PR description, note which sections of `client.py` map to which new files.

**Detection:** Not a runtime issue. Makes code review harder.

**Phase:** Applies to all three file splits.

---

### Pitfall 12: Vite/TypeScript Module Resolution After File Moves

**What goes wrong:** Vite caches module resolution. After moving `sceneSlice.ts`, the dev server may serve stale modules until the cache is cleared. TypeScript's incremental compilation (`.tsbuildinfo`) may also cache old paths.

**Why it happens:** Vite's HMR and TypeScript's incremental compilation both cache file paths. Moving a file does not always invalidate the cache.

**Prevention:**
1. Restart the Vite dev server after file moves.
2. Delete `node_modules/.vite` cache if stale imports persist.
3. Delete `tsconfig.tsbuildinfo` if TypeScript reports phantom errors.
4. This is a development-time annoyance, not a production issue (production builds from scratch).

**Detection:** Phantom "module not found" errors in dev mode that do not reproduce in a clean build.

**Phase:** Development workflow concern during all frontend splits.

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| `client.py` split | Circular imports between `ZnDraw` and `SocketManager` (Pitfall 4) | Use `from __future__ import annotations` (already present) + `TYPE_CHECKING` guards |
| `client.py` split | Broken imports in 20+ files (Pitfall 1) | Comprehensive `__init__.py` re-exports; run tests after first commit |
| `client.py` split | Private attribute coupling (Pitfall 9) | Accept coupling; do not introduce abstractions |
| `sceneSlice.ts` split | Cross-slice `get()` type narrowing (Pitfall 2) | All sub-slices must use `StateCreator<AppState, [], [], SubSlice>` |
| `sceneSlice.ts` split | Duplicate property keys from partial migration (Pitfall 6) | Atomic moves: remove from old slice and add to new in same commit |
| `sceneSlice.ts` split | Re-export chain in `store.tsx` (Pitfall 5) | Update `store.tsx` re-exports when helpers move |
| `useSocketManager.ts` split | Handler registration/cleanup mismatch (Pitfall 3) | Keep registration centralized even if handler definitions are extracted |
| `useSocketManager.ts` split | Dependency array divergence (Pitfall 7) | Preserve the current pattern: stable setters + two reactive values |
| All phases | Git blame loss (Pitfall 11) | Multi-commit strategy: structure first, move code second |

## Verification Checklist (run after each split)

**Python (`client.py`):**
```bash
# 1. Verify all imports resolve
uv run python -c "from zndraw.client import ZnDraw, APIManager, SocketManager, ZnDrawLock, RoomLockedError, ZnDrawError, NotConnectedError, Sessions, atoms_to_json_dict, _estimate_frame_size"

# 2. Run full test suite
uv run pytest tests/ -x

# 3. Type check
uv run pyright .
```

**TypeScript (`sceneSlice.ts`):**
```bash
# 1. Type check
cd frontend && bunx tsc --noEmit

# 2. Build
cd frontend && bun run build

# 3. Run E2E tests (especially geometry-drawing, editing specs)
cd frontend && bun run playwright test
```

**TypeScript (`useSocketManager.ts`):**
```bash
# 1. Type check
cd frontend && bunx tsc --noEmit

# 2. Build
cd frontend && bun run build

# 3. Run E2E tests (especially socket-sync, frame-invalidation specs)
cd frontend && bun run playwright test
```

## Sources

- Direct codebase analysis of `src/zndraw/client.py` (2297 lines, 4 classes, 6 helper functions)
- Direct codebase analysis of `frontend/src/hooks/useSocketManager.ts` (949 lines, 25 socket handlers, 28-item dependency array)
- Direct codebase analysis of `frontend/src/stores/slices/sceneSlice.ts` (596 lines, 50+ methods, cross-slice `get()` calls to `roomId`, `acquireLock`, `releaseLock`, `showSnackbar`)
- Import graph analysis via grep across `src/` and `tests/` directories
- `.planning/codebase/CONCERNS.md` (existing tech debt documentation)
- `.planning/codebase/TESTING.md` (test infrastructure documentation)
- Zustand slice pattern from `frontend/src/store.tsx` (5 slices composed via spread)

---

*Pitfalls research: 2026-03-05*
