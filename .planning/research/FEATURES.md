# Feature Landscape

**Domain:** Monolithic file refactoring (Python backend module, React hook, Zustand slice)
**Researched:** 2026-03-05

## Table Stakes

Features users expect. Missing = the refactor is incomplete or introduces regressions.

### client.py (2297 lines) -> `src/zndraw/client/` package

| Deliverable | Why Expected | Complexity | Notes |
|-------------|--------------|------------|-------|
| `__init__.py` re-exports `ZnDraw` and all public names | `from zndraw.client import ZnDraw` must keep working; `__init__.py` imports it | Low | Single source of truth for public API |
| Serialization helpers in own module (`serialization.py`) | `atoms_to_json_dict`, `json_dict_to_atoms`, `raw_frame_to_atoms`, chunking constants are pure functions with zero class coupling | Low | Lines 79-163. No state, no imports from other classes in this file |
| Exception classes in own module (`exceptions.py`) | `ZnDrawError`, `NotConnectedError`, `RoomLockedError` are standalone, referenced by multiple classes | Low | Lines 170-180. Already isolated by section comment |
| `ZnDrawLock` in own module (`lock.py`) | Self-contained dataclass; only depends on `APIManager` | Low | Lines 187-245. Clean boundary, single dependency |
| `APIManager` in own module (`api.py`) | 850+ lines of pure HTTP wrapper methods grouped by REST resource; clearly bounded by its class | Med | Lines 253-1136. Largest single class. Already has internal section headers marking concern groups |
| `SocketManager` in own module (`socket.py`) | Self-contained dataclass with 3 event handlers; depends on `ZnDraw` and `SyncClientWrapper` | Low | Lines 1143-1263. Small, clear boundary |
| `ZnDraw` main class in own module (`core.py`) | The MutableSequence facade that wires everything together | Med | Lines 1270-2298. Has forward reference to `SocketManager`, `APIManager`. Circular import risk needs `TYPE_CHECKING` |
| All 499 existing tests pass unchanged | Pure structural refactor -- zero behavioral change is the contract | High | Integration tests exercise the full stack; any import breakage fails loudly |
| `from zndraw import ZnDraw` keeps working | Public API contract through `__init__.py` | Low | Already just re-exports from `zndraw.client` |

### useSocketManager.ts (949 lines) -> handler group modules

| Deliverable | Why Expected | Complexity | Notes |
|-------------|--------------|------------|-------|
| Connection/lifecycle handlers extracted (`connectionHandlers.ts`) | `onConnect`, `onDisconnect`, `onConnectError`, room join logic -- cohesive lifecycle concern ~250 lines | Med | Most complex group: room creation on 404, version check, REST fetches for initial state |
| Frame handlers extracted (`frameHandlers.ts`) | `onFrameUpdate`, `onFramesInvalidate`, `onFrameSelectionUpdate` -- frame data concern | Low | ~80 lines, self-contained invalidation logic |
| Geometry handlers extracted (`geometryHandlers.ts`) | `onGeometriesInvalidate`, `onDefaultCameraInvalidate` -- geometry data concern | Low | ~110 lines, clear fetch-and-update pattern |
| Chat handlers extracted (`chatHandlers.ts`) | `onChatMessageNew`, `onChatMessageUpdated`, `onTyping` -- chat concern | Low | ~70 lines, direct query cache mutations |
| Scene invalidation handlers extracted (`sceneHandlers.ts`) | `onSelectionsInvalidate`, `onSelectionGroupsInvalidate`, `onBookmarksInvalidate`, `onFiguresInvalidate` -- bulk invalidation concern | Low | ~80 lines, all follow the `createInvalidateHandler` pattern or similar |
| Room/lock/progress handlers extracted (`roomHandlers.ts`) | `onRoomUpdate`, `onRoomDelete`, `onLockUpdate`, `onProgressStarted/Update/Complete` -- room-level events | Low | ~80 lines, straightforward store updates |
| Schema/invalidation handlers extracted (`queryHandlers.ts`) | `onInvalidate`, `onSchemaInvalidate` -- React Query cache concern | Low | ~20 lines, trivial |
| Main hook remains as orchestrator (`useSocketManager.ts`) | Registers/unregisters handlers, manages effect lifecycle, owns cleanup | Med | Becomes ~100-150 lines: imports handler factories, wires them to `socket.on/off`, owns the `useEffect` |
| Handler factory pattern (`createInvalidateHandler`) shared | Already exists in the file; must be importable by handler modules | Low | Move to a shared utils file or keep in main hook and pass as parameter |
| 13 Playwright E2E specs pass unchanged | Socket handler changes could break real-time sync in UI | High | E2E tests exercise socket flow end-to-end |

### sceneSlice.ts (596 lines) -> sub-slices

| Deliverable | Why Expected | Complexity | Notes |
|-------------|--------------|------------|-------|
| Geometry sub-slice (`geometrySlice.ts`) | `geometries`, `geometrySchemas`, `geometryDefaults`, `geometryUpdateSources`, `geometryFetchingStates`, `neededFrameKeys` + all geometry CRUD actions | Med | ~150 lines of state + ~80 lines of actions. Has cross-concern: `removeGeometry` also cleans selections, active curve, attached camera, editing callbacks |
| Selection sub-slice (`selectionSlice.ts`) | `selections`, `selectionGroups`, `updateSelections`, `updateSelectionForGeometry`, `loadSelectionGroup` | Low | ~80 lines. Makes REST calls to `updateSelectionAPI` |
| Editing sub-slice (`editingSlice.ts`) | `mode`, `transformMode`, `editingSelectedAxis`, `editingCombinedCentroid`, `editingCallbacks`, `pendingFrameEdits`, `editingFrameDataCount`, and enter/exit/cycle editing actions | Med | ~180 lines. Complex: enter/exit lock acquisition, frame edit save/auto-save, REST calls |
| Drawing sub-slice (`drawingSlice.ts`) | `drawingPointerPosition`, `drawingIsValid`, `activeCurveForDrawing`, enter/exit drawing mode | Med | ~100 lines. `enterDrawingMode` creates geometry via REST if none exist -- cross-concern with geometry slice |
| Camera sub-slice or merge into geometry | `attachedCameraKey`, `curveRefs`, `attachToCamera`, `registerCurveRef`, `unregisterCurveRef` | Low | ~40 lines. Small enough to merge into geometry slice |
| Shared helpers stay in scene barrel or utils | `getActiveCurves`, `selectPreferredCurve` -- re-exported from `store.tsx` | Low | Must remain importable from same path |
| `SceneSlice` interface composed from sub-interfaces | `SceneSlice = GeometrySlice & SelectionSlice & EditingSlice & DrawingSlice` | Low | Follows existing `AppState` composition pattern |
| `createSceneSlice` composes sub-creators | Matches how `AppState` composes `createConnectionSlice`, `createPlaybackSlice`, etc. | Med | Zustand slice composition with cross-slice access via `get()` |

## Differentiators

Features that go above and beyond a basic split. Not expected, but valuable.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Unit tests for extracted Python modules | Validates each module independently; catches import/wiring bugs before integration tests | Med | Focus on: `APIManager` method dispatch, `SocketManager` handler registration, serialization round-trips |
| TypeScript handler tests (vitest) | First frontend unit tests; proves handler logic in isolation without socket | Med | Project has zero frontend unit tests today. Would need vitest setup. High value for `onFramesInvalidate` invalidation logic |
| Typed handler signatures | Replace `data: any` parameters in socket handlers with typed interfaces | Low | Handlers currently use `any` for all socket payloads. Types exist in backend Pydantic models but are not shared |
| `APIManager` grouped by resource concern | Split `api.py` further: `api/frames.py`, `api/geometries.py`, `api/rooms.py`, etc. | Med | 850+ lines with 15+ resource groups. Each group is self-contained. Could use mixin pattern or separate classes |
| Index barrel files for handler groups | `hooks/handlers/index.ts` that re-exports all handler factories | Low | Cleaner imports in the orchestrator hook |
| Deprecation removal during split | `vis.log()` (deprecated in favor of `vis.chat.send()`) and `progress_bar` (deprecated in favor of `ZnDrawTqdm`) | Low | Two deprecated methods in ZnDraw class. Clean break is allowed per PROJECT.md |

## Anti-Features

Features to explicitly NOT build during this refactor.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Backwards-compatibility shims (`from zndraw.client import APIManager`) | PROJECT.md says "clean break allowed"; no external consumers of internal classes | Only maintain `from zndraw.client import ZnDraw` (the public API). Internal imports like `APIManager` move to `zndraw.client.api` |
| Functional changes to any module | Out of scope per PROJECT.md. Mixing behavior changes with structural refactor makes failures ambiguous | File issues for functional improvements; do them in separate PRs |
| Further splitting `APIManager` into per-resource files | Over-engineering for a refactor milestone. `api.py` at ~850 lines is within reasonable single-file size | Keep as single file; add section headers (already present). Split only if a future milestone demands it |
| Abstract base classes for handler groups | Handlers are not polymorphic; ABCs add indirection with zero benefit | Use plain functions or factory functions |
| Shared socket event type package (Python <-> TypeScript) | Requires codegen infrastructure, build pipeline changes, cross-repo concerns | Use inline TypeScript interfaces. Type sharing is a separate project |
| Moving `accessors.py` into the client package | Already extracted (628 lines in its own file); not part of this refactor scope | Leave `accessors.py` where it is |
| Splitting `useSocketManager` into multiple hooks | Multiple `useEffect` hooks with socket registration would cause registration ordering bugs and double-mount issues in React StrictMode | Keep single `useEffect` in one orchestrator hook; extract handler *functions*, not hooks |
| Converting Zustand slice pattern to separate stores | Would break the `AppState` composition and require all selectors to use different store references | Keep composing into single `AppState`; use sub-slice creators within `createSceneSlice` |
| Lazy imports within handler modules (frontend) | Tree-shaking handles dead code; lazy imports add complexity for zero benefit in bundled frontend code | Use standard top-level imports |
| Re-typing `geometries: Record<string, any>` during the split | Correct typing of geometry data is a separate concern; mixing it in makes the refactor PR unreviewable | Keep `any` types; file a follow-up for geometry typing |

## Feature Dependencies

```
client.py split:
  serialization.py        (no dependencies)
  exceptions.py           (no dependencies)
  lock.py                 -> api.py (uses APIManager)
  api.py                  -> exceptions.py (uses raise_for_status)
  socket.py               -> core.py (references ZnDraw via TYPE_CHECKING)
  core.py                 -> api.py, socket.py, exceptions.py, serialization.py
  __init__.py             -> core.py (re-exports ZnDraw)

useSocketManager.ts split:
  handler modules          -> shared types (handler factory type, store action types)
  useSocketManager.ts      -> all handler modules (imports and registers)

sceneSlice.ts split:
  geometrySlice.ts         (no cross-slice dependencies for state)
  selectionSlice.ts        -> geometrySlice (removeGeometry cleans selections)
  editingSlice.ts          -> lockSlice (acquireLock/releaseLock via get())
  drawingSlice.ts          -> geometrySlice (reads/creates geometries), lockSlice (acquireLock)
  createSceneSlice         -> all sub-slices (composes them)
```

## MVP Recommendation

Prioritize:
1. **client.py split** -- highest line count, most test coverage, lowest risk (Python module system is straightforward)
2. **sceneSlice.ts split** -- follows established Zustand slice pattern already used by other slices, clear concern boundaries
3. **useSocketManager.ts split** -- most nuanced (React effect lifecycle, cleanup ordering, StrictMode double-mount), but lowest structural risk since handlers are pure functions

Defer:
- **Further `APIManager` decomposition**: 850 lines is large but each method is a thin REST wrapper. Splitting by resource would create 10+ tiny files with identical patterns. Keep as single file unless a future milestone adds significant API surface.
- **Frontend unit test infrastructure**: Valuable but orthogonal to the split itself. Set up vitest in a separate milestone, then add tests for extracted handlers.
- **Deprecation removal**: `vis.log()` and `progress_bar` removal is trivial but should be a separate commit to keep the refactor diff clean.

## Sources

- Direct analysis of `src/zndraw/client.py` (2297 lines) -- HIGH confidence
- Direct analysis of `frontend/src/hooks/useSocketManager.ts` (949 lines) -- HIGH confidence
- Direct analysis of `frontend/src/stores/slices/sceneSlice.ts` (596 lines) -- HIGH confidence
- Existing slice pattern in `connectionSlice.ts` (53 lines), `lockSlice.ts` (177 lines), `playbackSlice.ts` (131 lines), `uiSlice.ts` (119 lines) -- HIGH confidence
- Existing `accessors.py` (628 lines) as precedent for extracted accessor pattern -- HIGH confidence
- `store.tsx` composition pattern (`AppState = ConnectionSlice & PlaybackSlice & SceneSlice & LockSlice & UISlice`) -- HIGH confidence
- PROJECT.md scope/constraints -- HIGH confidence
