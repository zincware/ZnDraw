# Architecture Patterns

**Domain:** Monolithic file decomposition in a FastAPI + React/Zustand application
**Researched:** 2026-03-05

## Recommended Architecture

Three monolithic files split into focused modules, each organized as a package (directory with `__init__.py` or `index.ts`) that re-exports the public API. The split preserves all existing import paths by re-exporting from the package init.

### File 1: `client.py` (2297 lines) -> `src/zndraw/client/` package

```
src/zndraw/client/
    __init__.py          # Re-exports: ZnDraw, APIManager, SocketManager,
                         #   ZnDrawLock, ZnDrawError, NotConnectedError,
                         #   RoomLockedError, atoms_to_json_dict,
                         #   json_dict_to_atoms, raw_frame_to_atoms,
                         #   _estimate_frame_size
    _serialization.py    # atoms_to_json_dict, json_dict_to_atoms,
                         #   raw_frame_to_atoms, _estimate_frame_size,
                         #   _TARGET_CHUNK_BYTES, _MAX_CHUNK_FRAMES
    _exceptions.py       # ZnDrawError, NotConnectedError, RoomLockedError
    _lock.py             # ZnDrawLock (edit lock context manager)
    _api.py              # APIManager (all REST API methods)
    _socket.py           # SocketManager (Socket.IO connection + handlers)
    _zndraw.py           # ZnDraw class (MutableSequence, properties, mount)
```

**Rationale:** Each module maps to exactly one class or one concern. The `_` prefix signals these are internal modules -- the public API is `from zndraw.client import ZnDraw`.

### File 2: `useSocketManager.ts` (949 lines) -> handler modules

```
frontend/src/hooks/
    useSocketManager.ts          # Orchestrator: registers/unregisters handlers,
                                 #   manages connect/disconnect, auth flow.
                                 #   Imports handler factories from subdirectory.
    socketHandlers/
        index.ts                 # Re-exports all handler factories
        connectionHandlers.ts    # onConnect, onDisconnect, onConnectError
        frameHandlers.ts         # onFrameUpdate, onFramesInvalidate,
                                 #   onFrameSelectionUpdate
        geometryHandlers.ts      # onGeometriesInvalidate,
                                 #   onDefaultCameraInvalidate,
                                 #   onSelectionsInvalidate
        roomHandlers.ts          # onRoomUpdate, onRoomDelete,
                                 #   onSelectionGroupsInvalidate,
                                 #   onBookmarksInvalidate
        chatHandlers.ts          # onChatMessageNew, onChatMessageUpdated,
                                 #   onTyping
        figureHandlers.ts        # onFiguresInvalidate
        lockHandlers.ts          # onLockUpdate
        progressHandlers.ts      # onProgressStarted, onProgressUpdate,
                                 #   onProgressComplete
        cacheHandlers.ts         # onInvalidate, onSchemaInvalidate
                                 #   (React Query invalidation)
```

**Rationale:** The grouping follows the domain resource each handler operates on. The orchestrator (`useSocketManager.ts`) shrinks from 949 to ~150 lines: just the `useEffect` skeleton that registers/unregisters all handlers plus the connect/auth sequence.

### File 3: `sceneSlice.ts` (596 lines) -> composable sub-slices

```
frontend/src/stores/slices/
    sceneSlice.ts                # Combines sub-slices, exports SceneSlice type.
                                 #   Re-exports getActiveCurves, selectPreferredCurve.
    sceneSlices/
        index.ts                 # Re-exports all sub-slice creators
        geometrySlice.ts         # Geometry CRUD state + actions:
                                 #   geometries, geometrySchemas, geometryDefaults,
                                 #   geometryUpdateSources, geometryFetchingStates,
                                 #   neededFrameKeys, set/update/remove geometry,
                                 #   curveRefs, hoveredGeometryInstance, particleCount,
                                 #   curveLength, loadedDynamicPositions, attachedCameraKey
        selectionSlice.ts        # Selection state + actions:
                                 #   selections, selectionGroups,
                                 #   set/update/load selections,
                                 #   updateSelectionForGeometry
        modeSlice.ts             # Mode transitions + drawing/editing state:
                                 #   mode, transformMode, editingSelectedAxis,
                                 #   drawingPointerPosition, drawingIsValid,
                                 #   editingCombinedCentroid, editingCallbacks,
                                 #   activeCurveForDrawing,
                                 #   enter/exit drawing/editing mode,
                                 #   cycleTransformMode, subscribeToEditing,
                                 #   notifyEditingChange
        frameEditSlice.ts        # Frame editing persistence:
                                 #   pendingFrameEdits, editingFrameDataCount,
                                 #   setPendingFrameEdit, clearPendingFrameEdits,
                                 #   saveFrameEdits, increment/decrement count
```

**Rationale:** Zustand slice composition is already the pattern this codebase uses (`AppState = A & B & C & ...`). The scene slice is itself composed of sub-slices using the same `StateCreator<AppState>` signature. The parent `sceneSlice.ts` spreads them together just like `store.tsx` spreads the top-level slices.

## Component Boundaries

### client.py decomposition

| Module | Responsibility | Depends On | Line Estimate |
|--------|---------------|------------|---------------|
| `_serialization.py` | Frame encode/decode (ase.Atoms <-> dict) | `asebytes`, `ase`, `base64` | ~90 |
| `_exceptions.py` | Client exception hierarchy | None | ~15 |
| `_lock.py` | Edit lock context manager with auto-refresh | `_api.py` (APIManager type) | ~60 |
| `_api.py` | All HTTP REST calls, error handling | `_exceptions.py`, `httpx`, server schemas | ~880 |
| `_socket.py` | Socket.IO connect/disconnect/handlers | `_api.py` (via ZnDraw ref), `socketio` | ~130 |
| `_zndraw.py` | Main client class, MutableSequence, properties | Everything above, `accessors.py` | ~620 |

**Key boundary rule:** `_api.py` has zero dependencies on `_zndraw.py` or `_socket.py`. `_socket.py` depends on `_zndraw.py` (via the `zndraw` back-reference). `_zndraw.py` composes `APIManager` and `SocketManager`.

### useSocketManager.ts decomposition

| Module | Responsibility | Depends On |
|--------|---------------|------------|
| `connectionHandlers.ts` | Auth, version check, room join/create, disconnect | `socket`, `connectWithAuth`, API client, all store setters |
| `frameHandlers.ts` | Frame navigation, cache invalidation | `queryClient`, `setCurrentFrame`, `setFrameCount` |
| `geometryHandlers.ts` | Geometry CRUD sync, selection sync | `queryClient`, API client, geometry store setters |
| `roomHandlers.ts` | Room metadata updates, bookmarks, selection groups | `queryClient`, API client, store setters |
| `chatHandlers.ts` | Chat message insertion, typing indicators | `queryClient`, store setters |
| `figureHandlers.ts` | Figure cache + window management | `queryClient`, `windowManagerStore` |
| `lockHandlers.ts` | Lock acquisition/release/expiry from other users | Store lock setters |
| `progressHandlers.ts` | Progress tracker lifecycle | Store progress setters |
| `cacheHandlers.ts` | Generic React Query invalidation (extension data, schemas) | `queryClient` |

**Key boundary rule:** Each handler module exports a factory function (or plain functions) that receives dependencies (roomId, queryClient, store setters) as parameters. No handler module imports from another handler module. The orchestrator is the only module that knows about all handlers.

### sceneSlice.ts decomposition

| Module | Responsibility | Cross-Slice Access |
|--------|---------------|-------------------|
| `geometrySlice.ts` | Geometry data, schemas, defaults, fetching states, curve refs | None outward |
| `selectionSlice.ts` | Per-geometry selections, selection groups, load/update | Reads `roomId` from connection slice |
| `modeSlice.ts` | View/drawing/editing mode transitions, lock acquisition | Reads `roomId`, calls `acquireLock`/`releaseLock` from lock slice, reads geometries |
| `frameEditSlice.ts` | Pending frame edits, save/clear | Reads `roomId`, calls `showSnackbar` from UI slice |

**Key boundary rule:** Sub-slices use `get()` to read cross-slice state (this is already the Zustand pattern). No sub-slice imports from another sub-slice directly -- they access sibling state through `AppState` via `get()`.

## Data Flow

### client.py internal data flow

```
User code
    |
    v
ZnDraw (_zndraw.py)
    |-- frame ops --> APIManager._api.py  --> HTTP --> Server
    |-- .selections --> Selections (accessors.py) --> APIManager --> HTTP
    |-- .mount() -----> SocketManager (_socket.py) --> Socket.IO --> Server
    |-- .run() -------> JobManager (zndraw-joblib) --> APIManager --> HTTP
    |-- .get_lock() --> ZnDrawLock (_lock.py) --> APIManager --> HTTP
```

Import direction is strictly: `_zndraw.py` -> `_socket.py` -> (no further internal deps).
`_api.py` is a leaf -- nothing inside `client/` depends on it except `_lock.py`, `_socket.py`, and `_zndraw.py`.

### useSocketManager.ts event flow

```
Socket.IO event arrives
    |
    v
useSocketManager.ts (orchestrator)
    |-- dispatches to --> handler function (from socketHandlers/*.ts)
        |-- calls --> API client function (REST fetch)
        |-- calls --> Zustand store setter (state update)
        |-- calls --> queryClient.invalidateQueries (cache bust)
```

All handlers are registered in one `useEffect`. The cleanup function unregisters all handlers. The orchestrator does NOT contain handler logic -- it delegates to imported handler functions.

### sceneSlice.ts state flow

```
Three.js component interaction
    |
    v
Zustand action (e.g., updateSelections, enterDrawingMode)
    |-- reads --> AppState via get() (cross-slice: roomId, acquireLock)
    |-- calls --> REST API (optimistic update pattern for selections)
    |-- sets  --> Zustand state via set() (triggers React re-render)
```

## Patterns to Follow

### Pattern 1: Handler Factory for Socket Events (useSocketManager)

Each handler module exports functions that accept dependencies, returning event handlers. This avoids closure-over-everything and makes handlers independently testable.

```typescript
// socketHandlers/frameHandlers.ts
import type { QueryClient } from "@tanstack/react-query";

interface FrameHandlerDeps {
  roomId: string | null;
  setCurrentFrame: (frame: number) => void;
  setFrameCount: (count: number) => void;
  queryClient: QueryClient;
}

export function createFrameUpdateHandler(deps: FrameHandlerDeps) {
  return (data: { frame: number }) => {
    const { setCurrentFrame } = deps;
    // ... handler logic
    setCurrentFrame(data.frame);
  };
}

export function createFramesInvalidateHandler(deps: FrameHandlerDeps) {
  return (data: FramesInvalidateEvent) => {
    // ... handler logic
  };
}
```

**Why:** Handlers become pure functions of their deps. Testable without React. The orchestrator (`useSocketManager.ts`) instantiates them with current deps and registers with `socket.on()`.

### Pattern 2: Zustand Sub-Slice Composition (sceneSlice)

Sub-slices use the same `StateCreator<AppState>` type and are spread into the parent slice.

```typescript
// sceneSlices/selectionSlice.ts
import type { StateCreator } from "zustand";
import type { AppState } from "../../../store";

export interface SelectionState {
  selections: Record<string, number[]>;
  selectionGroups: Record<string, Record<string, number[]>>;
  setSelections: (selections: Record<string, number[]>) => void;
  updateSelections: (geometryKey: string, id: number, isShiftPressed: boolean) => void;
  // ...
}

export const createSelectionSlice: StateCreator<AppState, [], [], SelectionState> = (
  set, get
) => ({
  selections: {},
  selectionGroups: {},
  setSelections: (selections) => set({ selections }),
  updateSelections: (geometryKey, id, isShiftPressed) => {
    const state = get();
    const roomId = state.roomId; // cross-slice access via get()
    // ... logic
  },
});
```

```typescript
// sceneSlice.ts (parent combiner)
import { createGeometrySlice } from "./sceneSlices/geometrySlice";
import { createSelectionSlice } from "./sceneSlices/selectionSlice";
import { createModeSlice } from "./sceneSlices/modeSlice";
import { createFrameEditSlice } from "./sceneSlices/frameEditSlice";

export type SceneSlice = GeometryState & SelectionState & ModeState & FrameEditState;

export const createSceneSlice: StateCreator<AppState, [], [], SceneSlice> = (
  ...a
) => ({
  ...createGeometrySlice(...a),
  ...createSelectionSlice(...a),
  ...createModeSlice(...a),
  ...createFrameEditSlice(...a),
});
```

**Why:** This is exactly how the top-level `store.tsx` already composes slices. Applying the same pattern one level deeper is natural and consistent.

### Pattern 3: Python Package with Re-export Init (client.py)

The `__init__.py` re-exports everything that was previously importable from `client.py`.

```python
# src/zndraw/client/__init__.py
"""ZnDraw Python client package."""

from zndraw.client._api import APIManager
from zndraw.client._exceptions import NotConnectedError, RoomLockedError, ZnDrawError
from zndraw.client._lock import ZnDrawLock
from zndraw.client._serialization import (
    _estimate_frame_size,
    _TARGET_CHUNK_BYTES,
    _MAX_CHUNK_FRAMES,
    atoms_to_json_dict,
    json_dict_to_atoms,
    raw_frame_to_atoms,
)
from zndraw.client._socket import SocketManager
from zndraw.client._zndraw import ZnDraw

__all__ = [
    "APIManager",
    "NotConnectedError",
    "RoomLockedError",
    "SocketManager",
    "ZnDraw",
    "ZnDrawError",
    "ZnDrawLock",
    "atoms_to_json_dict",
    "json_dict_to_atoms",
    "raw_frame_to_atoms",
]
```

**Why:** `from zndraw.client import ZnDraw` continues to work. `from zndraw import ZnDraw` continues to work (via `zndraw/__init__.py` which imports from `zndraw.client`). Zero import breakage.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Splitting by Size Instead of Cohesion
**What:** Splitting `_api.py` into `_api_frames.py`, `_api_geometries.py`, etc. because it is ~880 lines.
**Why bad:** APIManager is a single class with a single responsibility (HTTP client). Splitting it scatters related methods across files and forces readers to search for which file has `get_geometry()`. The 880-line size is fine because every method follows the same pattern (build request, call HTTP, check status, return).
**Instead:** Keep APIManager as one class in one file. It is large but simple -- every method is 5-10 lines of the same pattern.

### Anti-Pattern 2: Circular Dependencies Between Sub-Modules
**What:** `_socket.py` imports from `_zndraw.py`, and `_zndraw.py` imports from `_socket.py`.
**Why bad:** Circular imports cause `ImportError` at module load time.
**Instead:** `_socket.py` takes `ZnDraw` as a constructor parameter and accesses it via `self.zndraw`. The type annotation uses `TYPE_CHECKING` to avoid the circular import. This is already the pattern in the current code.

### Anti-Pattern 3: Handler Modules Importing Each Other
**What:** `geometryHandlers.ts` importing a utility from `frameHandlers.ts`.
**Why bad:** Creates coupling between handler modules. Changes to one handler ripple to others.
**Instead:** Extract shared utilities into a separate file (e.g., `socketHandlers/utils.ts` for `createInvalidateHandler` factory). Handler modules are leaves -- they import from utils and from external deps only.

### Anti-Pattern 4: Splitting the SceneSlice Interface Across Files Without a Combiner
**What:** Each sub-slice exports its own interface, and `store.tsx` imports all four interfaces.
**Why bad:** Leaks internal structure to the store composition layer. Adding a new sub-slice requires changing `store.tsx`.
**Instead:** `sceneSlice.ts` remains the single export point. It defines `SceneSlice = GeometryState & SelectionState & ModeState & FrameEditState` and exports `createSceneSlice`. `store.tsx` does not change at all.

## Dependency / Import Flow

### client.py package

```
External deps (ase, httpx, socketio, msgpack, pydantic, zndraw_joblib, zndraw_socketio)
    ^
    |
_serialization.py  (ase, asebytes, base64)
_exceptions.py     (standalone)
    ^
    |
_api.py            (_exceptions, httpx, zndraw.exceptions, zndraw.schemas)
    ^
    |
_lock.py           (_api type only, threading)
    ^
    |
_socket.py         (socketio, zndraw_socketio, zndraw.socket_events)
    ^               references _zndraw.ZnDraw via TYPE_CHECKING
    |
_zndraw.py         (_api, _socket, _lock, _serialization, _exceptions,
                     zndraw.accessors, zndraw_joblib)
    ^
    |
__init__.py         (re-exports from all modules)
```

**Direction:** Strictly top-down. No module imports from a module below it in this diagram. `_api.py` is the most independent (leaf). `_zndraw.py` is the most dependent (root).

### useSocketManager.ts handler modules

```
External deps (socket, API client, Zustand stores, React Query)
    ^
    |
socketHandlers/utils.ts          (createInvalidateHandler factory)
    ^
    |
socketHandlers/{domain}Handlers.ts  (each imports utils + external deps)
    ^
    |
socketHandlers/index.ts          (barrel re-export)
    ^
    |
useSocketManager.ts              (imports all handlers, registers with socket)
```

**Direction:** Strictly upward. Handler modules are leaves. The orchestrator is the root.

### sceneSlice.ts sub-slices

```
External deps (AppState type, API client functions, THREE types)
    ^
    |
sceneSlices/{domain}Slice.ts    (each uses StateCreator<AppState>)
    ^
    |
sceneSlice.ts                   (spreads sub-slices, exports SceneSlice type)
    ^
    |
store.tsx                       (unchanged -- still imports createSceneSlice)
```

**Direction:** Sub-slices access each other through `get()` at runtime (Zustand's cross-slice pattern), but have zero import-time dependencies on each other.

## Suggested Split Order

### Order: client.py FIRST, then sceneSlice.ts, then useSocketManager.ts

**Rationale:**

1. **client.py first** because:
   - It has the clearest class boundaries (each class is self-contained).
   - It has 499 backend tests providing instant verification.
   - The `accessors.py` separation is already done -- this extends the same pattern.
   - The re-export `__init__.py` pattern guarantees zero import breakage for tests.
   - No frontend build tooling needed to validate.

2. **sceneSlice.ts second** because:
   - Zustand sub-slice composition is a well-known pattern already used in this codebase.
   - The boundary between geometry/selection/mode/frameEdit is already visible in the interface definition (lines 27-127 of sceneSlice.ts).
   - Only `store.tsx` imports from `sceneSlice.ts` -- minimal blast radius.
   - No runtime behavior change -- just spreading sub-slices instead of inlining everything.

3. **useSocketManager.ts last** because:
   - It is the most complex split: handlers form closures over many dependencies (roomId, queryClient, 20+ store setters).
   - The handler factory pattern requires designing a dependency injection interface.
   - Testing socket handlers requires mocking Socket.IO events, query client, and store -- more test infrastructure than the other two.
   - The current code works fine as-is -- the 949 lines are all in one `useEffect` which is ugly but functional.

### Per-split verification strategy

| Split | Verification |
|-------|-------------|
| client.py | `uv run pytest tests/` (499 tests) + `uv run pyright .` |
| sceneSlice.ts | `bun run build` (TypeScript compilation) + 13 Playwright E2E specs |
| useSocketManager.ts | `bun run build` + 13 Playwright E2E specs |

## Cross-Cutting Concerns

### Shared Utility: `createInvalidateHandler` (useSocketManager)

The current `useSocketManager.ts` defines a `createInvalidateHandler` factory function (lines 85-99) used by two handlers. After the split, this should live in `socketHandlers/utils.ts` and be imported by handler modules that need it.

### Type Imports and Circular References (client.py)

`_socket.py` holds a reference to `ZnDraw` (the class from `_zndraw.py`). This is a forward reference that would be circular at import time. Solution: use `TYPE_CHECKING` guard for the type annotation and accept the `ZnDraw` instance as a runtime parameter (already the pattern in the current `SocketManager.__init__`).

### State Type (`AppState`) Access (sceneSlice.ts)

All sub-slices need the `AppState` type for the `StateCreator` generic parameter. This creates an import from sub-slices back to `store.tsx` which defines `AppState`. This is not circular because `AppState` is a type-only import (erased at runtime). The existing slices already do this.

## Sources

- Zustand slice pattern: existing codebase (`frontend/src/store.tsx`, existing 5 slices)
- Python package re-export pattern: existing codebase (`src/zndraw/geometries/__init__.py`)
- Handler factoring: common React pattern for extracting logic from hooks

---

*Architecture analysis: 2026-03-05*
