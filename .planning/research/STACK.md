# Technology Stack: Monolithic File Refactoring

**Project:** Split Monolithic Files (ZnDraw)
**Researched:** 2026-03-05

## Recommended Stack

This is a structural refactor, not a greenfield build. The "stack" here is the set of techniques, patterns, and tools for safely splitting three monolithic files into well-organized modules. No new libraries need to be added.

### Python: Converting `client.py` (2297 lines) to `client/` Package

| Technique | Purpose | Why |
|-----------|---------|-----|
| Module-to-package conversion | Replace `client.py` with `client/` directory | Python's `__init__.py` re-export mechanism preserves the public API (`from zndraw.client import ZnDraw`) with zero downstream changes. The existing section headers (`# ===`) already mark natural module boundaries. |
| `__init__.py` re-exports with `__all__` | Maintain backwards-compatible import paths | Every consumer (`__init__.py`, `cli.py`, `accessors.py`, `exceptions.py`, `cli_agent/`) imports from `zndraw.client`. Re-exporting from `__init__.py` means zero import changes outside the package. |
| Explicit `as` re-export syntax | Type checker compatibility | Use `from .zndraw import ZnDraw as ZnDraw` (tautological `as`) so pyright and mypy treat it as an intentional re-export, not an unused import. |
| `TYPE_CHECKING` for cross-module refs | Avoid circular imports between new modules | `ZnDraw` references `APIManager` and `SocketManager` at runtime (composition). Extracted modules that need each other's types only for annotations use `if TYPE_CHECKING:` guards. Already used in the codebase. |

**Confidence:** HIGH -- This is the standard Python pattern for module-to-package conversion. The codebase already uses it for `geometries/` and `extensions/`.

#### Proposed Module Structure

```
src/zndraw/client/
  __init__.py          # Re-exports: ZnDraw, APIManager, SocketManager, ZnDrawLock, exceptions
  serialization.py     # atoms_to_json_dict, json_dict_to_atoms, raw_frame_to_atoms, _estimate_frame_size, chunk constants
  exceptions.py        # ZnDrawError, NotConnectedError, RoomLockedError
  lock.py              # ZnDrawLock context manager
  api.py               # APIManager dataclass (lines 254-1137, ~880 lines)
  socket.py            # SocketManager dataclass (lines 1144-1264, ~120 lines)
  zndraw.py            # ZnDraw class (lines 1271-2297, ~1025 lines)
```

**Rationale for these boundaries:**

1. **`serialization.py`** -- Pure functions with no class dependencies. Already marked with a section header. Zero coupling to the rest of the module.
2. **`exceptions.py`** -- Exception classes are leaf nodes in the dependency graph. Other modules raise them but they depend on nothing internal.
3. **`lock.py`** -- `ZnDrawLock` only needs `APIManager` (passed in as constructor arg). Self-contained context manager.
4. **`api.py`** -- `APIManager` is the largest class (880 lines) and is a standalone `@dataclass` that wraps httpx. It depends only on `serialization.py` and `exceptions.py`. The `accessors.py` module already imports it via `TYPE_CHECKING`.
5. **`socket.py`** -- `SocketManager` is small (120 lines), depends on `APIManager` (passed as arg), and handles socket.io event callbacks.
6. **`zndraw.py`** -- The main `ZnDraw` class composes `APIManager`, `SocketManager`, and all accessors. It stays as a single module because its methods are tightly coupled through `self.api` and `self.socket`.

**Why NOT split `ZnDraw` further:** At ~1025 lines, it is large but every method accesses `self.api` or `self.socket`. Splitting it into mixins or partial classes would create artificial coupling. The class implements `MutableSequence[ase.Atoms]` which is a single responsibility (frame collection operations). If it grows further, the `MutableSequence` protocol methods (`__getitem__`, `__setitem__`, `__delitem__`, `insert`, `extend`) could be extracted to a base class, but that is premature now.

### TypeScript: Splitting `sceneSlice.ts` (596 lines) into Dedicated Slices

| Technique | Purpose | Why |
|-----------|---------|-----|
| Zustand slice pattern (existing) | Split one large slice into 3-4 smaller slices | The project already uses `StateCreator<AppState, [], [], SliceInterface>` for 5 slices. Adding more follows the exact same pattern. No new abstractions needed. |
| Cross-slice `get()` access | Allow new slices to read/call each other | Zustand's `get()` returns the full `AppState`, so a `drawingSlice` can call `get().acquireLock()` from `lockSlice`. Already used in `lockSlice` and `sceneSlice`. |
| Re-export helpers from `store.tsx` | Maintain import paths for `getActiveCurves`, `selectPreferredCurve` | Currently re-exported from `store.tsx`. After the split, update re-exports to point to the new slice file. |

**Confidence:** HIGH -- The project already implements this exact pattern for 5 slices. This is adding more of the same.

#### Proposed Slice Structure

```
frontend/src/stores/slices/
  geometrySlice.ts     # geometries, geometrySchemas, geometryDefaults, geometryUpdateSources, geometryFetchingStates, neededFrameKeys, curveRefs, loadedDynamicPositions, particleCount, curveLength
  selectionSlice.ts    # selections, selectionGroups, hoveredGeometryInstance + updateSelections, loadSelectionGroup
  editingSlice.ts      # mode, transformMode, editingSelectedAxis, editingCombinedCentroid, editingCallbacks, pendingFrameEdits, editingFrameDataCount + enter/exitEditingMode, subscribeToEditing, notifyEditingChange, saveFrameEdits
  drawingSlice.ts      # drawingPointerPosition, drawingIsValid, activeCurveForDrawing, attachedCameraKey + enter/exitDrawingMode, attachToCamera
```

**Rationale for these boundaries:**

1. **`geometrySlice`** -- All geometry CRUD state and geometry metadata (schemas, defaults, fetching states, frame keys, curve refs, dynamic positions). These are always accessed together when rendering or updating geometries.
2. **`selectionSlice`** -- Selection state (per-geometry selections, selection groups, hovered instance). Selections are read/written independently of geometry data.
3. **`editingSlice`** -- Editing mode state, transform controls, editing callbacks, and frame edit batching. These form a cohesive "editing session" that's entered/exited together via lock acquisition.
4. **`drawingSlice`** -- Drawing mode state, pointer position, validation, active curve. Drawing is a distinct mode from editing with different UI controls.

**Why this split, not others:** The `mode` field (view/drawing/editing) could live in any slice, but it belongs in `editingSlice` because mode transitions trigger lock acquire/release which is editing-specific logic. The drawing mode is simpler (no lock needed for the mode itself, only for submitting drawn geometry).

### TypeScript: Splitting `useSocketManager.ts` (949 lines) into Handler Modules

| Technique | Purpose | Why |
|-----------|---------|-----|
| Extract handler functions to separate files | Reduce the single useEffect from 949 lines | Handler functions are pure side-effect code that takes store actions and socket/queryClient as dependencies. They can be defined outside the hook. |
| `useEffectEvent` (React 19.2) | Stable callback refs without dependency array bloat | React 19.2 (which this project uses at `^19.2.0`) ships stable `useEffectEvent`. Handlers wrapped in `useEffectEvent` always see latest state without appearing in the `useEffect` dependency array. This eliminates the 30-item dependency array. |
| Factory functions for handlers | Pass dependencies explicitly | Each handler module exports a factory: `createGeometryHandlers(deps)` that returns `{ onGeometryInvalidate, onDefaultCameraInvalidate, onSchemaInvalidate }`. The deps object contains `roomId`, store setters, `queryClient`. |

**Confidence:** HIGH for handler extraction, MEDIUM for `useEffectEvent` adoption.

The `useEffectEvent` recommendation is MEDIUM confidence because while it is stable in React 19.2, it requires careful consideration: the current pattern (massive dependency array with individual selectors) is functional. Switching to `useEffectEvent` is cleaner but changes how the hook re-synchronizes. The extraction of handler functions into separate files can be done without `useEffectEvent` as a first step.

#### Proposed Handler Structure

```
frontend/src/hooks/socket/
  index.ts                    # Re-exports useSocketManager
  useSocketManager.ts         # Slimmed-down orchestrator: connect/disconnect lifecycle, registers handlers from other files
  connectionHandlers.ts       # onConnect (init sequence), onDisconnect, onConnectError
  frameHandlers.ts            # onFrameUpdate, onFrameSelectionUpdate, onFramesInvalidate
  geometryHandlers.ts         # onGeometriesInvalidate, onDefaultCameraInvalidate, onSchemaInvalidate, onFiguresInvalidate
  selectionHandlers.ts        # onSelectionsInvalidate, onSelectionGroupsInvalidate
  chatHandlers.ts             # onChatMessageNew, onChatMessageUpdated, onTyping (+ typingTimeouts Map)
  roomHandlers.ts             # onRoomUpdate, onRoomDelete, onBookmarksInvalidate, onInvalidate
  lockHandlers.ts             # onLockUpdate
  progressHandlers.ts         # onProgressStarted, onProgressUpdate, onProgressComplete
```

**Rationale:** Group by domain concern, not by technical similarity. Chat handlers share the `typingTimeouts` Map. Geometry handlers share the `getGeometry` fetch logic. Frame handlers share the `queryClient` invalidation pattern.

**Why NOT individual hooks per event:** Creating 25 separate `useXxxEvent()` hooks would cause 25 separate `useEffect` calls, each with its own subscribe/unsubscribe lifecycle. This is worse for performance (more effect registrations on mount/unmount) and harder to coordinate (the connect handler must run before event handlers register). The orchestrator-plus-handler-modules pattern keeps one `useEffect` for lifecycle while distributing handler logic.

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| Python module split | `client/` package with `__init__.py` | Keep as single file with section comments | 2297 lines is unmaintainable. Section comments don't provide import boundaries, test isolation, or IDE navigation. |
| Python module split | 6 focused modules | 2 modules (api + everything else) | The serialization helpers and exceptions are leaf dependencies. Extracting them first simplifies the dependency graph for the remaining split. |
| Python refactoring tool | Manual extraction | Rope automated refactoring | Rope's `Move` works for individual functions/classes but doesn't handle module-to-package conversion. The manual approach is straightforward: cut section, paste into new file, add imports. |
| Zustand slice split | 4 new slices from sceneSlice | Keep sceneSlice but organize with comments | Comments don't provide type boundaries. A component subscribing to `geometries` currently re-renders on `drawingPointerPosition` changes. Separate slices fix this via Zustand's selector-based subscriptions. |
| Zustand slice split | 4 slices (geometry, selection, editing, drawing) | 2 slices (geometry+selection, modes) | Selections and geometries have different access patterns. Selections change on user click; geometries change on server push. Separating them allows finer-grained subscriptions. |
| Socket handler split | Handler modules with factory functions | Individual `useXxxSocket` hooks | 25 separate useEffects is worse for performance and coordination. A single orchestrator useEffect with extracted handler factories is the standard pattern. |
| Socket handler split | Plain function extraction | `useEffectEvent` for all handlers | `useEffectEvent` is beneficial but optional. Start with plain extraction (Phase 1), adopt `useEffectEvent` later if the dependency array remains problematic (Phase 2). |
| Socket handler split | Domain-grouped handler files | One handler per event | 25 files for 25 events is over-split. Domain groups (chat, geometry, frames) reflect how these handlers share dependencies and state. |

## Tools and Commands

### Python Refactoring Workflow

```bash
# Verify before each step
uv run pytest tests/ -x --timeout=120

# Type check after each module extraction
uv run pyright .

# Format after changes
uv run ruff format .
uv run ruff check --select I --fix .
```

No new dependencies needed. The refactoring uses only Python's built-in package/module system.

### TypeScript Refactoring Workflow

```bash
# Type check after each slice extraction
cd frontend && bun run tsc --noEmit

# Format after changes
bun run format

# Lint
bun run lint
```

No new dependencies needed. Zustand's `StateCreator` and React's `useEffectEvent` are already available.

## Anti-Patterns to Avoid

### Python

| Anti-Pattern | Why It's Bad | What to Do Instead |
|--------------|-------------|-------------------|
| Mixin classes for `ZnDraw` | Creates diamond inheritance, hides method resolution order, makes `self` types ambiguous | Keep `ZnDraw` as one class. Composition (delegate to `APIManager`) is already used. |
| Lazy imports to avoid circular deps | Masks architectural problems. If two modules import each other at runtime, the abstraction boundary is wrong. | Use `TYPE_CHECKING` for annotation-only imports. Runtime dependencies should flow one direction: `zndraw.py` -> `api.py` -> `serialization.py`. |
| `__init__.py` that imports everything | Defeats the purpose of splitting by loading all modules on first import | Only re-export the public API: `ZnDraw`, `APIManager`, `SocketManager`, `ZnDrawLock`, exceptions. Internal helpers stay internal. |
| Keeping backwards-compat shims | PROJECT.md explicitly says "clean break allowed on imports" | Update all internal imports directly. No `warnings.warn("import from X instead")` shims. |

### TypeScript / Zustand

| Anti-Pattern | Why It's Bad | What to Do Instead |
|--------------|-------------|-------------------|
| Separate Zustand stores instead of slices | Breaks cross-slice `get()` access. `editingSlice` needs `lockSlice.acquireLock()`. Separate stores can't share `get()`. | Keep the single `create<AppState>()` store. Add slices, not stores. |
| `immer` middleware on new slices | Inconsistent with existing slices that use plain `set()`. Mixing patterns in one store confuses contributors. | Use plain `set()` like existing `connectionSlice`, `playbackSlice`, `uiSlice`. Only domain-specific stores like `geometryStore` use immer. |
| Moving `mode` to a separate `modeSlice` | The mode transitions (`enterEditingMode`, `exitDrawingMode`) touch lock state, geometry state, and selection state. Isolating `mode` creates a slice that must reach into 3 other slices for every action. | Keep `mode` in `editingSlice` since editing mode transitions are the most complex (lock acquire/release). Drawing mode actions live in `drawingSlice` but set `mode` via `set()` on the shared state. |
| Wrapping every handler in `useEffectEvent` immediately | Premature optimization. The current pattern works. Introducing `useEffectEvent` everywhere in one PR adds risk. | Extract handler functions first (pure refactor). Consider `useEffectEvent` as a follow-up optimization if the dependency array causes issues. |

## Verification Strategy

### Python: Confirming the Split Preserves Behavior

1. **Import compatibility test:** After creating `client/`, verify that `from zndraw.client import ZnDraw` still works. Run `uv run python -c "from zndraw import ZnDraw; print(ZnDraw)"`.
2. **Full test suite:** `uv run pytest tests/ -x` -- all 499 tests must pass. The tests exercise the client through the test server, so any broken import or missing method surfaces immediately.
3. **Type checking:** `uv run pyright .` -- catches missing re-exports, incorrect `TYPE_CHECKING` guards, and broken cross-module references.

### TypeScript: Confirming the Split Preserves Behavior

1. **Type checking:** `bun run tsc --noEmit` -- catches missing exports, incorrect `StateCreator` type params, and broken cross-slice `get()` calls.
2. **Build:** `bun run build` -- confirms Vite can resolve all new module paths.
3. **E2E tests:** The 13 Playwright specs exercise geometry editing, selections, and drawing modes, which are exactly the concerns being split.

## Sources

- [Python `__init__.py` re-exports](https://realpython.com/python-init-py/) -- MEDIUM confidence (community source)
- [Scientific Python exports guide](https://learn.scientific-python.org/development/patterns/exports/) -- HIGH confidence (official scientific Python guide, explains `as` re-export pattern)
- [Zustand slices pattern (DeepWiki)](https://deepwiki.com/pmndrs/zustand/7.1-slices-pattern) -- MEDIUM confidence (derived from official docs)
- [Zustand slices wiki](https://github.com/pmndrs/zustand/wiki/Splitting-the-store-into-separate-slices) -- HIGH confidence (official repo wiki)
- [React 19.2 useEffectEvent](https://react.dev/reference/react/useEffectEvent) -- HIGH confidence (official React docs)
- [React 19.2 release blog](https://react.dev/blog/2025/10/01/react-19-2) -- HIGH confidence (official React blog)
- [Rope refactoring library](https://rope.readthedocs.io/en/latest/overview.html) -- HIGH confidence (official docs, evaluated and rejected for this use case)
- Existing codebase patterns: `geometries/`, `extensions/`, `connectionSlice.ts`, `lockSlice.ts` -- highest confidence (verified in-repo)
