# Phase 2: Scene Sub-Slices - Research

**Researched:** 2026-03-06
**Domain:** Zustand v5 slice decomposition / TypeScript interface composition
**Confidence:** HIGH

## Summary

Phase 2 decomposes the monolithic `sceneSlice.ts` (596 lines, ~60 interface members) into four focused sub-slices: geometry, selection, editing, and drawing. The project already uses Zustand v5.0.8's `StateCreator<AppState, [], [], SliceType>` slice pattern throughout all five existing slices. The sub-slice decomposition follows the same pattern -- each sub-slice is a `StateCreator` function whose output is a sub-interface, and `createSceneSlice` becomes a thin composer that spreads all four sub-slice creators.

The key insight is that this is purely a structural refactor. The `AppState` type in `store.tsx` remains `ConnectionSlice & PlaybackSlice & SceneSlice & LockSlice & UISlice` -- only the internal composition of `SceneSlice` changes. Since all sub-slices receive `AppState` as their first generic parameter, cross-slice access via `get()` and `set()` continues to work identically. No component selectors change because `useAppStore` still exposes the same flat `AppState` shape.

**Primary recommendation:** Create a `stores/slices/scene/` directory with `geometrySubSlice.ts`, `selectionSubSlice.ts`, `editingSubSlice.ts`, `drawingSubSlice.ts`, and an `index.ts` barrel that re-exports the composed `SceneSlice` type and `createSceneSlice` function. This preserves the existing import path from `store.tsx`.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Camera state (`attachedCameraKey`, `attachToCamera`, `setAttachedCameraKey`) merges into the geometry sub-slice
- Cameras are a geometry type -- `attachToCamera` already reads from the `geometries` dict
- Keeps the sub-slice count at the planned 4 (geometry, selection, editing, drawing)
- `removeGeometry` stays in the geometry sub-slice -- it's fundamentally a geometry operation
- Cross-domain cleanup (selections, activeCurve, camera, editingCallbacks) remains atomic in a single `set()` call, accessing other slice state via `state` (since all sub-slices share `AppState`)
- `enterDrawingMode` / `exitDrawingMode` live in the drawing sub-slice
- `enterEditingMode` / `exitEditingMode` live in the editing sub-slice
- Each action "owns" its mode transition and uses `get()` to access cross-slice state (lock, geometries, selections)
- The `mode` field (`"view" | "drawing" | "editing"`) and `setMode` setter live in the geometry sub-slice
- Geometry acts as the "base" sub-slice; `mode` is scene-level state not exclusive to drawing or editing
- Drawing and editing sub-slices read/write `mode` via `get()` / `set()`

### Claude's Discretion
- Directory structure for sub-slice files (new `scene/` subdirectory or flat alongside existing slices)
- Placement of helper functions (`getActiveCurves`, `selectPreferredCurve`) and re-export adjustments in `store.tsx`
- Exact field-to-sub-slice mapping for unambiguous members (e.g., `curveRefs` -> geometry, `drawingPointerPosition` -> drawing)
- TypeScript composition mechanics (`StateCreator` generics for sub-slices)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| SCEN-01 | Geometry sub-slice extracted (`geometrySlice.ts`) | Member mapping complete (25 state + action members); StateCreator pattern verified |
| SCEN-02 | Selection sub-slice extracted (`selectionSlice.ts`) | Member mapping complete (8 members); cross-slice `roomId` access via `get()` verified |
| SCEN-03 | Editing sub-slice extracted (`editingSlice.ts`) | Member mapping complete (15 members); lock + snackbar cross-slice access verified |
| SCEN-04 | Drawing sub-slice extracted (`drawingSlice.ts`) | Member mapping complete (7 members); curve creation + lock cross-slice access verified |
| SCEN-05 | Camera concern handled (merged into geometry) | Locked decision: 3 camera members merge into geometry sub-slice |
| SCEN-06 | `SceneSlice` interface composed from sub-interfaces | TypeScript intersection type `GeometrySubSlice & SelectionSubSlice & EditingSubSlice & DrawingSubSlice` verified |
| SCEN-07 | `createSceneSlice` composes sub-creators following existing pattern | Spread-based composition `...createGeometrySubSlice(set, get, store)` verified against existing `store.tsx` pattern |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| zustand | 5.0.8 | State management with slices pattern | Already installed; `StateCreator<T, Mis, Mos, U>` is the composition primitive |
| typescript | 5.9.3 | Type-safe interface composition | Already installed; intersection types compose sub-slice interfaces |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| @playwright/test | 1.58.0 | E2E regression testing | Run existing 12 E2E specs to verify no behavioral changes |
| biome | (devDep) | Formatting and linting | `bun run format` and `bun run lint` after changes |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Sub-directory `scene/` | Flat files in `slices/` | Sub-directory keeps 4 new files grouped; flat would clutter `slices/` which currently has 5 files |

## Architecture Patterns

### Recommended Project Structure
```
frontend/src/stores/slices/
├── connectionSlice.ts       # existing (unchanged)
├── lockSlice.ts             # existing (unchanged)
├── playbackSlice.ts         # existing (unchanged)
├── uiSlice.ts               # existing (unchanged)
├── sceneSlice.ts            # DELETED after migration
└── scene/
    ├── index.ts             # barrel: exports SceneSlice, createSceneSlice, helpers
    ├── geometrySubSlice.ts  # GeometrySubSlice interface + createGeometrySubSlice
    ├── selectionSubSlice.ts # SelectionSubSlice interface + createSelectionSubSlice
    ├── editingSubSlice.ts   # EditingSubSlice interface + createEditingSubSlice
    └── drawingSubSlice.ts   # DrawingSubSlice interface + createDrawingSubSlice
```

### Pattern 1: Sub-Slice StateCreator Signature
**What:** Each sub-slice uses the same `StateCreator<AppState, [], [], SubSliceType>` signature as top-level slices, because Zustand's `StateCreator` first generic is always the full store type.
**When to use:** Every sub-slice file.
**Example:**
```typescript
// Source: existing pattern in connectionSlice.ts, lockSlice.ts, etc.
import type { StateCreator } from "zustand";
import type { AppState } from "../../../store";

export interface GeometrySubSlice {
  geometries: Record<string, any>;
  // ... other geometry members
  setGeometries: (geometries: Record<string, any>) => void;
  // ... other geometry actions
}

export const createGeometrySubSlice: StateCreator<
  AppState,
  [],
  [],
  GeometrySubSlice
> = (set, get) => ({
  geometries: {},
  // ... defaults and actions
});
```

### Pattern 2: Barrel Composition in scene/index.ts
**What:** The barrel file composes sub-slice interfaces via intersection and spreads sub-slice creators, then re-exports the combined type and creator using the original names (`SceneSlice`, `createSceneSlice`).
**When to use:** The single `scene/index.ts` file.
**Example:**
```typescript
// scene/index.ts
import type { StateCreator } from "zustand";
import type { AppState } from "../../../store";
import type { GeometrySubSlice } from "./geometrySubSlice";
import { createGeometrySubSlice } from "./geometrySubSlice";
import type { SelectionSubSlice } from "./selectionSubSlice";
import { createSelectionSubSlice } from "./selectionSubSlice";
import type { EditingSubSlice } from "./editingSubSlice";
import { createEditingSubSlice } from "./editingSubSlice";
import type { DrawingSubSlice } from "./drawingSubSlice";
import { createDrawingSubSlice } from "./drawingSubSlice";

// Composed type -- identical shape to original SceneSlice
export type SceneSlice = GeometrySubSlice &
  SelectionSubSlice &
  EditingSubSlice &
  DrawingSubSlice;

// Composed creator -- same spread pattern as store.tsx
export const createSceneSlice: StateCreator<AppState, [], [], SceneSlice> = (
  set,
  get,
  store,
) => ({
  ...createGeometrySubSlice(set, get, store),
  ...createSelectionSubSlice(set, get, store),
  ...createEditingSubSlice(set, get, store),
  ...createDrawingSubSlice(set, get, store),
});

// Re-export helpers
export { getActiveCurves, selectPreferredCurve } from "./geometrySubSlice";
```

### Pattern 3: Cross-Slice Access via get()/set()
**What:** Sub-slices access state owned by other sub-slices (or other top-level slices) through `get()` and `set()`, which always operate on the full `AppState`. This is the established Zustand pattern.
**When to use:** Any action that needs cross-domain state (e.g., `removeGeometry` clears selections, `enterEditingMode` calls `acquireLock`).
**Example:**
```typescript
// In geometrySubSlice.ts -- removeGeometry needs selection + editing state
removeGeometry: (key) =>
  set((state) => {
    const { [key]: removed, ...rest } = state.geometries;
    const { [key]: removedSelection, ...restSelections } = state.selections;
    // state includes ALL AppState members -- cross-slice access works
    return {
      geometries: rest,
      selections: restSelections,
      // ... other cross-domain cleanup
    };
  }),
```

### Pattern 4: Import Path Update in store.tsx
**What:** `store.tsx` currently imports from `./stores/slices/sceneSlice`. After refactoring, it imports from `./stores/slices/scene` (the barrel). The `AppState` type composition and `create()` call remain unchanged.
**When to use:** Single change in `store.tsx`.
**Example:**
```typescript
// store.tsx -- only import path changes
import type { SceneSlice } from "./stores/slices/scene";
import { createSceneSlice } from "./stores/slices/scene";

// Re-exports -- only source path changes
export {
  getActiveCurves,
  selectPreferredCurve,
} from "./stores/slices/scene";

// AppState type -- UNCHANGED
export type AppState = ConnectionSlice &
  PlaybackSlice &
  SceneSlice &  // still the same composed type
  LockSlice &
  UISlice;

// Store creation -- UNCHANGED
export const useAppStore = create<AppState>((...a) => ({
  ...createConnectionSlice(...a),
  ...createPlaybackSlice(...a),
  ...createSceneSlice(...a),  // internally spreads 4 sub-slices
  ...createLockSlice(...a),
  ...createUISlice(...a),
}));
```

### Anti-Patterns to Avoid
- **Importing sub-slices directly from components:** Components should only import from `store.tsx` (`useAppStore`). Sub-slices are internal to the store composition. If a component needs `GeometrySubSlice` type, expose it through the barrel, but selectors should use `AppState`.
- **Duplicating state across sub-slices:** Each field lives in exactly one sub-slice. Cross-slice access is via `get()`/`set()`, never by copying state.
- **Using separate Zustand stores for sub-slices:** Out of scope per REQUIREMENTS.md ("Converting Zustand slices to separate stores would break `AppState` composition and cross-slice `get()` calls").

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Type composition | Manual union type with all 60 members | `GeometrySubSlice & SelectionSubSlice & EditingSubSlice & DrawingSubSlice` | TypeScript intersection types are the idiomatic composition mechanism |
| Slice composition | Custom merge function | Spread operator `...createSubSlice(set, get, store)` | Matches existing `store.tsx` pattern exactly; Zustand docs recommend this |
| Cross-slice state access | Passing state between sub-slices as parameters | `get()` / `set()` on full `AppState` | Zustand's design -- `set`/`get` always operate on the complete store |

**Key insight:** This refactor requires zero new libraries or utilities. The existing Zustand `StateCreator` type and spread composition pattern handle everything. The work is purely taxonomic (classifying members) and mechanical (moving code).

## Member Classification

### GeometrySubSlice (~25 members)
**State fields:**
- `geometries`, `geometrySchemas`, `geometryDefaults`, `geometryUpdateSources`, `geometryFetchingStates`, `neededFrameKeys`
- `mode` (locked decision: geometry is the "base" sub-slice)
- `attachedCameraKey` (locked decision: camera merges into geometry)
- `curveRefs`
- `hoveredGeometryInstance`
- `particleCount`, `curveLength`
- `loadedDynamicPositions`

**Actions:**
- `setGeometries`, `setGeometrySchemas`, `setGeometryDefaults`, `updateGeometry`, `removeGeometry`
- `setGeometryFetching`, `removeGeometryFetching`, `registerFrameKeys`, `unregisterFrameKeys`, `getIsFetching`
- `setMode` (locked decision)
- `attachToCamera`, `setAttachedCameraKey` (locked decision)
- `registerCurveRef`, `unregisterCurveRef`
- `setHoveredGeometryInstance`
- `setParticleCount`, `setCurveLength`
- `registerLoadedDynamicPositions`, `unregisterLoadedDynamicPositions`

**Helpers (module-level functions):**
- `getActiveCurves`, `selectPreferredCurve`

**API imports needed:** `createGeometry`, `getGeometry`, `updateActiveCamera`

### SelectionSubSlice (~8 members)
**State fields:**
- `selections`, `selectionGroups`

**Actions:**
- `setSelections`, `updateSelectionForGeometry`, `setSelectionGroups`, `loadSelectionGroup`, `updateSelections`

**Cross-slice access:** Reads `roomId` from `ConnectionSlice` via `get()`

**API imports needed:** `updateSelection` (aliased as `updateSelectionAPI`)

### EditingSubSlice (~15 members)
**State fields:**
- `transformMode`, `editingSelectedAxis`, `editingCombinedCentroid`, `editingCallbacks`
- `pendingFrameEdits`, `editingFrameDataCount`

**Actions:**
- `enterEditingMode`, `exitEditingMode` (locked decision)
- `setTransformMode`, `cycleTransformMode`, `setEditingSelectedAxis`
- `setEditingCombinedCentroid`
- `subscribeToEditing`, `notifyEditingChange`
- `setPendingFrameEdit`, `clearPendingFrameEdits`, `saveFrameEdits`
- `incrementEditingFrameDataCount`, `decrementEditingFrameDataCount`

**Cross-slice access:** Reads `roomId`, `mode`, `selections` via `get()`; calls `acquireLock`, `releaseLock`, `showSnackbar`, `updateSelectionForGeometry`, `saveFrameEdits` via `get()`

**API imports needed:** `partialUpdateFrame`

### DrawingSubSlice (~7 members)
**State fields:**
- `drawingPointerPosition`, `drawingIsValid`, `activeCurveForDrawing`

**Actions:**
- `enterDrawingMode`, `exitDrawingMode` (locked decision)
- `setDrawingPointerPosition`, `setDrawingIsValid`, `setActiveCurveForDrawing`

**Cross-slice access:** Reads `geometries`, `mode`, `roomId` via `get()`; calls `acquireLock`, `releaseLock` via `get()`

**API imports needed:** `createGeometry`, `getGeometry`

## Common Pitfalls

### Pitfall 1: Circular Import Between Sub-Slices and store.tsx
**What goes wrong:** Sub-slice files import `AppState` from `../../store.tsx`, which imports `SceneSlice` from the barrel. If the barrel re-exports types from sub-slices that import `AppState`, a circular dependency forms.
**Why it happens:** TypeScript can handle circular type-only imports (since types are erased at runtime), but it can cause confusing tooling issues.
**How to avoid:** All sub-slice files import `type { AppState }` (type-only import). TypeScript erases these at compile time, breaking the runtime cycle. The existing slices already use this pattern (`import type { AppState } from "../../store"`).
**Warning signs:** IDE red squiggles on `AppState` import, or `undefined` errors at runtime.

### Pitfall 2: Forgetting to Pass `store` as Third Argument
**What goes wrong:** The `StateCreator` function receives `(set, get, store)` -- three arguments. If the barrel's `createSceneSlice` only forwards `(set, get)`, sub-slices that might use `store` (for subscriptions) will fail.
**Why it happens:** Most existing slice creators only destructure `(set, get)` and never use `store`, so it's easy to forget.
**How to avoid:** Always forward all three arguments: `...createGeometrySubSlice(set, get, store)`.
**Warning signs:** TypeScript error on the `StateCreator` return type, or runtime error in `store.subscribe()`.

### Pitfall 3: SceneSlice Type Mismatch After Decomposition
**What goes wrong:** The composed `GeometrySubSlice & SelectionSubSlice & EditingSubSlice & DrawingSubSlice` type must be structurally identical to the original `SceneSlice`. If any member is missing or has a different signature, TypeScript will error on the `AppState` type.
**Why it happens:** Member accidentally omitted from one sub-slice, or a typo in a member name.
**How to avoid:** Delete the old `SceneSlice` interface last. First, create all sub-slices and verify `tsc --noEmit` passes. The compiler enforces structural compatibility since `AppState` uses the composed type.
**Warning signs:** TypeScript errors in components that use `useAppStore` selectors for scene state.

### Pitfall 4: Helper Functions Losing Access After Move
**What goes wrong:** `getActiveCurves` and `selectPreferredCurve` are used inside `removeGeometry` (geometry sub-slice) and `enterDrawingMode` (drawing sub-slice), plus re-exported from `store.tsx` and used in `Canvas.tsx`.
**Why it happens:** When helpers move to `geometrySubSlice.ts`, the drawing sub-slice needs to import them. The re-export chain in `store.tsx` also needs updating.
**How to avoid:** Keep helpers in `geometrySubSlice.ts` (they operate on geometry data). Drawing sub-slice imports them from `./geometrySubSlice`. The barrel re-exports them. `store.tsx` updates its re-export source.
**Warning signs:** Compilation errors in `drawingSubSlice.ts` or `Canvas.tsx`.

### Pitfall 5: Mutating Shared References (Map/Set) Across Sub-Slices
**What goes wrong:** `editingCallbacks` is a `Map<string, Set<...>>`. The `removeGeometry` action (geometry sub-slice) deletes entries from it. If the geometry sub-slice creates a new Map copy but the editing sub-slice still holds the old reference, state becomes inconsistent.
**Why it happens:** Maps and Sets are reference types; Zustand's shallow merge won't catch reference equality issues.
**How to avoid:** `removeGeometry` already creates a new Map via `new Map(state.editingCallbacks)` before deleting. This pattern must be preserved exactly. The key insight: `removeGeometry` accesses `editingCallbacks` via the `state` parameter in `set((state) => ...)`, which always reflects current state.
**Warning signs:** Editing callbacks not being cleaned up when geometries are removed.

## Code Examples

### Complete Sub-Slice File Template
```typescript
// Source: derived from existing lockSlice.ts and sceneSlice.ts patterns
import type { StateCreator } from "zustand";
import type { AppState } from "../../../store";

export interface ExampleSubSlice {
  // State
  exampleField: string;
  // Actions
  setExampleField: (value: string) => void;
  crossSliceAction: () => void;
}

export const createExampleSubSlice: StateCreator<
  AppState,
  [],
  [],
  ExampleSubSlice
> = (set, get) => ({
  exampleField: "default",

  setExampleField: (value) => set({ exampleField: value }),

  crossSliceAction: () => {
    // Access other slices via get()
    const { roomId } = get();
    if (!roomId) return;
    // Modify own state + other state via set()
    set({ exampleField: "updated" });
  },
});
```

### Barrel File with Type Composition
```typescript
// Source: derived from store.tsx pattern
import type { StateCreator } from "zustand";
import type { AppState } from "../../../store";

import type { DrawingSubSlice } from "./drawingSubSlice";
import { createDrawingSubSlice } from "./drawingSubSlice";
import type { EditingSubSlice } from "./editingSubSlice";
import { createEditingSubSlice } from "./editingSubSlice";
import type { GeometrySubSlice } from "./geometrySubSlice";
import { createGeometrySubSlice } from "./geometrySubSlice";
import type { SelectionSubSlice } from "./selectionSubSlice";
import { createSelectionSubSlice } from "./selectionSubSlice";

export type SceneSlice = GeometrySubSlice &
  SelectionSubSlice &
  EditingSubSlice &
  DrawingSubSlice;

export const createSceneSlice: StateCreator<AppState, [], [], SceneSlice> = (
  set,
  get,
  store,
) => ({
  ...createGeometrySubSlice(set, get, store),
  ...createSelectionSubSlice(set, get, store),
  ...createEditingSubSlice(set, get, store),
  ...createDrawingSubSlice(set, get, store),
});

export { getActiveCurves, selectPreferredCurve } from "./geometrySubSlice";
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Single monolithic slice file | Sub-slice decomposition via barrel pattern | N/A -- project-specific refactor | Files go from 596 lines to ~100-200 each |
| Zustand v4 `StateCreator` | Zustand v5 `StateCreator` (same API, stricter types) | v5.0.0 (2024) | No functional change; types are compatible |
| `create()(...)` double-call for TS | `create<T>(...)` direct call | v5.0.0 | Project already uses v5 single-call pattern |

**Deprecated/outdated:**
- Zustand v4's `create()()` double-parentheses pattern is not used in this project. v5's `create<AppState>(...)` is already in place.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Playwright 1.58.0 |
| Config file | `frontend/playwright.config.ts` |
| Quick run command | `cd frontend && bunx playwright test --project=chromium` |
| Full suite command | `cd frontend && bunx playwright test` |

Note: There is no Vitest/Jest setup for frontend unit tests (deferred to v2 per TEST-01). Validation relies on TypeScript compilation (`tsc --noEmit`) and the 12 Playwright E2E specs.

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SCEN-01 | Geometry sub-slice extracted | compile | `cd frontend && bunx tsc --noEmit` | N/A (type check) |
| SCEN-02 | Selection sub-slice extracted | compile | `cd frontend && bunx tsc --noEmit` | N/A (type check) |
| SCEN-03 | Editing sub-slice extracted | compile + E2E | `cd frontend && bunx playwright test e2e/editing.spec.ts` | editing.spec.ts exists |
| SCEN-04 | Drawing sub-slice extracted | compile + E2E | `cd frontend && bunx playwright test e2e/geometry-drawing.spec.ts` | geometry-drawing.spec.ts exists |
| SCEN-05 | Camera concern handled | compile + E2E | `cd frontend && bunx playwright test e2e/camera-session.spec.ts` | camera-session.spec.ts exists |
| SCEN-06 | SceneSlice composed from sub-interfaces | compile | `cd frontend && bunx tsc --noEmit` | N/A (type check) |
| SCEN-07 | createSceneSlice composes sub-creators | compile + E2E | `cd frontend && bunx playwright test` | All 12 specs |

### Sampling Rate
- **Per task commit:** `cd frontend && bunx tsc --noEmit`
- **Per wave merge:** `cd frontend && bunx playwright test`
- **Phase gate:** Full E2E suite green + `tsc --noEmit` clean

### Wave 0 Gaps
None -- existing test infrastructure (TypeScript compiler + 12 Playwright E2E specs) covers all phase requirements. No new test files needed for this structural refactor.

## Open Questions

1. **E2E spec count: 12 vs 13**
   - What we know: CONTEXT.md and requirements mention "13 Playwright E2E specs" but the filesystem contains 12 `.spec.ts` files in `frontend/e2e/`.
   - What's unclear: Whether a 13th spec was planned but not yet created, or was removed.
   - Recommendation: Count passes by running `bunx playwright test` and verify all 12 pass. The number is not critical -- the gate is "all existing E2E specs pass unchanged."

2. **`docs-screenshots.spec.ts` relevance**
   - What we know: This spec takes documentation screenshots. It may have external dependencies (running server, specific room state).
   - What's unclear: Whether it's included in regular CI runs or is manual-only.
   - Recommendation: Include in the full suite run but note it may be skipped/flaky in local dev. Focus on the 11 functional E2E specs for regression.

## Sources

### Primary (HIGH confidence)
- `frontend/src/stores/slices/sceneSlice.ts` -- 596 lines, complete SceneSlice interface with 60 members, reviewed in full
- `frontend/src/store.tsx` -- AppState composition, create() call, re-exports, reviewed in full
- `frontend/src/stores/slices/connectionSlice.ts` -- established StateCreator pattern, reviewed in full
- `frontend/src/stores/slices/lockSlice.ts` -- cross-slice access pattern (mode, showSnackbar), reviewed in full
- `frontend/node_modules/zustand/vanilla.d.ts` -- `StateCreator<T, Mis, Mos, U>` type definition verified from installed v5.0.8
- `frontend/package.json` -- zustand@5.0.8, typescript@5.9.3, @playwright/test@1.58.0

### Secondary (MEDIUM confidence)
- [Zustand StateCreator Type Discussion](https://github.com/pmndrs/zustand/discussions/2571) -- confirmed 4 generic parameters and cross-slice access design
- [Zustand Slices Pattern (DeepWiki)](https://deepwiki.com/pmndrs/zustand/7.1-slices-pattern) -- spread-based composition is the documented pattern

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all libraries already installed and patterns established in codebase
- Architecture: HIGH -- sub-slice composition follows the exact same pattern as existing top-level slice composition; verified against installed Zustand types
- Member classification: HIGH -- every member of SceneSlice systematically classified; locked decisions from CONTEXT.md applied
- Pitfalls: HIGH -- all pitfalls derived from actual code analysis (Map/Set references, import chains, helper function dependencies)

**Research date:** 2026-03-06
**Valid until:** 2026-04-06 (stable -- Zustand v5 API is mature, no breaking changes expected)
