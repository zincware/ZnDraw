# Phase 2: Scene Sub-Slices - Context

**Gathered:** 2026-03-06
**Status:** Ready for planning

<domain>
## Phase Boundary

Decompose `sceneSlice.ts` (596 lines) into focused sub-slices (geometry, selection, editing, drawing) that compose back into the existing `SceneSlice` interface without changing `store.tsx` or any consuming components. Pure structural refactor — no functional changes. All 13 Playwright E2E specs must pass unchanged.

</domain>

<decisions>
## Implementation Decisions

### Camera placement
- Camera state (`attachedCameraKey`, `attachToCamera`, `setAttachedCameraKey`) merges into the geometry sub-slice
- Cameras are a geometry type — `attachToCamera` already reads from the `geometries` dict
- Keeps the sub-slice count at the planned 4 (geometry, selection, editing, drawing)

### Cross-cutting actions
- `removeGeometry` stays in the geometry sub-slice — it's fundamentally a geometry operation
- Cross-domain cleanup (selections, activeCurve, camera, editingCallbacks) remains atomic in a single `set()` call, accessing other slice state via `state` (since all sub-slices share `AppState`)

### Mode transitions
- `enterDrawingMode` / `exitDrawingMode` live in the drawing sub-slice
- `enterEditingMode` / `exitEditingMode` live in the editing sub-slice
- Each action "owns" its mode transition and uses `get()` to access cross-slice state (lock, geometries, selections)

### Mode state ownership
- The `mode` field (`"view" | "drawing" | "editing"`) and `setMode` setter live in the geometry sub-slice
- Geometry acts as the "base" sub-slice; `mode` is scene-level state not exclusive to drawing or editing
- Drawing and editing sub-slices read/write `mode` via `get()` / `set()`

### Claude's Discretion
- Directory structure for sub-slice files (new `scene/` subdirectory or flat alongside existing slices)
- Placement of helper functions (`getActiveCurves`, `selectPreferredCurve`) and re-export adjustments in `store.tsx`
- Exact field-to-sub-slice mapping for unambiguous members (e.g., `curveRefs` → geometry, `drawingPointerPosition` → drawing)
- TypeScript composition mechanics (`StateCreator` generics for sub-slices)

</decisions>

<code_context>
## Existing Code Insights

### Reusable Assets
- `store.tsx`: Spread-based composition pattern `...createSceneSlice(...a)` — sub-slices must follow this same pattern
- `store.tsx` re-exports `getActiveCurves`, `selectPreferredCurve` from sceneSlice — re-export path needs updating
- Other slices (`connectionSlice.ts`, `playbackSlice.ts`, `lockSlice.ts`, `uiSlice.ts`) are single files in `stores/slices/` — establishes the existing flat structure

### Established Patterns
- `StateCreator<AppState, [], [], SliceType>` generic signature used by all slice creators
- All sub-slices will share `AppState` via Zustand's `set`/`get` — cross-slice `get()` calls are idiomatic
- `SceneSlice` interface has ~60 members across geometry (~25), selection (~8), editing (~15), drawing (~7), camera (~3), and mode (~2) concerns

### Integration Points
- `store.tsx` line 8-9: imports `SceneSlice` type and `createSceneSlice` — must keep working or update to composite import
- `store.tsx` line 14-17: re-exports `getActiveCurves`, `selectPreferredCurve` — source path changes
- `frontend/src/myapi/client`: provides `acquireEditLock`, `releaseEditLock`, `createGeometry`, `getGeometry`, `partialUpdateFrame`, `updateActiveCamera`, `updateSelection` — used by various sub-slices
- Components consuming `useAppStore` selectors for scene state — no changes needed since `AppState` type composition is preserved

</code_context>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 02-scene-sub-slices*
*Context gathered: 2026-03-06*
