# Dockview UI Redesign — Design Spec

**Date:** 2026-04-02
**Status:** Draft (revised 2026-04-16)

## Overview

Replace ZnDraw's fixed sidebar + floating window layout with a VS Code-style two-tier architecture using `dockview-react`:

- **Tools** (Selections, Modifiers, Analysis, Geometries, Plots Browser, Chat) live in single-panel-per-bar zones controlled by activity bar icons. Three activity bars exist (Left, Right, Bottom). Icons can be dragged between bars.
- **Views** (3D Viewer, plot tabs) live in a splittable editor area in the center. Tabs can be split horizontally/vertically, dragged between groups, floated over the grid, or popped out to a separate browser window.

The design follows VS Code's well-established conventions: tools configure, views show content, and the two tiers have different layout primitives.

## Goals

- VS Code-like docking experience with three activity bars (Left/Right/Bottom) and a splittable editor area
- Activity bar icons can be dragged between bars (Left ↔ Right ↔ Bottom)
- Each bar shows at most one tool panel at a time (single-panel-per-bar behavior)
- Editor area supports splitting, tabbing, floating-over-grid, and pop-out windows for views
- Plots Browser panel lists all available plots; clicking opens them as tabs in the editor area
- 3D Viewer is pinned in the editor area — always present, cannot be closed
- New plots from the server auto-open as tabs in the active editor group
- Floating chat (current ChatGPT-style overlay) preserved via dockview's native floating groups
- Default startup matches the current UI feel: 3D viewer fills everything, sidebars empty

## Non-Goals (YAGNI)

- Layout persistence (localStorage / server) — deferred. dockview supports `toJSON()`/`fromJSON()`; cheap to add later.
- Multiple 3D viewports — architecturally supported but not implemented (no camera management for multiple viewports).
- VS Code native multi-panel integration — explicitly out of scope. The existing `vscode-zndraw` single-webview approach continues to work as-is. A separate spec may later cover postMessage-based VS Code enhancements.
- Bar icon ordering persistence — out of scope.
- Custom dockview themes beyond a light/dark mapping to MUI.
- Keyboard shortcuts for panel focus.
- Drag plots from the browser panel into the editor area as splits — for now, clicking opens as a tab in the active group.

## Backend Changes

**None.** The existing REST + Socket.IO surface covers everything this redesign needs:

| Panel | Existing surface |
|-------|------------------|
| Plots Browser | `GET /v1/rooms/{room_id}/figures` (list) + `GET /v1/rooms/{room_id}/figures/{key}` (data) + `figure_invalidate` socket event |
| Selections | existing selection-groups + frame-selection routes + `SelectionGroupsInvalidate` |
| Modifiers / Analysis | existing extension job system |
| Geometries | existing CRUD + `GeometryInvalidate` |
| Chat | existing message routes + `MessageNew` / `MessageEdited` |

This is a frontend-only refactor.

## Library

**dockview-react** (MIT, zero dependencies, daily maintenance, currently v5.2.0)

Used capabilities:

- Splittable grid + tab groups for the editor area
- `floatingGroupBounds="boundedWithinViewport"` for clamping
- `addPanel({ floating: true })` and `addFloatingGroup(group, { position })` for in-app floating panels
- Shift+drag-tab gesture for users to float a tab
- `addPopoutGroup(group, ...)` for separate browser windows
- Serializable layout state via `toJSON()`/`fromJSON()`
- React 19 compatible, TypeScript core

## Architecture

### Two Tiers

**Tier 1 — Tools (managed outside dockview):**

The three activity bars (Left, Right, Bottom) and their associated panel zones are managed by a Zustand slice (not dockview). This is because:

1. Each bar holds at most one panel at a time, not a splittable grid
2. Bar icons can be dragged between bars via HTML5 drag-and-drop
3. The interaction model differs fundamentally from dockview's grid

State shape:

```ts
type ActivityBarSlice = {
  leftBarIcons: PanelId[];       // ordered list
  rightBarIcons: PanelId[];
  bottomBarIcons: PanelId[];
  activeLeft: PanelId | null;    // currently visible panel in left zone
  activeRight: PanelId | null;
  activeBottom: PanelId | null;
};
```

Default state: all six tool icons on `leftBarIcons`, `activeLeft = null`, right/bottom empty.

**Tier 2 — Views (managed by dockview):**

A single `DockviewReact` instance manages the editor area between the AppBar and the FrameProgressBar. The 3D Viewer and any opened plot views live inside it.

### Root Layout

```
<App>
  <LandingPage>
    <AppBar />                         — outside dockview, always visible
    <Box class="workspace">            — flex row
      <ActivityBar position="left" />  — outside dockview
      <SidebarZone position="left" />  — outside dockview, shows activeLeft panel
      <DockviewReact />                — editor area (3D Viewer + plot views)
      <SidebarZone position="right" /> — outside dockview, shows activeRight panel
      <ActivityBar position="right" /> — outside dockview
    </Box>
    <BottomZone />                     — outside dockview, shows activeBottom panel
    <ActivityBar position="bottom" />  — outside dockview
    <FrameProgressBar />               — outside dockview, always visible
  </LandingPage>
</App>
```

`SidebarZone` and `BottomZone` are simple containers that mount the active panel's React component (Selections, Modifiers, etc.) for their bar's `active*` state.

### Default Layout (Startup)

```
┌──────────────────────────────────────────┐
│  AppBar                                  │
├────┬─────────────────────────────────┬───┤
│ A  │                                 │ A │
│ c  │                                 │ c │
│ t  │      3D Viewer (editor area)    │ t │
│ i  │           pinned tab            │ i │
│ v  │                                 │ v │
│ i  │                                 │ i │
│ t  │                                 │ t │
│ y  │                                 │ y │
├────┴─────────────────────────────────┴───┤
│  Bottom Activity Bar (empty)             │
├──────────────────────────────────────────┤
│  FrameProgressBar                        │
└──────────────────────────────────────────┘
```

Only the 3D Viewer is visible. Right activity bar exists but has no icons; bottom bar same. Identical feel to the current UI.

### With Plots Browser + Two Plots Open

```
┌──────────────────────────────────────────────────┐
│  AppBar                                          │
├────┬────────────┬───────────────┬────────────┬───┤
│ A  │ Plots      │   3D Viewer   │  Energy    │ A │
│ c  │ ─────────  │  (pinned tab) │  RDF (tab) │ c │
│ t  │ ● Energy   │               │            │ t │
│ i  │ ● RDF      │               │   chart    │ i │
│ v  │ ○ RMSD     │               │            │ v │
│ i  │ ○ Distance │               │            │ i │
│ t  │            │               │            │ t │
│ y  │ Open (2)   │               │            │ y │
│    │ ✕ Energy   │               │            │   │
│    │ ✕ RDF      │               │            │   │
├────┴────────────┴───────────────┴────────────┴───┤
│  Bottom Activity Bar                             │
├──────────────────────────────────────────────────┤
│  FrameProgressBar                                │
└──────────────────────────────────────────────────┘
```

Editor area is split: 3D Viewer (left group) + plot tabs (right group). Plots Browser sidebar shows which plots are open with filled-dot indicators.

## Activity Bar Behavior

Each activity bar (Left, Right, Bottom) is a strip of icons. Each icon represents a tool panel.

**Click an icon:**
- If the bar's `active*` is null → open this panel (set `active*` to the panel id)
- If the bar's `active*` already equals this panel → close it (set `active*` to null)
- If the bar's `active*` is a different panel → switch to this one (only one open at a time per bar)

**Drag an icon:**
- Drag-source: any activity bar icon
- Drop-target: any activity bar (Left, Right, Bottom)
- Drop-target invalid: the editor area (tools cannot become views)
- On drop: remove from source bar's `*BarIcons`, append to target bar's `*BarIcons`. If the panel was active in its source bar, set source `active*` to null.
- Visual: highlight valid drop zones, show "❌ Not a valid drop zone" overlay on the editor area.

**Reorder within a bar:** drag-and-drop within the same bar reorders icons.

**Active indicator:** the active icon in each bar gets a colored border (left border for Left bar, right border for Right bar, bottom border for Bottom bar).

## Editor Area (Dockview) Behavior

**Default content:** the 3D Viewer panel, registered as the only initial panel.

**3D Viewer panel:**
- Pinned (cannot be closed). Implemented by omitting the close button in its tab renderer or via dockview's panel API.
- Resizable by adjacent splits.
- Can be moved between groups; on re-parent, hook into dockview's `onDidLocationChange` to detect WebGL context loss and re-init the R3F renderer. Brief flash (<100ms) only on user-initiated moves.
- Can be made floating (Shift+drag) or popped out — both supported. R3F re-init handled the same way.

**Plot view panels:**
- One per opened plot. Panel id = `plot-{figureKey}`.
- Created when a user clicks an available plot in the Plots Browser, OR automatically when the server emits `figure_invalidate` with `operation="set"` for a key that is not currently open.
- New plots auto-open in the active editor group as a new tab.
- Closing a plot tab removes it from the editor area but does not delete the figure on the server. The plot remains in the Plots Browser as a closeable available entry.
- If the server emits `figure_invalidate` with `operation="delete"`, any open plot tab for that key is auto-closed.

**Editor-area drop rule:** activity bar icons cannot be dropped here. Plot tabs / 3D Viewer tab can be dropped anywhere within the editor area (split, stack, float).

## Tool Panel Specifications

### Panel: selections

- **Content:** existing `SelectionsPanel` component (unwrapped from current sidebar layout)
- **Default bar:** left

### Panel: modifiers

- **Content:** existing `SecondaryPanel` (modifiers variant)
- **Default bar:** left

### Panel: analysis

- **Content:** existing `SecondaryPanel` (analysis variant)
- **Default bar:** left

### Panel: geometries

- **Content:** existing `GeometryPanel`
- **Default bar:** left

### Panel: plots-browser (NEW)

- **Default bar:** left
- **Content:** new component `PlotsBrowserPanel.tsx`
- **Sections:**
  - **Available Plots** — every figure key from `GET /v1/rooms/{room_id}/figures`. Each row shows the key and a status dot (filled = currently open as a view, empty = not open). Click opens the plot as a view in the active editor group, or focuses its tab if already open.
  - **Currently Open** — the subset of plots that have an open editor tab. Each row has a × to close the editor tab.
- **Reactive updates:** subscribes to `figure_invalidate` to keep the list fresh.
- **Empty state:** "No plots available" placeholder.

### Panel: chat

- **Content:** existing `ChatWindow` component without the `react-rnd` wrapper
- **Default bar:** left
- **Floating mode:** users can Shift+drag the chat tab to float it over the 3D viewer. dockview handles drag/resize/clamping.

## Component Changes

### New Files

- `frontend/src/components/dockview/DockviewLayout.tsx` — Root dockview component, registers panel component factories, sets up the default layout (just 3D Viewer).
- `frontend/src/components/dockview/ActivityBar.tsx` — Three instances, parameterised by position (`left`/`right`/`bottom`). Renders icons from the activity bar slice, handles click + drag-and-drop.
- `frontend/src/components/dockview/SidebarZone.tsx` — Container that renders the active panel for a Left or Right bar.
- `frontend/src/components/dockview/BottomZone.tsx` — Container that renders the active panel for the Bottom bar.
- `frontend/src/components/PlotsBrowserPanel.tsx` — New tool panel listing available + currently open plots.
- `frontend/src/components/dockview/PlotView.tsx` — Thin wrapper that adapts existing plot rendering (Plotly + figure data fetch) to dockview's `IContentRenderer` interface.
- `frontend/src/components/dockview/ViewerView.tsx` — Thin wrapper for the 3D viewer panel, hooks into dockview lifecycle for WebGL context handling.
- `frontend/src/stores/slices/activityBarSlice.ts` — Zustand slice as described above.

### Modified Files

- `frontend/src/pages/landingPage.tsx` — Replace the current fixed-sidebar + Canvas + WindowManager layout with `<ActivityBar position="left">`, `<SidebarZone position="left">`, `<DockviewLayout>`, `<SidebarZone position="right">`, `<ActivityBar position="right">`, `<BottomZone>`, `<ActivityBar position="bottom">`. Drop the `drag-boundary-container` div.
- `frontend/src/pages/landingPage.tsx` AppBar (currently inline in this file, lines ~222–380) — Remove the AddPlot button (replaced by Plots Browser) and the chat icon (replaced by activity bar icon). Keep the rest.
- `frontend/src/components/ChatWindow.tsx` — Remove `react-rnd` wrapper. Export the chat content directly. Drop the dual fullscreen/draggable mode logic; floating is now handled by dockview.
- `frontend/src/components/FigureWindow.tsx` — Remove `react-rnd` wrapper and the figure-key dropdown selector. Each plot view renders one figure directly.

### Removed

- `frontend/src/components/WindowManager.tsx`
- `frontend/src/stores/windowManagerStore.ts`
- `frontend/src/components/AddPlotButton.tsx` (functionality moves to Plots Browser; auto-open replaces the explicit button)
- `frontend/src/components/SideBar.tsx` (replaced by the new ActivityBar + SidebarZone components)
- `frontend/src/components/PrimaryDrawer.tsx` (replaced by the new ActivityBar component)
- `formStore.ts`'s `selectedCategory` state — superseded by the activity bar slice
- `react-rnd` dependency

### Unchanged

- All Three.js components in `frontend/src/components/three/`
- Zustand store slices: connection, playback, scene, lock, UI (chatOpen / chatUnreadCount remain — they drive activity bar state for the chat icon)
- Socket.IO communication layer
- All backend code

## Three.js / WebGL Considerations

- The R3F `<Canvas>` lives inside a normal dockview panel. Resizing adjacent panels triggers R3F's built-in resize handling automatically.
- WebGL context is only at risk when the panel's DOM node is re-parented. Re-parent triggers: dragging the 3D Viewer to a different group, making it floating, popping it out.
- On context loss: hook into dockview's `onDidLocationChange` (and `onWebGLContextLost` on the canvas), unmount/remount the R3F tree. Brief flash (<100ms), only on user-initiated moves.
- Normal interaction (orbiting, zooming, frame scrubbing, playback) does not cause re-parenting and is unaffected.
- Browser WebGL context limit (8–16) is not a concern with a single 3D viewport.

## Migration Strategy

**Feature-flagged, incremental.** The big-bang approach was reconsidered as too risky for a refactor that touches the entire layout, removes a dependency, and adds a new library.

**Plan:**

1. Implement the new layout end-to-end behind a feature flag. Possible mechanisms (pick one during implementation):
   - URL query param: `?dockview=1`
   - User setting (server-side preference)
   - `localStorage` toggle
2. Both the old layout (`SideBar` + `WindowManager`) and the new (`ActivityBar` + `DockviewLayout`) coexist during this period.
3. Iterate based on feedback while the old layout remains the default.
4. Once the new layout is stable, flip the default. Keep the old layout reachable via the flag for one release as a safety valve.
5. Remove the old layout, `SideBar`, `PrimaryDrawer`, `WindowManager`, `windowManagerStore`, `react-rnd` dependency, and `formStore.selectedCategory` in a follow-up PR.

This sequence allows shipping the new UI to opt-in users early without destabilising the default experience.

## Dependencies

### Add

- `dockview-react` (latest, currently v5.2.0)

### Remove (in the cleanup PR after the flag flip)

- `react-rnd` (currently v10.5.2) — fully replaced by dockview's floating groups

## Open Questions (to settle during implementation)

- **Mobile fallback:** dockview is desktop-oriented. Decide whether mobile gets a stripped-down layout (e.g., the 3D viewer fullscreen with a single bottom drawer for tools) or simply uses the existing layout regardless of the flag.
- **Layout reset:** UX for "reset to default layout" — top-of-bar menu item, command palette entry, or dev-only?
- **AppBar reconciliation:** the AppBar still has buttons that interact with panels (e.g., the existing chat-unread badge). Confirm whether the badge moves to the activity bar chat icon or stays on the AppBar as a separate indicator.

## Future Considerations (Not In Scope)

- Layout persistence via `dockview.toJSON()` → server-side preference
- Multiple 3D viewports
- VS Code postMessage integration for theme sync, command palette, file-open from VS Code
- Custom dockview themes synced with MUI light/dark
- Keyboard shortcuts (Ctrl+1/2/3 to focus zones)
- Drag plots from the Plots Browser into the editor area as splits (currently click opens as a tab)
