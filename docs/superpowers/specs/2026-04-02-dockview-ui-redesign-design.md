# Dockview UI Redesign — Design Spec

**Date:** 2026-04-02
**Status:** Draft (revised 2026-04-16 — full redesign; prior draft superseded)

## Overview

Replace ZnDraw's fixed sidebar + floating window layout with a VS Code-style two-tier architecture using `dockview-react`:

- **Tools** (Selections, Modifiers, Analysis, Geometries, Plots Browser, Rooms, Filesystem, Chat) live in single-panel-per-bar zones controlled by activity bar icons. Three activity bars exist (Left, Right, Bottom). Icons can be dragged between bars.
- **Views** (3D Viewer, plot tabs) live in a splittable editor area in the center. Tabs can be split horizontally/vertically, dragged between groups, floated over the grid, or popped out to a separate browser window.

The design follows VS Code's well-established conventions: tools configure, views show content, and the two tiers have different layout primitives. The 3D Viewer's lifecycle is bound to the current room — closing it leaves the room (cascade close of all plot tabs).

## Goals

- VS Code-like docking experience with three activity bars (Left/Right/Bottom) and a splittable editor area
- Activity bar icons can be dragged between bars (Left ↔ Right ↔ Bottom)
- Each bar shows at most one tool panel at a time (single-panel-per-bar behavior)
- Editor area supports splitting, tabbing, floating-over-grid, and pop-out windows for views
- **Plots Browser** panel lists all available plots; click opens as a tab in the active group, drag opens at any dockview drop target (split/stack/float)
- **Rooms** panel replaces the AppBar rooms menu; clicking a row switches rooms in-place (History API, no page reload)
- **Filesystem** panel replaces the full-page `/filesystem` route
- Closing the 3D Viewer leaves the room (cascade closes all plot tabs); switching rooms via the Rooms panel has identical cascade semantics
- New plots from the server auto-open as tabs in the active editor group
- Chat moves from floating overlay to the right activity bar (default closed)
- Single frontend source of truth for every panel: `PANELS` registry. Adding a new tool/view = one registry entry + one component file.
- **Full replacement, no artifacts.** No feature flag, no coexistence, no dead code. Final verification sweep confirms zero remnants.

## Non-Goals (YAGNI)

- **Pydantic-configurable panel layout.** Evaluated and dropped. Panel placement is UI chrome, not a visualization feature — the Pydantic-non-negotiable rule is about computation/API features (geometries, representations, analysis), not ephemeral UI preferences. If Python control is ever needed, it will be a Socket.IO command layer, not a config-endpoint problem — designed fresh when requested.
- **Layout persistence** (localStorage / server). Dockview supports `toJSON()`/`fromJSON()`; cheap to add later.
- **Per-user persisted layout.**
- **Multiple 3D viewports.** Architecturally supported by dockview but not implemented (no camera management for multiple viewports).
- **VS Code native multi-panel integration.** Existing `vscode-zndraw` single-webview approach continues to work as-is. A separate spec may later cover postMessage-based VS Code enhancements.
- **Bar icon ordering persistence.**
- **Custom dockview themes beyond a light/dark mapping to MUI.**
- **Keyboard shortcuts for panel focus.**
- **Mobile fallback.** ZnDraw is desktop-first (playback bar is already desktop-only). Deferred.

## Backend Changes

**None.** The existing REST + Socket.IO surface covers everything this redesign needs:

| Panel | Existing backend surface |
|-------|--------------------------|
| Plots Browser | `GET /v1/rooms/{room_id}/figures` + `figure_invalidate` socket event |
| Selections | existing selection-groups + frame-selection routes + `SelectionGroupsInvalidate` |
| Modifiers / Analysis | existing extension job system |
| Geometries | existing CRUD + `GeometryInvalidate` |
| Chat | existing message routes + `MessageNew` / `MessageEdited` |
| Rooms | existing room management routes + invalidation events |
| Filesystem | existing filesystem-browser routes |

This is a frontend-only refactor.

## Library

**dockview-react** (MIT, zero dependencies, currently v5.2.x)

Architectural decision, verified via context7:

- A single `DockviewReact` instance manages the **editor area only**.
- Activity bars and sidebar/bottom zones are **not** dockview — a DockviewReact grid cannot enforce "single panel per bar." Custom containers + a small Zustand slice manage them.
- Nested DockviewReact instances were considered and rejected as overkill for sidebars that hold at most one panel.

Used dockview capabilities:

- Splittable grid + tab groups for the editor area
- `floatingGroupBounds="boundedWithinViewport"` for clamping
- `addPanel({ floating: true })` and `addFloatingGroup(group, { position })` for in-app floating panels
- Shift+drag-tab gesture for users to float a tab
- `addPopoutGroup(group, ...)` for separate browser windows
- External drag-drop integration (plot browser rows → editor area) via `onUnhandledDragOverEvent` + `onDidDrop`
- React 19 compatible, TypeScript core

## Architecture

### Two Tiers

**Tier 1 — Tools (custom, outside dockview):**

Three activity bars (Left, Right, Bottom) and their associated panel zones are managed by a Zustand slice + plain React containers. Reasons:

1. Each bar holds at most one panel — a DockviewReact grid would allow splitting and closing arbitrarily, defeating the constraint.
2. Bar icons must be dragged between bars (HTML5 drag-and-drop); dockview's drag system is built around panels, not icons.
3. The interaction model differs fundamentally from dockview's grid.

Zustand slice:

```ts
type ActivityBarSlice = {
  leftBarIcons: PanelId[];       // ordered, panel ids from PANELS registry
  rightBarIcons: PanelId[];
  bottomBarIcons: PanelId[];
  activeLeft: PanelId | null;    // which panel is visible in the left zone
  activeRight: PanelId | null;
  activeBottom: PanelId | null;
  // actions
  moveIconToBar(id: PanelId, bar: 'left' | 'right' | 'bottom', index?: number): void;
  toggleActive(bar: 'left' | 'right' | 'bottom', id: PanelId): void;
  resetLayout(): void;  // restore defaults from PANELS registry
};
```

Initial state is derived from the `default` entry on each panel in the `PANELS` registry (see below).

**Tier 2 — Views (managed by dockview):**

A single `DockviewReact` instance manages the editor area between the AppBar and the FrameProgressBar. The 3D Viewer and any opened plot views live inside it.

### Panel Registry (frontend source of truth)

Single file, `frontend/src/panels/registry.ts`:

```ts
import type { ComponentType } from 'react';
import type { IDockviewPanelProps } from 'dockview-react';
import type { SvgIconComponent } from '@mui/icons-material';

export type PanelId =
  | 'selections' | 'modifiers' | 'analysis' | 'geometries'
  | 'plots-browser' | 'rooms' | 'filesystem' | 'chat'
  | 'viewer';

type PanelDefault = {
  bar: 'left' | 'right' | 'bottom' | 'editor';
  active?: boolean;
  order?: number;
};

type ToolPanelDef = {
  kind: 'tool';
  icon: SvgIconComponent;
  label: string;
  component: ComponentType;
  default: PanelDefault;
};

type ViewPanelDef = {
  kind: 'view';
  title: string;
  component: ComponentType<IDockviewPanelProps>;
  closable: boolean;
  cascadeOnClose?: (string | RegExp)[];   // panel ids to close when this one is closed
  onClose?: () => void;                   // side-effects: e.g., viewer-close triggers leave-room
};

export type PanelDef = ToolPanelDef | ViewPanelDef;

export const PANELS: Record<PanelId, PanelDef> = {
  selections:      { kind: 'tool', icon: FilterCenterFocus, label: 'Selections',    component: SelectionsPanel,    default: { bar: 'left',  order: 0 } },
  modifiers:       { kind: 'tool', icon: Build,             label: 'Modifiers',     component: ModifiersPanel,     default: { bar: 'left',  order: 1 } },
  analysis:        { kind: 'tool', icon: Analytics,         label: 'Analysis',      component: AnalysisPanel,      default: { bar: 'left',  order: 2 } },
  geometries:      { kind: 'tool', icon: Category,          label: 'Geometries',    component: GeometryPanel,      default: { bar: 'left',  order: 3 } },
  'plots-browser': { kind: 'tool', icon: ShowChart,         label: 'Plots',         component: PlotsBrowserPanel,  default: { bar: 'left',  order: 4 } },
  rooms:           { kind: 'tool', icon: MeetingRoom,       label: 'Rooms',         component: RoomsPanel,         default: { bar: 'left',  order: 5 } },
  filesystem:      { kind: 'tool', icon: Folder,            label: 'Files',         component: FilesystemPanel,    default: { bar: 'left',  order: 6 } },
  chat:            { kind: 'tool', icon: Chat,              label: 'Chat',          component: ChatPanel,          default: { bar: 'right', order: 0 } },

  viewer:          {
    kind: 'view', title: '3D Viewer', component: ViewerView, closable: true,
    cascadeOnClose: [/^plot-.+$/],
    onClose: () => leaveCurrentRoom(),   // drops room, clears URL, shows welcome screen
  },
};
```

Adding a new static panel = add to `PanelId` union + one registry entry + one component file. TypeScript catches missing bindings at compile time.

Plot views are dynamic (`plot-{figureKey}`) and are not in the static registry — they are created via a factory when a user opens a plot or the server emits `figure_invalidate op=set`.

### Default Layout Table

| id | kind | default bar | order | active on start | closable | Notes |
|---|---|---|---|---|---|---|
| `selections` | tool | left | 0 | no | yes | existing `SelectionsPanel` |
| `modifiers` | tool | left | 1 | no | yes | existing `SecondaryPanel` (modifiers variant) |
| `analysis` | tool | left | 2 | no | yes | existing `SecondaryPanel` (analysis variant) |
| `geometries` | tool | left | 3 | no | yes | existing `GeometryPanel` |
| `plots-browser` | tool | left | 4 | no | yes | **new** — lists available + open plots |
| `rooms` | tool | left | 5 | no | yes | **new** — in-place room switching via History API |
| `filesystem` | tool | left | 6 | no | yes | **new** — replaces full-page `/filesystem` route |
| `chat` | tool | right | 0 | no | yes | moved from AppBar floating overlay |
| `viewer` | view | editor | — | yes (pinned initial tab) | **yes (close = leave room)** | R3F canvas, always present while in a room |
| `plot-{key}` | view (dynamic) | editor | — | auto-open on `figure_invalidate op=set` | yes | one per open figure, auto-close on `op=delete` |

**Bottom bar:** empty by default. Users drag any tool icon into it.

### Root Layout

```
<App>
  <LandingPage>
    <AppBar />                         — outside dockview, always visible
    <Box class="workspace">            — flex row
      <ActivityBar position="left" />  — outside dockview
      <SidebarZone position="left" />  — outside dockview, mounts activeLeft panel
      <DockviewLayout />               — editor area (viewer + plot views)
      <SidebarZone position="right" /> — outside dockview, mounts activeRight panel
      <ActivityBar position="right" /> — outside dockview
    </Box>
    <BottomZone />                     — outside dockview, mounts activeBottom panel
    <ActivityBar position="bottom" />  — outside dockview
    <FrameProgressBar />               — outside dockview, always visible, pinned
  </LandingPage>
</App>
```

`SidebarZone` and `BottomZone` render `PANELS[activeX].component` for their bar's current `active*` state, or nothing when `active*` is null.

### Default Startup Layout

Only visible chrome is the AppBar, left activity bar (7 icons), right activity bar (1 icon: chat), bottom activity bar (empty), viewer filling the editor area, and FrameProgressBar. Identical feel to the current UI; the left sidebar opens on first icon click.

```
┌──────────────────────────────────────────┐
│  AppBar                                  │
├────┬─────────────────────────────────┬───┤
│ A  │                                 │ A │
│ c  │                                 │ c │
│ t  │      3D Viewer (editor area)    │ t │
│ i  │                                 │ i │
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

### With Plots Browser Open + Two Plots Open

```
┌──────────────────────────────────────────────────┐
│  AppBar                                          │
├────┬────────────┬───────────────┬────────────┬───┤
│ A  │ Plots      │   3D Viewer   │  Energy    │ A │
│ c  │ ─────────  │  (tab)        │  RDF (tab) │ c │
│ t  │ ● Energy   │               │            │ t │
│ i  │ ● RDF      │               │            │ i │
│ v  │ ○ RMSD     │               │   chart    │ v │
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

## Activity Bar Behavior

Each activity bar (Left, Right, Bottom) is a strip of icons. Each icon represents a tool panel.

**Click an icon:**
- If the bar's `active*` is null → open this panel
- If the bar's `active*` already equals this panel → close it (set `active*` to null)
- If a different panel is active → switch (only one open at a time per bar)

**Drag an icon:**
- Drag-source: any activity bar icon
- Drop-target valid: any activity bar
- Drop-target invalid: the editor area (tools cannot become views)
- On drop: remove from source bar's icon list, append to target bar's icon list. If the panel was active in its source bar, set source `active*` to null.
- Visual feedback: highlight valid drop zones; show "❌ Not a valid drop zone" overlay on the editor area.

**Reorder within a bar:** drag-and-drop within the same bar reorders icons.

**Active indicator:** a colored border (left border for Left bar, right for Right, top for Bottom) on the active icon.

## Editor Area Behavior

**Initial content:** the 3D Viewer panel, the only initial panel. No plots open by default.

**3D Viewer panel:**
- `closable: true`. The tab has an X.
- **Close = leave current room (F1 semantic):**
  1. If any `plot-*` tabs are open, show a confirmation toast: `"Leave room? N plot(s) will close. [Leave] [Cancel]"`. Empty editor = no confirm.
  2. On confirm: dockview closes all `plot-*` panels (via `cascadeOnClose: [/^plot-.+$/]`), then closes the viewer panel.
  3. `onClose` side-effect triggers `leaveCurrentRoom()`: disconnects the room socket, clears room state in Zustand, updates URL to `/` via `history.pushState`.
  4. Editor area shows a welcome placeholder: `"No room selected. Open the Rooms panel to pick one."` — rendered by `DockviewLayout` when no panels are present.
- Resizable by adjacent splits.
- Can be moved between groups, floated (Shift+drag), or popped out. On re-parent, hook into dockview's `onDidLocationChange` to detect WebGL context loss and re-init the R3F renderer. Brief flash (<100ms) only on user-initiated moves.

**Plot view panels:**
- Panel id: `plot-{figureKey}`.
- Created when:
  - A user clicks a row in the Plots Browser → panel opens as a tab in the **active editor group**
  - A user drags a row from the Plots Browser onto the editor area → panel opens at the drop target (split/stack/float — any dockview drop zone)
  - The server emits `figure_invalidate op=set` for a key that is not currently open → auto-opens as a tab in the active editor group
- Closing a plot tab removes it from the editor area but does **not** delete the figure on the server. The plot remains available in the Plots Browser.
- If the server emits `figure_invalidate op=delete`, any open plot tab for that key is auto-closed.
- Plot tabs are reparentable, splittable, floatable, and popout-able.

**Room switching semantics (same cascade as F1):**
- When the user clicks a room row in the Rooms panel (different from current): close all `plot-*` tabs, `history.pushState('/rooms/{new_id}')`, viewer re-inits for the new room, plots auto-open for any existing figures in the new room.
- No confirmation toast for switching — the action is as explicit as closing.

**Editor-area drop rule:** activity bar icons cannot be dropped here; plots-browser rows can; plot tabs / viewer tab can be dropped anywhere within the editor area.

## Plots Browser Panel

- **Default bar:** left
- **Sections:**
  - **Available Plots** — every figure key from `GET /v1/rooms/{room_id}/figures`. Each row shows the key and a status dot (filled = open, empty = available). Click → open in active editor group, or focus if already open. Drag → any dockview drop target.
  - **Currently Open** — subset with open editor tabs. Each row has a × to close.
- **Reactive updates:** subscribes to `figure_invalidate` to refresh.
- **Empty state:** "No plots available" placeholder.
- **Drag integration:** rows emit a dockview-compatible drag payload (`application/x-dockview-panel` or similar) consumed via `onUnhandledDragOverEvent` + `onDidDrop` on `DockviewReact`.

## Rooms Panel

- **Default bar:** left
- **Content:** list of accessible rooms, current one highlighted, per-row actions (star, copy link, leave/delete).
- **Row click:** F1 cascade (close viewer + plots) → load new room → `history.pushState('/rooms/{id}')`.
- **Subscribes to** room-list invalidation events.
- **Empty state:** "No rooms available — create one via the Filesystem panel."

## Filesystem Panel

- **Default bar:** left
- **Content:** reuses existing `FilesystemBrowser` component logic, stripped of page-level chrome (no standalone title/header).
- **File action:** opens `LoadFileDialog` → creates a new room → F1 cascade (close current viewer/plots) → switch to the new room.
- **Removes:** `frontend/src/pages/filesystemBrowser.tsx` and its route. The `/filesystem` router entry redirects to current room + opens the filesystem panel.

## Chat Panel

- **Default bar:** right
- **Content:** existing `ChatWindow` content (message list + input), unwrapped from the `react-rnd` overlay and renamed to `ChatPanel`.
- **Unread indicator:** the chat icon in the right activity bar shows the existing unread badge.
- **No floating mode in v1.** The sidebar container is not dockview, so Shift+drag-to-float does not apply. Users who want chat on top of the viewer can move the chat icon to the bottom bar (still docked) or wait for a future "detach to floating group" action.

## FrameProgressBar

**Unchanged.** Fixed at the bottom of the viewport, always visible. Not part of dockview, not draggable, not closable.

## AppBar Changes

**Removed:** `AddPlotButton` (functionality → Plots Browser), chat icon + badge (→ right-bar chat icon + badge), `RoomManagementMenu` (→ Rooms panel).

**Kept:** title, mode buttons (drawing/editing), theme toggle, connection info, screenshot, file upload, user profile, "..." menu.

**"..." menu additions:** `Reset layout` — calls `activityBarSlice.resetLayout()` to restore defaults from `PANELS`, clears dockview to initial (viewer only), closes floating/popout panels.

## Three.js / WebGL Considerations

- R3F `<Canvas>` lives inside a normal dockview panel. Resizing adjacent panels triggers R3F's built-in resize handling automatically.
- WebGL context is only at risk when the panel's DOM node is re-parented: dragging the viewer to another group, floating, or popping out.
- On context loss: hook `onDidLocationChange` + `onWebGLContextLost`, unmount/remount the R3F tree. Brief flash (<100ms), user-initiated moves only.
- Normal interaction (orbiting, zooming, frame scrubbing, playback) does not cause re-parenting and is unaffected.
- Browser WebGL context limit (8–16) is not a concern with a single viewport.

## Component Changes

### New Files

- `frontend/src/panels/registry.ts` — `PANELS` registry, `PanelId` union, types.
- `frontend/src/panels/ActivityBar.tsx` — parameterized by position; renders icons from the Zustand slice, handles click + drag-and-drop.
- `frontend/src/panels/SidebarZone.tsx` — container for left/right zones.
- `frontend/src/panels/BottomZone.tsx` — container for bottom zone.
- `frontend/src/panels/DockviewLayout.tsx` — root dockview component. Registers `viewer` + `plotView` component factories, sets up initial layout (viewer only), handles welcome state when empty.
- `frontend/src/panels/ViewerView.tsx` — thin wrapper for the 3D viewer; hooks dockview lifecycle for WebGL context handling; owns the `onClose` → `leaveCurrentRoom` side-effect.
- `frontend/src/panels/PlotView.tsx` — plot view adapter; fetches figure data, renders Plotly.
- `frontend/src/panels/plotViewFactory.ts` — `createPlotView(figureKey)` for dynamic instantiation.
- `frontend/src/panels/PlotsBrowserPanel.tsx` — lists available + open plots; drag source.
- `frontend/src/panels/RoomsPanel.tsx` — lists rooms; row click → F1 cascade + `history.pushState`.
- `frontend/src/panels/FilesystemPanel.tsx` — unwrapped filesystem browser; file action → new-room flow.
- `frontend/src/panels/ChatPanel.tsx` — unwrapped chat content.
- `frontend/src/stores/slices/activityBarSlice.ts` — Zustand slice.
- `frontend/src/hooks/useLeaveRoom.ts` (or similar) — room-leave cascade logic.

### Modified Files

- `frontend/src/pages/landingPage.tsx` — replace the fixed-sidebar + Canvas + WindowManager layout with the root layout described above. Drop `drag-boundary-container`. Remove AppBar buttons (AddPlot, chat, rooms menu). Add "Reset layout" to "..." menu.

### Removed Files / Symbols

- `frontend/src/components/SideBar.tsx`
- `frontend/src/components/PrimaryDrawer.tsx`
- `frontend/src/components/WindowManager.tsx`
- `frontend/src/components/FigureWindow.tsx`
- `frontend/src/components/AddPlotButton.tsx`
- `frontend/src/stores/windowManagerStore.ts`
- `frontend/src/formStore.ts` (`selectedCategory`, `selectedExtensions` — superseded by activity bar slice)
- `frontend/src/pages/filesystemBrowser.tsx`
- `chatOpen` / `setChatOpen` / `chatUnread` setters in `uiSlice` that become unused after wiring chat as an activity-bar panel (verify: the unread **count** state stays; the open/close toggle is gone)
- `react-rnd` dependency (`package.json`)
- `ChatWindow.tsx` — replaced by `ChatPanel.tsx` (the content logic moves; the react-rnd wrapper file is deleted)

### Unchanged

- All Three.js components in `frontend/src/components/three/`
- Zustand slices: connection, playback, scene, lock, UI (minus chat-open toggle)
- Socket.IO communication layer
- All backend code

## Migration Plan (big bang, split commits)

**One PR, one merge. Split into reviewable commits. Final commit is a verification sweep.**

1. Add `dockview-react`. Create `frontend/src/panels/` with registry types, empty component stubs.
2. Implement `ActivityBar`, `SidebarZone`, `BottomZone`, `activityBarSlice`.
3. Implement `DockviewLayout`, `ViewerView`, `PlotView`, `plotViewFactory`.
4. Implement new panels: `PlotsBrowserPanel`, `RoomsPanel`, `FilesystemPanel`.
5. Unwrap `ChatWindow` → `ChatPanel`; drop `react-rnd` imports from chat.
6. Wire F1 close-cascade + URL history (leave-room hook, cascade on `cascadeOnClose`, welcome placeholder).
7. Replace `landingPage.tsx` layout. Update AppBar (remove chat/figures/rooms buttons, add "Reset layout").
8. **Delete removed files/symbols** (see Removed Files / Symbols). Remove `react-rnd` from `package.json`; regenerate lockfile.
9. **Verification sweep** — must pass before merge:
   - `bun run typecheck` → 0 errors
   - `bun run lint` → 0 errors
   - Search must return 0 matches:
     ```
     rg -i 'react-rnd|WindowManager|PrimaryDrawer|selectedCategory|openWindow|FigureWindow|AddPlotButton|SideBar\.tsx|windowManagerStore|formStore|filesystemBrowser\.tsx'
     ```
   - `bun run test` (frontend) + `uv run pytest` → all pass or tests updated
   - Manual smoke: default startup, open each left-bar panel, drag icon to right/bottom, open a plot (click + drag variants), split a plot tab, float a tab, popout a tab, close viewer (confirm leave-room cascade), switch rooms, reset layout

## Dependencies

### Add

- `dockview-react` (~v5.2.x, MIT, zero runtime dependencies)

### Remove

- `react-rnd` (~v10.5.2) — fully replaced by dockview's floating groups for chat/plot floating

## Future Considerations (Not In Scope)

- **Layout persistence** via `dockview.toJSON()` → user settings / server-side preference.
- **Per-user persisted layout.**
- **Python runtime panel control** — Socket.IO command layer (`UIOpenPanel`, `UICloseTab`, `UISwitchActive`). Designed when first concretely requested.
- **Multiple 3D viewports.**
- **VS Code postMessage integration** — theme sync, command palette, file-open from VS Code.
- **Custom dockview themes synced with MUI light/dark.**
- **Keyboard shortcuts** (Ctrl+1/2/3 to focus zones, etc.).
- **Mobile fallback.**

## Open Questions (to settle during implementation)

- **`chatOpen` audit:** confirm whether any UI other than the old floating overlay reads `chatOpen`/`setChatOpen` (e.g., AppBar badge). If only the overlay, remove; if read elsewhere, migrate the readers to derive from activity bar active state.
- **Welcome placeholder copy / design:** plain text first pass; polish later.
- **Reset Layout confirmation:** show a "Reset layout?" confirm toast before clearing, or silent reset? Default: silent (consistent with dockview's drag interactions being instant).
