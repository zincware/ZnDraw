# Dockview UI Follow-Up — Design Spec

**Date:** 2026-04-16
**Status:** Draft
**Parent spec:** [`2026-04-02-dockview-ui-redesign-design.md`](2026-04-02-dockview-ui-redesign-design.md)

## Overview

Follow-up to the main dockview redesign addressing the 7 issues in the post-merge smoke test (`ISSUES.md` at repo root) plus one layout bug found during brainstorming. Scope is frontend-only; no backend changes.

Summary of the direction:

- **Reclaim chrome space** with auto-hiding activity bars that stay reachable as drop targets via a 4px edge sliver that lights up during a panel-icon drag.
- **Theme dockview to MUI** and sync light/dark mode.
- **Fix the bottom-drawer shrink regression** so opening a bottom panel actually shrinks the 3D viewer.
- **Add popout + maximize** actions on dockview groups.
- **Open new plots beside the viewer**, not on top of it.
- **Simplify the plots browser** to one list.
- **Auto-resize plotly** when a panel is resized.

## Goals

- Default startup reclaims ≥30px vertical + ≥40px horizontal for the viewer (measured relative to the current main branch).
- Dockview panel chrome visually matches MUI (light and dark).
- Opening a bottom panel shrinks the viewer, same as opening a left/right panel.
- Every panel has popout + maximize affordances.
- Plot tabs open adjacent to the viewer by default, never on top of it.
- No loss of existing functionality (drag-icon-between-bars still works; chat is still reachable; the 3D viewer still closeable with leave-room cascade).

## Non-Goals (YAGNI)

- Layout persistence (already deferred in the parent spec).
- Customising dockview theme beyond the light/dark pair synced to MUI.
- Replacing dockview's native tab design or splitter gutters.
- Keyboard shortcuts for panel activation.
- A "close empty room" flow beyond what parent spec defined.

## Issue-by-issue design

### 1. Plotly auto-resize on panel size change

Dockview fires `IDockviewPanelProps.api.onDidDimensionsChange` with the new `{ width, height }` whenever the panel's container resizes (splitter drag, group move, popout, maximize). `PlotView` subscribes and calls `Plotly.Plots.resize(container)`. Debounce is unnecessary — dockview already throttles.

### 2. Panel headers follow MUI light/dark mode

Two parts:

1. **Theme sync.** Read MUI's color mode via `useColorScheme()` and toggle `className` on `DockviewReact` between `dockview-theme-light` and `dockview-theme-dark` so the base palette flips with MUI.
2. **MUI palette override.** New stylesheet `frontend/src/panels/dockview-mui.css` imported by `DockviewLayout.tsx`. Sets the subset of dockview CSS variables to MUI tokens:
   - `--dv-group-view-background-color` → `var(--mui-palette-background-default)`
   - `--dv-tabs-and-actions-container-background-color` → `var(--mui-palette-background-paper)`
   - `--dv-tab-background-color` → `var(--mui-palette-background-paper)`
   - `--dv-activegroup-visiblepanel-tab-background-color` → `var(--mui-palette-background-default)`
   - `--dv-tab-header-foreground-color` → `var(--mui-palette-text-primary)`
   - `--dv-inactivegroup-visiblepanel-tab-background-color` → `var(--mui-palette-action-hover)`
   - `--dv-separator-border` → `var(--mui-palette-divider)`
   - `--dv-paneview-active-outline-color` → `var(--mui-palette-primary-main)`

   The list above is the authoritative set for this change. Any additional `--dv-*` variables discovered during implementation that visibly diverge from MUI are added in the same stylesheet with a one-line comment stating the observed mismatch.

### 3. Popout button

Dockview renders one instance of `rightHeaderActionsComponent` per group automatically, so the popout button appears on every group (viewer group, plot groups, any group the user splits off). Clicking calls `api.addPopoutGroup(group)`. Popouts reuse the same MUI theme by serving the HTML shell from the existing static assets.

### 4. Maximize / fullscreen button

Second button in the same per-group `rightHeaderActionsComponent`. Toggles `api.maximizeGroup(group)` / `api.exitMaximizedGroup()` based on `api.hasMaximizedGroup(group)`. Icon flips between maximize and un-maximize. No additional state management required.

### 5. Chrome economy — auto-hiding activity bars with sliver drop zones

**Design:** `ActivityBar` has three visual states:

| State | When | Size | Visuals |
|---|---|---|---|
| Full | bar has ≥1 icon | 48×viewport (or 40×viewport for bottom) | icons rendered |
| Sliver | bar is empty AND no drag is in progress | 4px | transparent (invisible) |
| Hot | bar is empty AND a panel-icon drag is in progress | 4px | blue tint + primary-color border, pulsing |

State tracking: each `ActivityBar` instance listens on `window` for `dragstart`/`dragend` and sets local `isDragActive` when the drag payload includes `application/x-zndraw-panel-id` (the activity-bar-icon MIME — distinct from dockview's internal tab-drag MIME and from the plots-browser drag MIME, so it only trips on icon-between-bars drags). The existing `moveIconToBar` handler is unchanged — it's still the drop consumer.

**Chat default bar:** changes from `right` to `left`. Registry updates, no other callers affected.

**Welcome placeholder:** kept, restyled to use MUI `text.secondary` so it respects light/dark mode. Copy unchanged.

### 6. PlotsBrowser simplification

Drop the "Currently Open" section entirely. One list with the status dot (filled = currently open as a tab) conveys the same info. The close `×` button moves inline onto each open row's hover state:

```
● energy              — open (filled dot); hover reveals × to close
○ rdf                 — available (outlined dot); click to open
```

### 7. Plot opening position

`openPlotTab(api, figureKey, opts)` resolves the target group deterministically:

1. If an existing `plot-*` panel is open, use its group. When multiple groups contain plots, pick the group containing the panel whose id sorts first lexicographically among all `plot-*` panels in the layout (stable across sessions and independent of insertion order).
2. Else, find the viewer's group and place the new plot tab in a new group **to the right of it** (`position: { referenceGroup: viewerGroup, direction: "right" }`).
3. Else (no viewer — rare), fall back to `addPanel` with no position (dockview default placement).

Same logic used by `figureHandlers.onFiguresInvalidate` for the `op=set` auto-open and by the `onDidDrop` path (the latter respects the user's explicit drop target, so it overrides this rule).

### 8. Fix the bottom-drawer shrink bug (NEW)

**Symptom:** opening a bottom panel does not shrink the 3D viewer. Left/right panels do shrink it.

**Root cause:** `DockviewLayout.tsx:81` wraps `DockviewReact` in `<Box sx={{ flexGrow: 1, position: "relative", minWidth: 0, minHeight: 0 }}>`. `DockviewReact` renders inline inside this box, so its rendered height is determined by its content, not the flex-parent's current computed height. Horizontal shrinks still reach the component because the row-flex above it re-measures children every change. Vertical shrinks propagated from a sibling of the row-flex (BottomZone) don't retrigger a size query on a block-sized child.

**Fix:** wrap `DockviewReact` in an absolutely-positioned box so its size is decoupled from children:

```tsx
<Box sx={{ flexGrow: 1, position: "relative", minWidth: 0, minHeight: 0 }}>
  <Box sx={{ position: "absolute", inset: 0 }}>
    <DockviewReact ... />
  </Box>
</Box>
```

`DockviewReact`'s built-in `ResizeObserver` fires on the absolute wrapper; panels reflow; R3F canvas picks up the new size.

## Component changes

### New files

- `frontend/src/panels/dockview-mui.css` — CSS variable overrides mapping `--dv-*` to MUI tokens.
- `frontend/src/panels/groupActions.tsx` — popout + maximize action buttons (consumed by `DockviewReact`'s `rightHeaderActionsComponent` prop).

### Modified files

- `frontend/src/panels/registry.tsx` — chat `default.bar: "right"` → `"left"`.
- `frontend/src/panels/ActivityBar.tsx` — add three-state render (full/sliver/hot), listen on `window` drag events.
- `frontend/src/panels/SidebarZone.tsx` — already null-on-no-active; no change.
- `frontend/src/panels/BottomZone.tsx` — already null-on-no-active; no change.
- `frontend/src/panels/DockviewLayout.tsx` — absolute wrapper around `DockviewReact`; wire `rightHeaderActionsComponent`; import `dockview-mui.css`; read MUI color scheme and toggle theme class.
- `frontend/src/panels/PlotView.tsx` — subscribe to `onDidDimensionsChange` + `Plotly.Plots.resize`.
- `frontend/src/panels/PlotsBrowserPanel.tsx` — drop "Currently Open" section; inline `×` on hover.
- `frontend/src/panels/plotViewFactory.ts` — update `openPlotTab` to resolve target group per the rule above.
- `frontend/src/hooks/socketHandlers/figureHandlers.ts` — pass the updated positioning through (no behavior change here; it already calls `openPlotTab`).

### Unchanged

- Backend.
- `useLeaveRoom` and cascade logic.
- Socket event handling.
- All R3F scene components.

## Migration & validation

One PR. Commits split by issue where natural (shrink fix separately from theme sync, etc.). Validation:

- Typecheck + lint clean.
- Existing E2E suite still passes (the 13 currently failing pre-existing specs are out of scope — fix in a separate pass).
- New E2E assertions added to `e2e/dockview-layout.spec.ts`:
  - Bottom panel activation shrinks the viewer's rendered height (`evaluate` the `[data-testid=viewer-view]` bounding box before/after).
  - Dragging an icon to the right edge reveals a hot drop zone and placing it succeeds.
  - Popout button triggers `addPopoutGroup` (detectable via a new window).
  - Maximize button fills the dockview area and the un-maximize button restores it.
  - Opening a plot with no plot group present places it to the right of the viewer's group (assert `[data-testid=viewer-view]` bounding box `x + width ≤ plot-panel` bounding box `x`).
  - Theme sync: toggling MUI color mode flips `DockviewReact`'s class between `dockview-theme-light` and `dockview-theme-dark`.
  - Plotly resize: resizing a plot panel via splitter drag produces a matching `.plotly` SVG `viewBox` update (or a Plotly layout width change).
  - PlotsBrowser simplification: only one list is rendered; no element with a "Currently Open" heading.
- Manual smoke via playwright-cli with a screenshot to `/tmp/dockview-validation/followup-*.png`.

## Space reclaimed (design targets)

| Change | Vertical | Horizontal |
|---|---|---|
| Chat icon moves to left bar, right bar collapses to 4px sliver | — | −44px |
| Bottom activity bar collapses to 4px sliver when empty | −36px | — |
| **Net default startup** | **−36px** | **−44px** |

The viewer tab bar remains visible so popout + maximize buttons (and the viewer `×` close affordance) are always reachable. The MUI theme sync makes the tab bar read as part of the app, not as foreign chrome, which was the core complaint — the visual fix solves the feel without the geometry fix being necessary.

## Open questions

- Do we need a user-facing indicator that a hidden activity bar CAN be summoned (e.g., a tiny chevron at the edge), or is drag-reveal enough? Default: drag-reveal only, no persistent indicator. Revisit after one week of use.
- Should maximize persist across tab switches within the maximized group, or reset on tab change? Default: persist (matches dockview's native behavior).
