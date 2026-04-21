# Dockview UI Fixes — Design Spec

**Date:** 2026-04-17
**Status:** Draft
**Parent spec:** [`2026-04-02-dockview-ui-redesign-design.md`](2026-04-02-dockview-ui-redesign-design.md) (via [`2026-04-16-dockview-followup-design.md`](2026-04-16-dockview-followup-design.md))
**Companion specs:**
- [`2026-04-17-internal-providers-design.md`](2026-04-17-internal-providers-design.md) — covers issue #5 (default filesystem via `@internal` providers)
- [`2026-04-17-dockview-perf-audit-design.md`](2026-04-17-dockview-perf-audit-design.md) — covers issue #2 (performance regression)

## Overview

Second follow-up to the dockview redesign. Addresses three issues from `ISSUES.md` lines 9, 11-16, 17. Frontend-only.

- **Resizable sidebars & bottom zone** (issue #1) — replace the hard-coded `width: 320` / `height: 260` with per-bar sizes stored in `activityBarSlice`; 4 px drag handles on the inner edges.
- **Rooms panel feature parity** (issue #3) — VS-Code-style compact list with search, header actions (new empty, upload), per-row `⋮` menu (set template, duplicate, lock, copy link, download, delete). Reactive via existing socket handlers; no refresh button.
- **Wider snap zones during drag** (issue #4) — empty activity bars widen from 4 px to 56 px the moment a panel-icon drag begins, with distinct "hot" and "over-zone" visual states.

Performance (issue #2) and default filesystem (issue #5) are scoped to the companion specs so this PR stays small and reviewable.

## Goals

- Sidebars and bottom zone can be resized by dragging their inner edge; widths survive close/open of the same panel within a session.
- Rooms panel exposes all room-management actions available on the `/rooms/` page (minus the DataGrid chrome) with the same socket-driven live updates.
- A dragged panel icon reliably lands on a target bar without pixel-hunting a 4 px strip.

## Non-Goals (YAGNI)

- **Layout persistence** across page reloads (localStorage or server). Widths reset on reload — consistent with the parent spec's non-goal.
- **Retiring the `/rooms/` page.** Deferred; decided explicitly.
- **Adding byte-size to the `Room` schema.** Verified with product — only frame counts are available.
- **Performance regressions.** Covered in the perf-audit spec.
- **Filesystem registration.** Covered in the internal-providers spec.
- **Pydantic-configurable panel sizes.** Layout is UI chrome, not a computation feature — parent spec's carve-out still applies.
- **Keyboard shortcuts** for resize or snap.
- **Mobile fallback.**

## Issue-by-issue design

### Issue #1 — Resizable sidebar & bottom zones

**State addition to `activityBarSlice.ts`:**

```ts
type ActivityBarSlice = /* existing */ & {
  leftWidth: number;    // px, default 320
  rightWidth: number;   // px, default 320
  bottomHeight: number; // px, default 260
  setBarSize(bar: BarPosition, px: number): void;  // clamped to [MIN, MAX]
};
```

Bounds are single constants, not Pydantic-configurable:

```ts
export const SIDEBAR_MIN_PX = 200;
export const SIDEBAR_MAX_PX = 640;
export const BOTTOM_MIN_PX = 120;
export const BOTTOM_MAX_PX = 560;
```

**Persistence:** zustand, session-only. `setBarSize` touches the slice only; nothing is written to localStorage or the server. A page reload resets to the defaults above.

**Rationale for session-only:** closing and re-opening the same panel feels stable (matches user intuition), but we don't invest in cross-session persistence until it's asked for — consistent with the parent spec's "Layout persistence" non-goal.

**Drag handle:**

- Each zone renders a 1 px visible handle on its inner edge (right edge for `SidebarZone position="left"`, left edge for `position="right"`, top edge for `BottomZone`). The handle's hit area is expanded to 8 px via a transparent wrapper so the target is easy to grab without making the visible chrome noisy.
- Cursor: `col-resize` for sidebars, `row-resize` for bottom.
- Implementation: plain `pointerdown` / `pointermove` / `pointerup` with `setPointerCapture` on the handle. No new npm dependency.
- `pointermove` calls `setBarSize(bar, newPx)`. Zustand shallow equality means only the zone re-renders; dockview's internal `ResizeObserver` picks up the flex-parent change and re-lays out the viewer.

**Interaction with the bottom-drawer shrink fix** (follow-up #1, line 108 of that spec): that fix wraps `DockviewReact` in an absolutely-positioned box so vertical sibling resizes propagate. This design reuses that wrapper unchanged; the new splitter is outside it.

### Issue #3 — Rooms panel feature parity

**Shape:** option A from mockup `/tmp/dockview-validation/07-rooms-mockups.png`. VS-Code-style compact list.

**Row layout:**

- Line 1 (weight 500, primary text): `room.description ?? room.id`
- Line 2 (12 px, `text.secondary`): `<id-prefix-8> · <frame_count> frames` — byte-size is not available on the `Room` schema and is out of scope to add.
- Right edge: lock icon (green open / red closed), hover-revealed `⋮` menu button.
- Current room: primary-tint background + 2 px left border in primary color.

**Header buttons** (replaces today's single list):

- `+ New empty room` — `createRoom({ room_id: crypto.randomUUID() })`, then `navigate(/rooms/<id>)`. No file, no `copy_from`.
- `⇪ Upload file` — hidden `<input type="file">` ref; on change, creates a new room + uploads the file (reuses the existing `handleFiles` flow from `pages/roomList.tsx`).
- **No refresh button.** `roomsStore` is already live-updated via the existing `onRoomUpdate` / `onRoomDelete` socket handlers in `frontend/src/hooks/socketHandlers/roomHandlers.ts`. Adding a manual refresh would imply the store can get stale, which it can't.

**Search:** `<TextField>` under the header. Client-side filter against `roomsArray` on `room.id` and `room.description` (case-insensitive substring, same semantics as the `/rooms/` page search box).

**Per-row `⋮` menu** (MUI `Menu`):

- **Set as template / Remove template** — `setDefaultRoom(isDefault ? null : roomId)`.
- **Duplicate room** — opens the existing `DuplicateRoomDialog` component.
- **Lock / Unlock** — `updateRoom(id, { locked: !locked })`.
- **Copy link** — `navigator.clipboard.writeText(window.location.origin + '/rooms/' + id)`.
- **Download current frame / all frames** — `downloadFrames({ roomId, indices? })`.
- **Delete room** — new. Requires a `DELETE /v1/rooms/{room_id}` on the backend; check at implementation time. If the endpoint exists: confirmation dialog → call → `roomsStore.removeRoom` fires via the socket handler. If not: disable the menu item with a tooltip "not yet supported".

**Drag-and-drop:** dropping a trajectory file onto the Rooms panel body creates a new room and uploads, same as the `/rooms/` page.

**Create-copy-from flow** (issue #3 bullet "create copy from room") — reuses the existing `DuplicateRoomDialog` triggered from the per-row menu. No new UI surface.

### Issue #4 — Wider snap zones during drag

Matches mockup `/tmp/dockview-validation/08-snap-mockups.png`.

**State machine in `ActivityBar.tsx`:**

| `isDragActive` | icons | cursor over bar | rendered size | visuals |
|---|---|---|---|---|
| false | 0 | — | 4 px | transparent sliver |
| false | ≥1 | — | 48 / 40 px | normal bar |
| true | 0 | false | **56 px** | `--mui-palette-primary-main` at 15 % fill + 2 px dashed border + pulse |
| true | 0 | **true** | **56 px** | 35 % fill + 2 px solid border + hint label |
| true | ≥1 | false | 48 / 40 px | dashed outline; no resize |
| true | ≥1 | **true** | 48 / 40 px | solid outline + hint label |

**Drag detection:** existing `window` `dragstart`/`dragend` listeners in `ActivityBar` already set `isDragActive` on `application/x-zndraw-panel-id`. No change.

**Cursor-over-self:** new local state, set in the bar element's `onDragOver`, cleared in `onDragLeave` / `onDragEnd`.

**Animation:** existing `pulse` keyframes (from follow-up #1) — reused on the wider zone.

**Hint label:** absolutely-positioned 11 px `primary.main` text at the zone's midpoint, visible only in the over-zone state.

- Left bar: `Drop to dock left`
- Right bar: `Drop to dock right`
- Bottom bar: `Drop to dock bottom`

**Scope:** hot zones appear on **empty bars only**. Already-populated bars get an outline during drag but don't resize, so the existing icons stay in place and the user can still drop onto a specific icon.

**Fallback:** if the browser doesn't fire drag events reliably (rare), the click-to-toggle + drag-on-icon paths are untouched — the UI reverts to today's 4 px target.

## Component changes

### Modified

- `frontend/src/stores/slices/activityBarSlice.ts` — add `leftWidth`, `rightWidth`, `bottomHeight`, `setBarSize`.
- `frontend/src/panels/SidebarZone.tsx` — read width from slice, render drag handle, `pointermove` handler.
- `frontend/src/panels/BottomZone.tsx` — same for height.
- `frontend/src/panels/ActivityBar.tsx` — expand state machine for 56 px hot / 35 %-fill over-zone; cursor-over-self local state; hint labels.
- `frontend/src/panels/RoomsPanel.tsx` — full rewrite to compact list with search + header buttons + per-row menu. Borrows logic from `pages/roomList.tsx` (create, upload, lock toggle, template toggle, duplicate dialog).
- `frontend/src/panels/registry.tsx` — export size constants.

### New

- `frontend/src/panels/roomsHeaderActions.tsx` — `+ New empty` and `⇪ Upload file` buttons.
- `frontend/src/panels/roomRowMenu.tsx` — `⋮` MUI `Menu` with the six actions.

### Unchanged

- `DockviewLayout.tsx` (absolute-wrapper fix from follow-up #1 stays).
- Socket handlers (`roomHandlers.ts` already provides the live updates the panel needs).
- Backend.

## Migration & validation

One PR. Commits split by issue where natural (resize separate from rooms-panel rewrite, etc.).

**Typecheck + lint:** zero errors.

**Existing E2E suite:** still passes.

**New E2E assertions in `frontend/tests/e2e/dockview-layout.spec.ts`:**

- **Resize:**
  - Drag the inner edge of the left sidebar by 100 px; assert `SidebarZone` bounding box widens by 100 px (± 1).
  - Resize below `SIDEBAR_MIN_PX` → clamped.
  - Resize above `SIDEBAR_MAX_PX` → clamped.
  - Resize, close the panel, re-open: width preserved.
  - Resize, reload page: width reset to default.
- **Rooms panel:**
  - Search filter narrows the list.
  - `+ New empty` creates a room; socket update populates the list.
  - Per-row `⋮` opens the menu; `Set as template` toggles the star; `Lock/Unlock` flips the icon.
  - Drag a `.xyz` file onto the panel body → new room created with the file uploaded.
- **Snap zones:**
  - Start a panel-icon drag; assert empty right-bar width grows to 56 px.
  - Move cursor into the hot zone; assert it carries `data-over-zone="true"` (or the equivalent attribute) and the hint label is visible.
  - Drop in the hot zone → panel moved to the target bar.
  - Cancel drag (Escape / drop outside) → bar shrinks back to 4 px.

**Manual smoke** via playwright-cli with screenshots written to `/tmp/dockview-validation/ui-fixes-*.png`.

## Open questions

- **Delete-room endpoint existence:** confirm `DELETE /v1/rooms/{room_id}` during implementation. If it doesn't exist, the `⋮` menu item is disabled + tooltip; no new backend work in scope.
