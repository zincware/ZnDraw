# Dockview UI Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement three independent UI fixes on the `spec/dockview-ui-redesign` branch: resizable sidebar/bottom zones, VS Code-style rooms panel with feature parity, and wider snap-drop zones during panel-icon drags.

**Architecture:** Frontend-only. Custom `pointermove`-based resize handles on `SidebarZone`/`BottomZone`, with sizes in a zustand slice (session-only, no localStorage). `ActivityBar` state machine expanded from 3 states (full/sliver/hot) to 4 (adds "over-zone" hover during drag). `RoomsPanel` is a full rewrite into a compact list with a socket-driven store, a header-action bar, and an MUI `Menu` per row.

**Tech Stack:** TypeScript, React 19, MUI 7, Zustand, dockview-react 5, Playwright (E2E). No new npm dependencies.

**Reference spec:** [`docs/superpowers/specs/2026-04-17-dockview-ui-fixes-design.md`](../specs/2026-04-17-dockview-ui-fixes-design.md).

---

## File Structure

| File | Action | Responsibility |
|---|---|---|
| `frontend/src/panels/registry.tsx` | Modify (add export) | Export size constants (min/max). |
| `frontend/src/stores/slices/activityBarSlice.ts` | Modify | Add `leftWidth`/`rightWidth`/`bottomHeight` + `setBarSize`. |
| `frontend/src/panels/SidebarZone.tsx` | Modify | Read width from slice, render drag handle, pointer-capture drag. |
| `frontend/src/panels/BottomZone.tsx` | Modify | Same for height. |
| `frontend/src/panels/ActivityBar.tsx` | Modify | Expand drag-state machine, hot 56 px zone, over-zone hint label. |
| `frontend/src/panels/RoomsPanel.tsx` | Rewrite | Compact list + search + header actions + per-row menu. |
| `frontend/src/panels/roomsHeaderActions.tsx` | Create | `+ New empty` and `⇪ Upload file` buttons for the panel header. |
| `frontend/src/panels/roomRowMenu.tsx` | Create | `⋮` MUI `Menu` with the six per-row actions. |
| `frontend/e2e/dockview-layout.spec.ts` | Modify (extend) | New E2E assertions per issue. |

No file deletions. No new npm dependencies. `/rooms/` page and `RoomManagementMenu.tsx` are untouched.

---

## Preflight

### Task 0: Verify baseline

Confirm the branch builds and the existing E2E suite passes before touching anything. This is the "nothing to blame me for" checkpoint.

**Files:** none.

- [ ] **Step 0.1: Ensure working tree is clean on `spec/dockview-ui-redesign`.**

```bash
git status --short
git rev-parse --abbrev-ref HEAD
```

Expected: HEAD is `spec/dockview-ui-redesign`, no uncommitted changes except untracked (`CLAUDE.md`, `ISSUES.md`, `frontend/.playwright-cli/`, etc.).

- [ ] **Step 0.2: Typecheck and lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

Expected: zero errors.

- [ ] **Step 0.3: Install frontend deps + run one existing dockview E2E to sanity-check.**

```bash
cd frontend && bun install
```

Start the servers in two terminals / background:

```bash
uv run zndraw --no-browser &
cd frontend && bun run dev &
```

Seed the test room once (see top of `frontend/e2e/dockview-layout.spec.ts`):

```bash
uv run zndraw-cli rooms create --room dockview-test 2>/dev/null || true
uv run python -c "
from zndraw import ZnDraw
import plotly.graph_objects as go
vis = ZnDraw(url='http://localhost:8000', room='dockview-test')
for key in ['energy']:
    if key not in vis.figures:
        vis.figures[key] = go.Figure(go.Scatter(x=[0,1,2], y=[0,1,4], mode='lines'))
"
```

Run one existing test:

```bash
cd frontend && bun x playwright test e2e/dockview-layout.spec.ts -g "default startup shows activity bars"
```

Expected: PASS.

If any of 0.1–0.3 fail: stop and reconcile before proceeding. Do not "just ignore" pre-existing failures — if they're pre-existing, commit them into the baseline notes here before moving on.

---

## Issue #1 — Resizable sidebar & bottom zones

### Task 1: Export size constants from `registry.tsx`

The zustand slice and the zones need to share `MIN`/`MAX`/default size constants. Put them in one place so future consumers import from there.

**Files:**
- Modify: `frontend/src/panels/registry.tsx` (append exports at the top-level, below the `PANELS` map).

- [ ] **Step 1.1: Open `frontend/src/panels/registry.tsx` and add the constants at the end of the file.**

Exact code to append (just before the final newline):

```typescript
// Sidebar / bottom-zone sizes (px). Session-only — stored in activityBarSlice.
export const SIDEBAR_DEFAULT_PX = 320;
export const SIDEBAR_MIN_PX = 200;
export const SIDEBAR_MAX_PX = 640;

export const BOTTOM_DEFAULT_PX = 260;
export const BOTTOM_MIN_PX = 120;
export const BOTTOM_MAX_PX = 560;
```

- [ ] **Step 1.2: Typecheck.**

```bash
cd frontend && bun --bun tsc --noEmit
```

Expected: zero errors.

- [ ] **Step 1.3: Commit.**

```bash
git add frontend/src/panels/registry.tsx
git commit -m "feat(panels): export sidebar/bottom size constants"
```

---

### Task 2: Add width/height state to `activityBarSlice`

Add three new fields + one setter. Keep the setter simple and clamp at the boundary.

**Files:**
- Modify: `frontend/src/stores/slices/activityBarSlice.ts`.

- [ ] **Step 2.1: Extend the `ActivityBarSlice` interface at the top of the file.**

Replace:

```typescript
export interface ActivityBarSlice {
	leftBarIcons: PanelId[];
	rightBarIcons: PanelId[];
	bottomBarIcons: PanelId[];
	activeLeft: PanelId | null;
	activeRight: PanelId | null;
	activeBottom: PanelId | null;

	moveIconToBar: (id: PanelId, bar: BarPosition, index?: number) => void;
	toggleActive: (bar: BarPosition, id: PanelId) => void;
	resetLayout: () => void;
}
```

With:

```typescript
export interface ActivityBarSlice {
	leftBarIcons: PanelId[];
	rightBarIcons: PanelId[];
	bottomBarIcons: PanelId[];
	activeLeft: PanelId | null;
	activeRight: PanelId | null;
	activeBottom: PanelId | null;

	leftWidth: number;
	rightWidth: number;
	bottomHeight: number;

	moveIconToBar: (id: PanelId, bar: BarPosition, index?: number) => void;
	toggleActive: (bar: BarPosition, id: PanelId) => void;
	resetLayout: () => void;
	setBarSize: (bar: BarPosition, px: number) => void;
}
```

- [ ] **Step 2.2: Import the new constants and seed initial state.**

At the top, add to the existing import block (replace the existing registry import):

```typescript
import {
	type BarPosition,
	BOTTOM_DEFAULT_PX,
	BOTTOM_MAX_PX,
	BOTTOM_MIN_PX,
	getDefaultsForBar,
	type PanelId,
	PANELS,
	SIDEBAR_DEFAULT_PX,
	SIDEBAR_MAX_PX,
	SIDEBAR_MIN_PX,
} from "../../panels/registry";
```

Replace `initialState()`:

```typescript
function initialState() {
	return {
		leftBarIcons: getDefaultsForBar("left"),
		rightBarIcons: getDefaultsForBar("right"),
		bottomBarIcons: getDefaultsForBar("bottom"),
		activeLeft: null as PanelId | null,
		activeRight: null as PanelId | null,
		activeBottom: null as PanelId | null,
		leftWidth: SIDEBAR_DEFAULT_PX,
		rightWidth: SIDEBAR_DEFAULT_PX,
		bottomHeight: BOTTOM_DEFAULT_PX,
	};
}
```

- [ ] **Step 2.3: Add the `setBarSize` action at the end of the slice body, just before the closing `});`.**

```typescript
	setBarSize: (bar, px) =>
		set(() => {
			if (bar === "bottom") {
				const clamped = Math.min(
					BOTTOM_MAX_PX,
					Math.max(BOTTOM_MIN_PX, px),
				);
				return { bottomHeight: clamped };
			}
			const clamped = Math.min(SIDEBAR_MAX_PX, Math.max(SIDEBAR_MIN_PX, px));
			return bar === "left"
				? { leftWidth: clamped }
				: { rightWidth: clamped };
		}),
```

- [ ] **Step 2.4: Typecheck.**

```bash
cd frontend && bun --bun tsc --noEmit
```

Expected: zero errors. If the slice is referenced anywhere that enumerates its keys (`Object.keys(state)` etc.), update those sites too.

- [ ] **Step 2.5: Commit.**

```bash
git add frontend/src/stores/slices/activityBarSlice.ts
git commit -m "feat(panels): add resizable sidebar/bottom state"
```

---

### Task 3: Sidebar zone drag handle

Replace the hard-coded `WIDTH = 320` with a zustand-backed width + render a drag handle on the inner edge.

**Files:**
- Modify: `frontend/src/panels/SidebarZone.tsx`.

- [ ] **Step 3.1: Replace the file body with the resizable version.**

Full file:

```tsx
import { Box } from "@mui/material";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import {
	type BarPosition,
	PANELS,
	SIDEBAR_MAX_PX,
	SIDEBAR_MIN_PX,
} from "./registry";

interface SidebarZoneProps {
	position: Exclude<BarPosition, "bottom">;
}

export function SidebarZone({ position }: SidebarZoneProps) {
	const active = useAppStore((s) =>
		position === "left" ? s.activeLeft : s.activeRight,
	);
	const width = useAppStore((s) =>
		position === "left" ? s.leftWidth : s.rightWidth,
	);
	const setBarSize = useAppStore((s) => s.setBarSize);
	const zoneRef = useRef<HTMLDivElement | null>(null);

	const onPointerDown = useCallback(
		(e: React.PointerEvent<HTMLDivElement>) => {
			e.preventDefault();
			(e.target as HTMLElement).setPointerCapture(e.pointerId);
			const startX = e.clientX;
			const startWidth = zoneRef.current?.getBoundingClientRect().width ?? width;

			const onMove = (ev: PointerEvent) => {
				const delta = ev.clientX - startX;
				const next = position === "left" ? startWidth + delta : startWidth - delta;
				setBarSize(position, next);
			};
			const onUp = (ev: PointerEvent) => {
				(e.target as HTMLElement).releasePointerCapture(ev.pointerId);
				window.removeEventListener("pointermove", onMove);
				window.removeEventListener("pointerup", onUp);
			};
			window.addEventListener("pointermove", onMove);
			window.addEventListener("pointerup", onUp);
		},
		[position, setBarSize, width],
	);

	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	// Clamp on render too, in case SSR-ish state had stale value.
	const safeWidth = Math.min(SIDEBAR_MAX_PX, Math.max(SIDEBAR_MIN_PX, width));

	return (
		<Box
			ref={zoneRef}
			data-testid={`sidebar-zone-${position}`}
			sx={{
				position: "relative",
				width: safeWidth,
				flexShrink: 0,
				borderColor: "divider",
				borderRight: position === "left" ? 1 : 0,
				borderLeft: position === "right" ? 1 : 0,
				display: "flex",
				flexDirection: "column",
				overflow: "hidden",
				bgcolor: "background.default",
			}}
		>
			<Component />
			<Box
				data-testid={`sidebar-resize-${position}`}
				onPointerDown={onPointerDown}
				sx={{
					position: "absolute",
					top: 0,
					bottom: 0,
					width: 8,
					// inner edge: right for left sidebar, left for right sidebar.
					[position === "left" ? "right" : "left"]: -4,
					cursor: "col-resize",
					zIndex: 1,
					"&:hover": { bgcolor: "action.hover" },
					touchAction: "none",
				}}
			/>
		</Box>
	);
}
```

- [ ] **Step 3.2: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

Expected: zero errors.

- [ ] **Step 3.3: Smoke test manually.**

Open `http://localhost:5173/rooms/dockview-test`, click Geometries icon, then drag the right edge of the sidebar outward. The sidebar should grow smoothly. Drag past `SIDEBAR_MAX_PX` — it should clamp.

- [ ] **Step 3.4: Commit.**

```bash
git add frontend/src/panels/SidebarZone.tsx
git commit -m "feat(panels): drag handle on sidebar inner edge"
```

---

### Task 4: Bottom zone drag handle

Analogous to Task 3 but for the top edge.

**Files:**
- Modify: `frontend/src/panels/BottomZone.tsx`.

- [ ] **Step 4.1: Replace the file body.**

Full file:

```tsx
import { Box } from "@mui/material";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { BOTTOM_MAX_PX, BOTTOM_MIN_PX, PANELS } from "./registry";

export function BottomZone() {
	const active = useAppStore((s) => s.activeBottom);
	const height = useAppStore((s) => s.bottomHeight);
	const setBarSize = useAppStore((s) => s.setBarSize);
	const zoneRef = useRef<HTMLDivElement | null>(null);

	const onPointerDown = useCallback(
		(e: React.PointerEvent<HTMLDivElement>) => {
			e.preventDefault();
			(e.target as HTMLElement).setPointerCapture(e.pointerId);
			const startY = e.clientY;
			const startHeight =
				zoneRef.current?.getBoundingClientRect().height ?? height;

			const onMove = (ev: PointerEvent) => {
				const delta = ev.clientY - startY;
				setBarSize("bottom", startHeight - delta); // dragging up grows the zone
			};
			const onUp = (ev: PointerEvent) => {
				(e.target as HTMLElement).releasePointerCapture(ev.pointerId);
				window.removeEventListener("pointermove", onMove);
				window.removeEventListener("pointerup", onUp);
			};
			window.addEventListener("pointermove", onMove);
			window.addEventListener("pointerup", onUp);
		},
		[setBarSize, height],
	);

	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	const safeHeight = Math.min(BOTTOM_MAX_PX, Math.max(BOTTOM_MIN_PX, height));

	return (
		<Box
			ref={zoneRef}
			data-testid="bottom-zone"
			sx={{
				position: "relative",
				height: safeHeight,
				flexShrink: 0,
				borderTop: 1,
				borderColor: "divider",
				display: "flex",
				flexDirection: "column",
				overflow: "hidden",
				bgcolor: "background.default",
			}}
		>
			<Component />
			<Box
				data-testid="bottom-resize"
				onPointerDown={onPointerDown}
				sx={{
					position: "absolute",
					left: 0,
					right: 0,
					top: -4,
					height: 8,
					cursor: "row-resize",
					zIndex: 1,
					"&:hover": { bgcolor: "action.hover" },
					touchAction: "none",
				}}
			/>
		</Box>
	);
}
```

- [ ] **Step 4.2: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

- [ ] **Step 4.3: Commit.**

```bash
git add frontend/src/panels/BottomZone.tsx
git commit -m "feat(panels): drag handle on bottom-zone top edge"
```

---

### Task 5: E2E for resize

Add four focused assertions.

**Files:**
- Modify: `frontend/e2e/dockview-layout.spec.ts`.

- [ ] **Step 5.1: Append the resize assertions to the end of `test.describe("dockview layout", () => { ... })`, just before the closing brace.**

```typescript
	test("sidebar resize: dragging the inner edge widens the zone", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-geometries").click();
		const zone = page.getByTestId("sidebar-zone-left");
		const handle = page.getByTestId("sidebar-resize-left");
		const before = await zone.boundingBox();
		if (!before) throw new Error("zone bounding box missing");

		const hBox = await handle.boundingBox();
		if (!hBox) throw new Error("handle bounding box missing");
		await page.mouse.move(hBox.x + hBox.width / 2, hBox.y + hBox.height / 2);
		await page.mouse.down();
		await page.mouse.move(hBox.x + 100, hBox.y + hBox.height / 2, { steps: 10 });
		await page.mouse.up();

		const after = await zone.boundingBox();
		if (!after) throw new Error("zone bounding box missing after");
		expect(after.width).toBeGreaterThan(before.width + 50);
	});

	test("sidebar resize clamps at SIDEBAR_MAX_PX", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-geometries").click();
		const zone = page.getByTestId("sidebar-zone-left");
		const handle = page.getByTestId("sidebar-resize-left");
		const hBox = await handle.boundingBox();
		if (!hBox) throw new Error("handle bounding box missing");

		// Drag way past the max; should clamp to 640 px.
		await page.mouse.move(hBox.x + hBox.width / 2, hBox.y + hBox.height / 2);
		await page.mouse.down();
		await page.mouse.move(hBox.x + 2000, hBox.y + hBox.height / 2, { steps: 20 });
		await page.mouse.up();

		const after = await zone.boundingBox();
		if (!after) throw new Error("zone bounding box missing after");
		expect(after.width).toBeLessThanOrEqual(641);
		expect(after.width).toBeGreaterThanOrEqual(639);
	});

	test("sidebar resize persists across close/open of the same panel", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-geometries").click();
		const handle = page.getByTestId("sidebar-resize-left");
		const hBox = await handle.boundingBox();
		if (!hBox) throw new Error("handle bounding box missing");
		await page.mouse.move(hBox.x + hBox.width / 2, hBox.y + hBox.height / 2);
		await page.mouse.down();
		await page.mouse.move(hBox.x + 150, hBox.y + hBox.height / 2, { steps: 10 });
		await page.mouse.up();
		const widthAfterResize = (
			await page.getByTestId("sidebar-zone-left").boundingBox()
		)?.width;
		if (widthAfterResize === undefined) throw new Error("no width");

		// Close and re-open.
		await page.getByTestId("activity-icon-geometries").click();
		await expect(page.getByTestId("sidebar-zone-left")).toBeHidden();
		await page.getByTestId("activity-icon-geometries").click();
		const widthAfterReopen = (
			await page.getByTestId("sidebar-zone-left").boundingBox()
		)?.width;
		expect(widthAfterReopen).toBe(widthAfterResize);
	});

	test("sidebar resize resets on page reload (session-only)", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-geometries").click();
		const handle = page.getByTestId("sidebar-resize-left");
		const hBox = await handle.boundingBox();
		if (!hBox) throw new Error("handle bounding box missing");
		await page.mouse.move(hBox.x + hBox.width / 2, hBox.y + hBox.height / 2);
		await page.mouse.down();
		await page.mouse.move(hBox.x + 150, hBox.y + hBox.height / 2, { steps: 10 });
		await page.mouse.up();

		await page.reload();
		await page.getByTestId("activity-icon-geometries").click();
		const box = await page.getByTestId("sidebar-zone-left").boundingBox();
		// Default is 320 ± 2 (border, etc.)
		expect(box?.width).toBeGreaterThanOrEqual(318);
		expect(box?.width).toBeLessThanOrEqual(322);
	});
```

- [ ] **Step 5.2: Run the new tests.**

```bash
cd frontend && bun x playwright test e2e/dockview-layout.spec.ts -g "sidebar resize"
```

Expected: 4 tests PASS.

- [ ] **Step 5.3: Commit.**

```bash
git add frontend/e2e/dockview-layout.spec.ts
git commit -m "test(panels): e2e for sidebar resize"
```

---

## Issue #4 — Wider snap zones during drag

### Task 6: Expand `ActivityBar` state machine (hot 56 px + over-zone)

Current file has two visual states: `full` / `sliver` / `hot`. We need a fourth — over-zone — and a wider 56 px hot zone. We also need to track whether the user's cursor is currently over *this* bar.

**Files:**
- Modify: `frontend/src/panels/ActivityBar.tsx`.

- [ ] **Step 6.1: Replace the `SLIVER_SX` and related constants near the top of the file.**

Find:

```typescript
const SLIVER_SX = {
	left: { width: 4 },
	right: { width: 4 },
	bottom: { height: 4, width: "100%" },
} as const;
```

Replace with:

```typescript
const SLIVER_IDLE_SX = {
	left: { width: 4 },
	right: { width: 4 },
	bottom: { height: 4, width: "100%" },
} as const;

// Hot zone (dragging) widens the empty bar to 56 px so it's an easy drop target.
const SLIVER_HOT_SX = {
	left: { width: 56 },
	right: { width: 56 },
	bottom: { height: 56, width: "100%" },
} as const;
```

- [ ] **Step 6.2: Extend the `SliverState` union.**

Find:

```typescript
type SliverState = "full" | "sliver" | "hot";
```

Replace with:

```typescript
type SliverState = "full" | "sliver" | "hot" | "over-zone";
```

- [ ] **Step 6.3: Track cursor-over-self with a local state and the bar's own `onDragOver` / `onDragLeave`.**

In the `ActivityBar` function body, **just after** `const [isDragActive, setIsDragActive] = useState(false);`, add:

```typescript
	const [isOverZone, setIsOverZone] = useState(false);
```

In the existing `useEffect` that listens for `dragstart`/`dragend`, extend the `onGlobalDragEnd` cleanup to also clear `isOverZone`:

```typescript
		const onGlobalDragEnd = () => {
			setIsDragActive(false);
			setIsOverZone(false);
		};
```

Make the `dragstart`/`dragend` listeners **passive** (`{ passive: true }`) — small perf win, free:

```typescript
		window.addEventListener("dragstart", onGlobalDragStart, { passive: true });
		window.addEventListener("dragend", onGlobalDragEnd, { passive: true });
```

- [ ] **Step 6.4: Replace the `state` derivation line (the ternary one near `icons.length === 0`).**

Find:

```typescript
	const state: SliverState =
		icons.length === 0 ? (isDragActive ? "hot" : "sliver") : "full";
```

Replace with:

```typescript
	const state: SliverState =
		icons.length === 0
			? isDragActive
				? isOverZone
					? "over-zone"
					: "hot"
				: "sliver"
			: "full";
```

- [ ] **Step 6.5: Update the empty-bar render branch (the `if (state !== "full")` block) to use the new states and zone size.**

Replace the entire `if (state !== "full") { return ( ... ); }` block with:

```tsx
	if (state !== "full") {
		const sizeSx = state === "sliver" ? SLIVER_IDLE_SX[position] : SLIVER_HOT_SX[position];
		const isHot = state === "hot";
		const isOver = state === "over-zone";
		const hintByPosition: Record<BarPosition, string> = {
			left: "Drop to dock left",
			right: "Drop to dock right",
			bottom: "Drop to dock bottom",
		};
		return (
			<Box
				data-testid={`activity-bar-${position}`}
				data-sliver-state={state}
				onDragOver={(e) => {
					onDragOver(e);
					if (isDragActive) setIsOverZone(true);
				}}
				onDragLeave={() => setIsOverZone(false)}
				onDrop={(e) => {
					setIsOverZone(false);
					onDropOnBar(e);
				}}
				sx={{
					display: "flex",
					flexShrink: 0,
					alignItems: "center",
					justifyContent: "center",
					position: "relative",
					...sizeSx,
					bgcolor: isOver
						? "rgba(25, 118, 210, 0.35)"
						: isHot
							? "rgba(25, 118, 210, 0.15)"
							: "transparent",
					borderStyle: isHot || isOver ? "solid" : "none",
					borderColor: "primary.main",
					borderWidth:
						isHot || isOver
							? position === "left"
								? "0 2px 0 0"
								: position === "right"
									? "0 0 0 2px"
									: "2px 0 0 0"
							: 0,
					// Dashed while hot (pulsing), solid when over-zone.
					borderStyle: isOver
						? "solid"
						: isHot
							? "dashed"
							: "none",
					animation: isHot ? `${pulse} 1s ease-in-out infinite` : "none",
				}}
			>
				{isOver && (
					<Box
						component="span"
						sx={{
							fontSize: 11,
							fontWeight: 500,
							color: "primary.main",
							whiteSpace: "nowrap",
							// Rotate for left/right so the text reads vertically along the bar.
							transform:
								position === "left"
									? "rotate(-90deg)"
									: position === "right"
										? "rotate(90deg)"
										: "none",
							pointerEvents: "none",
						}}
					>
						{hintByPosition[position]}
					</Box>
				)}
			</Box>
		);
	}
```

Note: the outer `sx` has two `borderStyle` lines — the second wins. That's fine; TypeScript will fold them. If the linter complains, delete the first `borderStyle` line.

- [ ] **Step 6.6: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

Fix any duplicate-key warnings. If the linter objects to the double `borderStyle`, keep only the second (conditional) definition.

- [ ] **Step 6.7: Commit.**

```bash
git add frontend/src/panels/ActivityBar.tsx
git commit -m "feat(panels): wider snap zones with over-zone state"
```

---

### Task 7: E2E for snap zones

Two tests: hot state widens, over-zone hover triggers the extra class.

**Files:**
- Modify: `frontend/e2e/dockview-layout.spec.ts`.

- [ ] **Step 7.1: Append to the existing `describe` block, after the resize tests.**

```typescript
	test("dragging widens empty bars to 56 px", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		const right = page.getByTestId("activity-bar-right");

		// Dispatch the synthetic dragstart (same pattern as the existing
		// "lights up empty bars as hot drop zones" test).
		await page.evaluate(() => {
			const dt = new DataTransfer();
			dt.setData("application/x-zndraw-panel-id", "selections");
			window.dispatchEvent(
				new DragEvent("dragstart", { dataTransfer: dt, bubbles: true }),
			);
		});

		await expect(right).toHaveAttribute("data-sliver-state", "hot");
		const box = await right.boundingBox();
		if (!box) throw new Error("right bar bounding box missing");
		expect(box.width).toBeGreaterThanOrEqual(54);
		expect(box.width).toBeLessThanOrEqual(58);

		await page.evaluate(() =>
			window.dispatchEvent(new DragEvent("dragend", { bubbles: true })),
		);
		await expect(right).toHaveAttribute("data-sliver-state", "sliver");
	});

	test("over-zone state appears when cursor enters a hot bar", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		const right = page.getByTestId("activity-bar-right");

		await page.evaluate(() => {
			const dt = new DataTransfer();
			dt.setData("application/x-zndraw-panel-id", "selections");
			window.dispatchEvent(
				new DragEvent("dragstart", { dataTransfer: dt, bubbles: true }),
			);
		});
		await expect(right).toHaveAttribute("data-sliver-state", "hot");

		// Dispatch dragover directly on the bar — Playwright's mouse.move
		// alone wouldn't carry the DataTransfer payload.
		await right.evaluate((el) => {
			const dt = new DataTransfer();
			dt.setData("application/x-zndraw-panel-id", "selections");
			el.dispatchEvent(
				new DragEvent("dragover", {
					dataTransfer: dt,
					bubbles: true,
					cancelable: true,
				}),
			);
		});
		await expect(right).toHaveAttribute("data-sliver-state", "over-zone");
		await expect(right.getByText(/Drop to dock right/i)).toBeVisible();

		await page.evaluate(() =>
			window.dispatchEvent(new DragEvent("dragend", { bubbles: true })),
		);
	});
```

- [ ] **Step 7.2: Run the new tests.**

```bash
cd frontend && bun x playwright test e2e/dockview-layout.spec.ts -g "widens empty|over-zone"
```

Expected: 2 tests PASS.

- [ ] **Step 7.3: Commit.**

```bash
git add frontend/e2e/dockview-layout.spec.ts
git commit -m "test(panels): e2e for wider snap + over-zone state"
```

---

## Issue #3 — Rooms panel feature parity

### Task 8: Header actions component

Two buttons: `+ New empty` and `⇪ Upload file`. Extracted into its own file because the panel is about to get busy.

**Files:**
- Create: `frontend/src/panels/roomsHeaderActions.tsx`.

- [ ] **Step 8.1: Create the file with exactly this content.**

```tsx
import AddIcon from "@mui/icons-material/Add";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import { Box, IconButton, Tooltip } from "@mui/material";
import { isAxiosError } from "axios";
import { useRef } from "react";
import { useNavigate } from "react-router-dom";
import { createRoom, uploadTrajectory } from "../myapi/client";
import { useAppStore } from "../store";

export function RoomsHeaderActions() {
	const navigate = useNavigate();
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const fileInputRef = useRef<HTMLInputElement | null>(null);

	const onNewEmpty = async () => {
		const id = crypto.randomUUID();
		try {
			await createRoom({ room_id: id });
			navigate(`/rooms/${id}`);
		} catch (err) {
			const detail = isAxiosError(err)
				? (err.response?.data?.detail ?? err.message)
				: "Failed to create room";
			showSnackbar(detail, "error");
		}
	};

	const onUploadClick = () => fileInputRef.current?.click();

	const onFileChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
		const files = e.target.files;
		if (!files || files.length === 0) return;
		const id = crypto.randomUUID();
		try {
			await createRoom({ room_id: id });
			for (const f of Array.from(files)) {
				await uploadTrajectory(id, f);
			}
			showSnackbar(`Room created with ${files.length} file(s)`, "success");
			navigate(`/rooms/${id}`);
		} catch (err) {
			const detail = isAxiosError(err)
				? (err.response?.data?.detail ?? err.message)
				: "Upload failed";
			showSnackbar(detail, "error");
		} finally {
			if (fileInputRef.current) fileInputRef.current.value = "";
		}
	};

	return (
		<Box sx={{ display: "flex", gap: 0.25 }}>
			<Tooltip title="New empty room">
				<IconButton
					size="small"
					data-testid="rooms-new-empty"
					onClick={onNewEmpty}
				>
					<AddIcon fontSize="small" />
				</IconButton>
			</Tooltip>
			<Tooltip title="Upload file → new room">
				<IconButton
					size="small"
					data-testid="rooms-upload"
					onClick={onUploadClick}
				>
					<UploadFileIcon fontSize="small" />
				</IconButton>
			</Tooltip>
			<input
				type="file"
				ref={fileInputRef}
				style={{ display: "none" }}
				onChange={onFileChange}
				data-testid="rooms-upload-input"
			/>
		</Box>
	);
}
```

- [ ] **Step 8.2: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

- [ ] **Step 8.3: Commit.**

```bash
git add frontend/src/panels/roomsHeaderActions.tsx
git commit -m "feat(rooms-panel): header actions component"
```

---

### Task 9: Per-row `⋮` menu component

Six items: set template / duplicate / lock / copy link / download frames / delete (disabled with tooltip for now — no backend endpoint exists).

**Files:**
- Create: `frontend/src/panels/roomRowMenu.tsx`.

- [ ] **Step 9.1: Create the file.**

```tsx
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import DeleteIcon from "@mui/icons-material/Delete";
import DownloadIcon from "@mui/icons-material/Download";
import DuplicateIcon from "@mui/icons-material/FileCopy";
import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import StarIcon from "@mui/icons-material/Star";
import {
	IconButton,
	ListItemIcon,
	ListItemText,
	Menu,
	MenuItem,
	Tooltip,
} from "@mui/material";
import MoreVertIcon from "@mui/icons-material/MoreVert";
import { useState } from "react";
import DuplicateRoomDialog from "../components/DuplicateRoomDialog";
import {
	type Room,
	downloadFrames,
	setDefaultRoom,
	updateRoom,
} from "../myapi/client";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";

interface Props {
	room: Room;
}

export function RoomRowMenu({ room }: Props) {
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const rooms = useRoomsStore((s) => s.roomsArray);
	const [anchor, setAnchor] = useState<HTMLElement | null>(null);
	const [duplicateOpen, setDuplicateOpen] = useState(false);
	const open = Boolean(anchor);

	const onSetTemplate = async () => {
		setAnchor(null);
		try {
			await setDefaultRoom(room.is_default ? null : room.id);
			useRoomsStore
				.getState()
				.updateRoom(room.id, { is_default: !room.is_default });
			showSnackbar(
				room.is_default ? "Template cleared" : "Set as template",
				"success",
			);
		} catch {
			showSnackbar("Failed to update template", "error");
		}
	};

	const onToggleLock = async () => {
		setAnchor(null);
		try {
			await updateRoom(room.id, { locked: !room.locked });
			showSnackbar(room.locked ? "Room unlocked" : "Room locked", "success");
		} catch {
			showSnackbar("Failed to update lock", "error");
		}
	};

	const onDuplicate = () => {
		setAnchor(null);
		setDuplicateOpen(true);
	};

	const onCopyLink = async () => {
		setAnchor(null);
		await navigator.clipboard.writeText(
			`${window.location.origin}/rooms/${room.id}`,
		);
		showSnackbar("Link copied", "success");
	};

	const onDownload = () => {
		setAnchor(null);
		downloadFrames({ roomId: room.id });
		showSnackbar("Downloading all frames", "success");
	};

	return (
		<>
			<IconButton
				size="small"
				data-testid={`room-row-menu-${room.id}`}
				onClick={(e) => setAnchor(e.currentTarget)}
				onMouseDown={(e) => e.stopPropagation()}
				aria-label="Room actions"
			>
				<MoreVertIcon fontSize="small" />
			</IconButton>
			<Menu anchorEl={anchor} open={open} onClose={() => setAnchor(null)}>
				<MenuItem onClick={onSetTemplate}>
					<ListItemIcon>
						{room.is_default ? <StarBorderIcon /> : <StarIcon />}
					</ListItemIcon>
					<ListItemText>
						{room.is_default ? "Remove template" : "Set as template"}
					</ListItemText>
				</MenuItem>
				<MenuItem onClick={onDuplicate}>
					<ListItemIcon>
						<DuplicateIcon />
					</ListItemIcon>
					<ListItemText>Duplicate room</ListItemText>
				</MenuItem>
				<MenuItem onClick={onToggleLock}>
					<ListItemIcon>
						{room.locked ? <LockOpenIcon /> : <LockIcon />}
					</ListItemIcon>
					<ListItemText>{room.locked ? "Unlock" : "Lock"}</ListItemText>
				</MenuItem>
				<MenuItem onClick={onCopyLink}>
					<ListItemIcon>
						<ContentCopyIcon />
					</ListItemIcon>
					<ListItemText>Copy link</ListItemText>
				</MenuItem>
				<MenuItem onClick={onDownload}>
					<ListItemIcon>
						<DownloadIcon />
					</ListItemIcon>
					<ListItemText>Download frames</ListItemText>
				</MenuItem>
				<Tooltip title="Deleting rooms is not yet supported by the backend.">
					<span>
						<MenuItem disabled>
							<ListItemIcon>
								<DeleteIcon />
							</ListItemIcon>
							<ListItemText>Delete</ListItemText>
						</MenuItem>
					</span>
				</Tooltip>
			</Menu>
			<DuplicateRoomDialog
				open={duplicateOpen}
				sourceRoomId={room.id}
				sourceDescription={room.description ?? room.id}
				existingRoomIds={rooms.map((r) => r.id)}
				onClose={() => setDuplicateOpen(false)}
			/>
		</>
	);
}
```

- [ ] **Step 9.2: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

If `DuplicateIcon` import fails (`FileCopy` not exported at that path), use `ContentCopyIcon` as a second reference and rename the import to distinguish from the copy-link icon.

- [ ] **Step 9.3: Commit.**

```bash
git add frontend/src/panels/roomRowMenu.tsx
git commit -m "feat(rooms-panel): per-row action menu"
```

---

### Task 10: Rewrite `RoomsPanel.tsx`

Header actions + search + list with description/frame-count/⋮. Row click still switches rooms with the leave-room cascade.

**Files:**
- Modify: `frontend/src/panels/RoomsPanel.tsx` (full rewrite).

- [ ] **Step 10.1: Replace the entire file.**

```tsx
import SearchIcon from "@mui/icons-material/Search";
import {
	Box,
	InputAdornment,
	List,
	ListItemButton,
	ListItemText,
	TextField,
	Typography,
} from "@mui/material";
import { useEffect, useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";
import { getDockviewApi } from "./DockviewLayout";
import { RoomsHeaderActions } from "./roomsHeaderActions";
import { RoomRowMenu } from "./roomRowMenu";

/**
 * VS-Code-style compact rooms list with search, header actions, and
 * a per-row actions menu. Reactive via the existing roomsStore.
 */
export function RoomsPanel() {
	const rooms = useRoomsStore((s) => s.roomsArray);
	const loading = useRoomsStore((s) => s.loading);
	const fetchRooms = useRoomsStore((s) => s.fetchRooms);
	const currentRoomId = useAppStore((s) => s.roomId);
	const leaveRoom = useLeaveRoom({ api: getDockviewApi() });
	const navigate = useNavigate();
	const [query, setQuery] = useState("");

	useEffect(() => {
		fetchRooms();
	}, [fetchRooms]);

	const filtered = useMemo(() => {
		const q = query.trim().toLowerCase();
		if (!q) return rooms;
		return rooms.filter(
			(r) =>
				r.id.toLowerCase().includes(q) ||
				(r.description && r.description.toLowerCase().includes(q)),
		);
	}, [rooms, query]);

	const switchToRoom = async (id: string) => {
		if (id === currentRoomId) return;
		await leaveRoom({ skipConfirm: true });
		navigate(`/rooms/${id}`);
	};

	return (
		<Box
			data-testid="rooms-panel"
			sx={{ display: "flex", flexDirection: "column", height: "100%" }}
		>
			<Box
				sx={{
					display: "flex",
					alignItems: "center",
					justifyContent: "space-between",
					px: 1.5,
					pt: 1,
					pb: 0.5,
				}}
			>
				<Typography
					variant="overline"
					sx={{ color: "text.secondary", letterSpacing: "0.08em" }}
				>
					Rooms
				</Typography>
				<RoomsHeaderActions />
			</Box>
			<Box sx={{ px: 1.5, pb: 1 }}>
				<TextField
					size="small"
					fullWidth
					value={query}
					onChange={(e) => setQuery(e.target.value)}
					placeholder="Search rooms…"
					data-testid="rooms-search"
					slotProps={{
						input: {
							startAdornment: (
								<InputAdornment position="start">
									<SearchIcon fontSize="small" />
								</InputAdornment>
							),
						},
					}}
				/>
			</Box>
			{loading && (
				<Box sx={{ px: 2 }}>
					<Typography variant="body2" color="text.secondary">
						Loading rooms…
					</Typography>
				</Box>
			)}
			{!loading && filtered.length === 0 && (
				<Box sx={{ px: 2, pt: 1 }}>
					<Typography variant="body2" color="text.secondary">
						{query
							? "No rooms match the search."
							: "No rooms available — create one above."}
					</Typography>
				</Box>
			)}
			<List dense sx={{ flexGrow: 1, overflow: "auto", pt: 0 }}>
				{filtered.map((r) => (
					<ListItemButton
						key={r.id}
						data-testid={`rooms-row-${r.id}`}
						selected={r.id === currentRoomId}
						onClick={() => switchToRoom(r.id)}
						sx={{
							alignItems: "center",
							pl: 1.5,
							pr: 0.5,
							borderLeft: 2,
							borderColor:
								r.id === currentRoomId ? "primary.main" : "transparent",
						}}
					>
						<ListItemText
							primary={r.description ?? r.id}
							secondary={`${r.id.slice(0, 8)} · ${r.frame_count} frame${
								r.frame_count === 1 ? "" : "s"
							}`}
							primaryTypographyProps={{
								variant: "body2",
								fontWeight: 500,
								noWrap: true,
							}}
							secondaryTypographyProps={{
								variant: "caption",
								color: "text.secondary",
							}}
						/>
						<RoomRowMenu room={r} />
					</ListItemButton>
				))}
			</List>
		</Box>
	);
}
```

- [ ] **Step 10.2: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

- [ ] **Step 10.3: Smoke test manually.** Open the rooms panel; confirm list renders with description + id-prefix · frame count, `⋮` menu opens, search filters.

- [ ] **Step 10.4: Commit.**

```bash
git add frontend/src/panels/RoomsPanel.tsx
git commit -m "feat(rooms-panel): VS-Code-style list with search + menu"
```

---

### Task 10b: Drag-and-drop trajectory onto the rooms panel

Spec calls for dropping a trajectory file onto the panel body to create a new room and upload — same as the `/rooms/` page. Wrap the existing hook into `RoomsPanel`.

**Files:**
- Modify: `frontend/src/panels/RoomsPanel.tsx`.

- [ ] **Step 10b.1: Import the hook + a local `handleFiles` matching `pages/roomList.tsx`.**

At the top of `RoomsPanel.tsx`, add these imports alongside the existing ones:

```tsx
import { useCallback } from "react";
import { isAxiosError } from "axios";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import { createRoom, uploadTrajectory } from "../myapi/client";
```

Inside the `RoomsPanel` component body, right after `const [query, setQuery] = useState("");`, add:

```tsx
	const showSnackbar = useAppStore((s) => s.showSnackbar);

	const handleFiles = useCallback(
		async (files: File[]) => {
			const newRoomId = crypto.randomUUID();
			try {
				await createRoom({ room_id: newRoomId });
				for (const file of files) {
					await uploadTrajectory(newRoomId, file);
				}
				showSnackbar(`Room created with ${files.length} file(s)`, "success");
				navigate(`/rooms/${newRoomId}`);
			} catch (error) {
				const detail = isAxiosError(error)
					? (error.response?.data?.detail ?? error.message)
					: "Upload failed";
				showSnackbar(detail, "error");
			}
		},
		[navigate, showSnackbar],
	);

	const {
		isDragging,
		handleDragOver,
		handleDragEnter,
		handleDragLeave,
		handleDrop,
	} = useDragAndDrop(handleFiles);
```

- [ ] **Step 10b.2: Wire the handlers onto the outer `<Box data-testid="rooms-panel">` and render a drop overlay.**

Replace the outer `<Box data-testid="rooms-panel" ... >` opening tag with:

```tsx
		<Box
			data-testid="rooms-panel"
			onDragOver={handleDragOver}
			onDragEnter={handleDragEnter}
			onDragLeave={handleDragLeave}
			onDrop={handleDrop}
			sx={{ display: "flex", flexDirection: "column", height: "100%", position: "relative" }}
		>
```

Inside the box, as the last child before the closing `</Box>`, add:

```tsx
			{isDragging && (
				<Box
					data-testid="rooms-drop-overlay"
					sx={{
						position: "absolute",
						inset: 0,
						bgcolor: "rgba(25, 118, 210, 0.2)",
						border: 2,
						borderStyle: "dashed",
						borderColor: "primary.main",
						display: "flex",
						alignItems: "center",
						justifyContent: "center",
						pointerEvents: "none",
						zIndex: 2,
					}}
				>
					<Typography variant="body2" color="primary.main" fontWeight={500}>
						Drop to create a new room
					</Typography>
				</Box>
			)}
```

- [ ] **Step 10b.3: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

- [ ] **Step 10b.4: Commit.**

```bash
git add frontend/src/panels/RoomsPanel.tsx
git commit -m "feat(rooms-panel): drag-drop trajectory to create new room"
```

---

### Task 11: E2E for rooms panel

Three focused assertions: header actions exist, search filters, `⋮` menu opens.

**Files:**
- Modify: `frontend/e2e/dockview-layout.spec.ts`.

- [ ] **Step 11.1: Append to the `describe` block.**

```typescript
	test("rooms panel header has new-empty + upload buttons", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-rooms").click();
		await expect(page.getByTestId("rooms-new-empty")).toBeVisible();
		await expect(page.getByTestId("rooms-upload")).toBeVisible();
	});

	test("rooms panel search filters the list", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-rooms").click();
		const rows = page.locator('[data-testid^="rooms-row-"]');
		const initialCount = await rows.count();
		expect(initialCount).toBeGreaterThan(0);

		await page.getByTestId("rooms-search").fill("zz-no-such-room-zz");
		await expect(rows).toHaveCount(0);
		await expect(
			page.getByText(/No rooms match the search/i),
		).toBeVisible();
	});

	test("rooms panel shows drop overlay when dragging a file in", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-rooms").click();
		const panel = page.getByTestId("rooms-panel");
		// Synthetic dragenter with a files payload — we don't need a real file
		// since we only verify the overlay fires on the drag state.
		await panel.evaluate((el) => {
			const dt = new DataTransfer();
			el.dispatchEvent(
				new DragEvent("dragenter", {
					dataTransfer: dt,
					bubbles: true,
					cancelable: true,
				}),
			);
		});
		await expect(page.getByTestId("rooms-drop-overlay")).toBeVisible();
	});

	test("rooms panel row menu opens and exposes template/lock items", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-rooms").click();
		const menuBtn = page.getByTestId(`room-row-menu-${ROOM}`);
		await expect(menuBtn).toBeVisible();
		await menuBtn.click();
		await expect(page.getByRole("menuitem", { name: /Set as template|Remove template/i })).toBeVisible();
		await expect(page.getByRole("menuitem", { name: /^Lock$|^Unlock$/ })).toBeVisible();
		await expect(page.getByRole("menuitem", { name: /Duplicate room/i })).toBeVisible();
		// Delete is disabled until a backend endpoint lands.
		const deleteItem = page.getByRole("menuitem", { name: /Delete/i });
		await expect(deleteItem).toHaveAttribute("aria-disabled", "true");
	});
```

- [ ] **Step 11.2: Run the new tests.**

```bash
cd frontend && bun x playwright test e2e/dockview-layout.spec.ts -g "rooms panel"
```

Expected: 3 tests PASS.

- [ ] **Step 11.3: Commit.**

```bash
git add frontend/e2e/dockview-layout.spec.ts
git commit -m "test(rooms-panel): e2e for header + search + menu"
```

---

## Final verification

### Task 12: Full suite + manual smoke

Make sure nothing else broke.

**Files:** none.

- [ ] **Step 12.1: Typecheck + lint.**

```bash
cd frontend && bun run lint && bun --bun tsc --noEmit
```

Expected: zero errors.

- [ ] **Step 12.2: Run the full dockview E2E suite.**

```bash
cd frontend && bun x playwright test e2e/dockview-layout.spec.ts
```

Expected: all tests PASS (existing ones + the 9 new ones added in tasks 5, 7, 11).

- [ ] **Step 12.3: Manual smoke via a live browser session.**

Start servers:

```bash
uv run zndraw --no-browser &
cd frontend && bun run dev &
```

Walk through this checklist by hand at `http://localhost:5173/rooms/dockview-test`:

- [ ] Open Geometries panel → drag the right edge → sidebar grows.
- [ ] Close + re-open Geometries → width preserved.
- [ ] Reload page → Geometries opens at 320 px.
- [ ] Drag a panel icon off the left bar — right bar widens to 56 px with dashed border + pulse; cursor over right bar → solid border + "Drop to dock right" label; drop → panel moves.
- [ ] Open Rooms panel → see compact list, search filters, `⋮` menu opens with 6 items (Delete disabled). Click `+ New empty` → new room created + navigated.

- [ ] **Step 12.4: No new commit here unless manual smoke turns up an issue.**

---

## Out-of-scope reminders

- **Don't touch `pages/roomList.tsx`.** The `/rooms/` page stays exactly as it is.
- **Don't add a `DELETE /v1/rooms/{id}` endpoint** or a `deleteRoom` client method in this plan — that's a separate backend PR. The menu item is disabled.
- **Don't add localStorage persistence** for sidebar widths — session-only by design.
- **Don't rewrite panels for perf** — that's spec C.
- **Don't register default filesystem providers** — that's spec A.
