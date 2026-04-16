# Dockview UI Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace ZnDraw's fixed sidebar + floating window layout with a VS Code-style two-tier architecture: three activity bars (Left/Right/Bottom) driving single-panel-per-bar zones for tools, and a `dockview-react`-managed editor area for views (3D viewer + plot tabs). Full replacement — no feature flag, no coexistence, no dead code.

**Architecture:** Activity bars + sidebar/bottom zones are plain React + a Zustand `activityBarSlice` (dockview cannot enforce "single panel per bar"). The editor area is a single `DockviewReact` instance hosting the viewer and any open plot tabs. A frontend `PANELS` registry (`frontend/src/panels/registry.ts`) is the single source of truth; adding a new panel = one entry + one component file. The 3D Viewer's lifecycle is bound to the current room — closing it triggers `leaveCurrentRoom()` (cascade closes plot tabs, clears URL). Switching rooms via the Rooms panel has identical cascade semantics.

**Tech Stack:** React 19, TypeScript, Zustand, MUI v7, `dockview-react` (NEW, ~v5.2.x), Vite 7. Removes `react-rnd`. No backend changes.

**Spec:** [`docs/superpowers/specs/2026-04-02-dockview-ui-redesign-design.md`](../specs/2026-04-02-dockview-ui-redesign-design.md)

---

## Dev & Validation Workflow

**Every task that changes UI code MUST be validated via playwright-cli against a live dev server before commit.** Do not skip this.

### Start the stack (once, keep running across tasks)

Open three terminals (or use tmux / background runs):

```bash
# Terminal 1 — backend on :8000
cd /Users/fzills/tools/zndraw-fastapi
uv run zndraw
```

```bash
# Terminal 2 — frontend dev (HMR, :5173, proxies /v1 /api /socket.io → :8000)
cd /Users/fzills/tools/zndraw-fastapi/frontend
bun install
bun run dev
```

```bash
# Terminal 3 — CLI room for test data (once per session)
cd /Users/fzills/tools/zndraw-fastapi
uv run zndraw-cli rooms list                         # verify server up
uv run zndraw-cli rooms create --id dockview-test    # or use existing
export ZNDRAW_ROOM=dockview-test
uv run zndraw-cli frames extend test_water.xyz       # load sample data
```

### Seed plots for panels that need them (run once, reuse)

```bash
uv run zndraw-cli extensions list                    # discover analysis extensions
uv run zndraw-cli extensions run "@internal:selections:All" --wait
# Pick an analysis extension that produces a figure, e.g.:
uv run zndraw-cli extensions run "@internal:analysis:Distance" --wait
uv run zndraw-cli figures list                       # confirm figures exist
```

If no analysis extension produces a figure, create one directly:

```bash
uv run python -c "
from zndraw import ZnDraw
import plotly.graph_objects as go
vis = ZnDraw(room='dockview-test')
for key in ['energy', 'rdf', 'rmsd', 'distance']:
    fig = go.Figure(go.Scatter(x=[0,1,2,3], y=[0,1,4,9], mode='lines'))
    fig.update_layout(title=key)
    vis.figures[key] = fig
"
```

### Validate via playwright-cli (after every UI task)

Use `:5173` (the Vite dev server) so HMR picks up edits. Using `:8000` would require rebuilding.

```bash
playwright-cli open http://localhost:5173/rooms/dockview-test
playwright-cli snapshot --filename=<task-name>-before.yaml
# ... click/drag/type to exercise the change ...
playwright-cli snapshot --filename=<task-name>-after.yaml
playwright-cli screenshot --filename=<task-name>.png
playwright-cli console                 # MUST be empty — no React errors, no warnings from our code
playwright-cli close
```

**Rules for playwright-cli validation:**

1. Always take a snapshot or screenshot showing the change worked. Save to `/tmp/dockview-validation/` or inline as artifacts.
2. Always run `playwright-cli console` — browser errors block commit.
3. Prefer `playwright-cli snapshot` (returns a YAML accessibility tree with `e<N>` refs) over screenshot for interactive asserts — you can verify specific nodes exist/don't exist.
4. Use `playwright-cli -s=dockview` for a reusable session across tasks, closing only at the end.

### Validate via zndraw CLI (for data-dependent changes)

Use the `/zndraw` skill when a task touches figures, rooms, filesystem, or extensions. Examples:

```bash
uv run zndraw-cli figures set rdf --data '{"type":"plotly","data":"{\"data\":[{\"x\":[0,1],\"y\":[0,1]}]}"}'
# → watch browser: plot should auto-open as a dockview tab

uv run zndraw-cli figures delete rdf
# → watch browser: the rdf tab should auto-close

uv run zndraw-cli rooms create --id dockview-switch-test
# → click the new room in the Rooms panel: should cascade-close plots + switch room + update URL
```

### Per-commit checks (run before every `git commit`)

```bash
cd frontend && bun run lint && bunx tsc --noEmit
```

The project does not ship unit tests (only Playwright E2E under `frontend/e2e/`). Tasks that need new tests add them there. Tasks that only affect UI rendering add a playwright-cli validation rather than a `.spec.ts` unless the behavior is worth protecting from regressions.

### E2E test runner

```bash
cd frontend
bunx playwright install chromium     # one-time
ZNDRAW_URL=http://localhost:8000 bunx playwright test e2e/dockview-activity-bar.spec.ts
```

The E2E tests run against the **built** static frontend on `:8000`, not the Vite dev server. Run `cd frontend && bun run build` before E2E runs (this rebuilds into `src/zndraw/static/`).

---

## File Structure

New files (all under `frontend/src/panels/` unless noted):

| File | Responsibility |
|---|---|
| `panels/registry.ts` | `PanelId` union, `PANELS` record, types. Single source of truth. |
| `panels/ActivityBar.tsx` | Icon strip for one bar. Click to toggle, drag to move between bars. |
| `panels/SidebarZone.tsx` | Container for Left/Right zones. Mounts `PANELS[active].component`. |
| `panels/BottomZone.tsx` | Same for the bottom zone (horizontal layout). |
| `panels/DockviewLayout.tsx` | Root `DockviewReact`, registers view components, welcome state, onUnhandled/onDidDrop. |
| `panels/ViewerView.tsx` | Wraps `MyScene`; hooks dockview lifecycle + owns `onClose → leaveCurrentRoom()`. |
| `panels/PlotView.tsx` | Dockview-native plot tab; migrates FigureWindow rendering (no Rnd). |
| `panels/plotViewFactory.ts` | `createPlotView(figureKey)` for dynamic instantiation. |
| `panels/PlotsBrowserPanel.tsx` | Left-bar panel listing available + open plots. Click opens, drag drops into dockview. |
| `panels/RoomsPanel.tsx` | Left-bar panel listing rooms. Row click = cascade + `history.pushState`. |
| `panels/FilesystemPanel.tsx` | Left-bar panel. Unwraps `FilesystemBrowser` from page-level chrome. |
| `panels/ChatPanel.tsx` | Right-bar panel. Unwraps chat content from `react-rnd`. |
| `stores/slices/activityBarSlice.ts` | Zustand slice: `leftBarIcons`, `rightBarIcons`, `bottomBarIcons`, `activeLeft/Right/Bottom`, `moveIconToBar`, `toggleActive`, `resetLayout`. |
| `hooks/useLeaveRoom.ts` | F1 cascade: confirm → close plots → close viewer → disconnect → `history.pushState('/')`. |
| `panels/resetLayoutMenu.tsx` | "Reset layout" item added to AppBar "..." menu; calls slice `resetLayout()` + clears dockview to viewer-only. |

Modified files:

| File | Change |
|---|---|
| `pages/landingPage.tsx` | Replace sidebar + Canvas + WindowManager layout with ActivityBar / SidebarZone / DockviewLayout / BottomZone root. Remove AddPlot, chat toggle, RoomManagementMenu from AppBar. |
| `App.tsx` | `/rooms/:roomId/files` redirects to `/rooms/:roomId?panel=filesystem` (opens filesystem panel on load). |
| `stores/slices/uiSlice.ts` | Drop `chatOpen`, `setChatOpen`. Keep `chatUnreadCount` + `incrementChatUnread` + `resetChatUnread` (bar icon shows the badge now). |
| `hooks/socketHandlers/figureHandlers.ts` | Replace `useWindowManagerStore` calls with dockview API (`openPlotTab(key)` / `closePlotTab(key)`). |
| `hooks/socketHandlers/chatHandlers.ts` | `incrementChatUnread` check gates on "is chat panel the active panel in its bar", not `chatOpen`. |
| `frontend/package.json` | Add `dockview-react`; remove `react-rnd`. |

Removed files (final verification task must find zero references):

- `src/components/SideBar.tsx`
- `src/components/PrimaryDrawer.tsx`
- `src/components/WindowManager.tsx`
- `src/components/FigureWindow.tsx`
- `src/components/AddPlotButton.tsx`
- `src/components/ChatWindow.tsx`
- `src/stores/windowManagerStore.ts`
- `src/formStore.ts`
- `src/pages/filesystemBrowser.tsx`

Unchanged: all `src/components/three/**`, `PlaybackSlice`, `SceneSlice`, `LockSlice`, `ConnectionSlice`, socket transport, backend.

---

## Task 0: Baseline & Setup

**Files:**
- Modify: `frontend/package.json`
- Modify: `frontend/src/panels/` (directory create)
- Test: — (baseline snapshot only)

- [ ] **Step 1: Verify dev environment**

```bash
cd /Users/fzills/tools/zndraw-fastapi
uv run zndraw-cli rooms list
cd frontend && bun install && bun run lint && bunx tsc --noEmit
```

Expected: server reachable, lint/type both exit 0.

- [ ] **Step 2: Take pre-migration baseline screenshot**

Start backend + frontend dev server (see "Dev & Validation Workflow" above). Then:

```bash
playwright-cli -s=dockview open http://localhost:5173/rooms/dockview-test
playwright-cli -s=dockview snapshot --filename=baseline-layout.yaml
playwright-cli -s=dockview screenshot --filename=baseline-layout.png
playwright-cli -s=dockview close
```

Expected: snapshot shows existing AppBar, SideBar icons, MyScene canvas, FrameProgressBar. Save both files to `/tmp/dockview-validation/`.

- [ ] **Step 3: Install dockview-react**

```bash
cd frontend
bun add dockview-react@^5.2.0
```

- [ ] **Step 4: Verify install**

```bash
cd frontend
grep '"dockview-react"' package.json
ls node_modules/dockview-react/dist
```

Expected: package listed in `package.json` dependencies; `dist` folder present.

- [ ] **Step 5: Create panels directory skeleton**

```bash
mkdir -p /Users/fzills/tools/zndraw-fastapi/frontend/src/panels
mkdir -p /Users/fzills/tools/zndraw-fastapi/frontend/src/hooks
```

- [ ] **Step 6: Typecheck passes with dockview installed**

```bash
cd frontend && bunx tsc --noEmit
```

Expected: 0 errors.

- [ ] **Step 7: Commit**

```bash
git add frontend/package.json frontend/bun.lockb
git commit -m "chore(frontend): add dockview-react dependency"
```

---

## Task 1: Panel Registry Types & Stubs

**Files:**
- Create: `frontend/src/panels/registry.ts`
- Create: `frontend/src/panels/stubs.tsx` (placeholder components — real ones replace in later tasks)

- [ ] **Step 1: Write placeholder components**

Create `frontend/src/panels/stubs.tsx`:

```tsx
import { Box, Typography } from "@mui/material";
import type { IDockviewPanelProps } from "dockview-react";

const Stub = ({ label }: { label: string }) => (
	<Box sx={{ p: 2 }}>
		<Typography variant="h6">{label}</Typography>
		<Typography variant="body2" color="text.secondary">
			Placeholder — implementation in a later task.
		</Typography>
	</Box>
);

export const StubSelections = () => <Stub label="Selections" />;
export const StubModifiers = () => <Stub label="Modifiers" />;
export const StubAnalysis = () => <Stub label="Analysis" />;
export const StubGeometries = () => <Stub label="Geometries" />;
export const StubPlotsBrowser = () => <Stub label="Plots Browser" />;
export const StubRooms = () => <Stub label="Rooms" />;
export const StubFilesystem = () => <Stub label="Files" />;
export const StubChat = () => <Stub label="Chat" />;

export const StubViewerView = (_props: IDockviewPanelProps) => (
	<Stub label="3D Viewer (stub)" />
);
```

- [ ] **Step 2: Write the registry**

Create `frontend/src/panels/registry.ts`:

```ts
import Analytics from "@mui/icons-material/Analytics";
import Build from "@mui/icons-material/Build";
import Category from "@mui/icons-material/Category";
import Chat from "@mui/icons-material/Chat";
import FilterCenterFocus from "@mui/icons-material/FilterCenterFocus";
import Folder from "@mui/icons-material/Folder";
import MeetingRoom from "@mui/icons-material/MeetingRoom";
import ShowChart from "@mui/icons-material/ShowChart";
import type { SvgIconComponent } from "@mui/icons-material";
import type { IDockviewPanelProps } from "dockview-react";
import type { ComponentType } from "react";
import {
	StubAnalysis,
	StubChat,
	StubFilesystem,
	StubGeometries,
	StubModifiers,
	StubPlotsBrowser,
	StubRooms,
	StubSelections,
	StubViewerView,
} from "./stubs";

export type PanelId =
	| "selections"
	| "modifiers"
	| "analysis"
	| "geometries"
	| "plots-browser"
	| "rooms"
	| "filesystem"
	| "chat"
	| "viewer";

export type BarPosition = "left" | "right" | "bottom";

type PanelDefault = {
	bar: BarPosition | "editor";
	active?: boolean;
	order?: number;
};

type ToolPanelDef = {
	kind: "tool";
	icon: SvgIconComponent;
	label: string;
	component: ComponentType;
	default: PanelDefault;
};

type ViewPanelDef = {
	kind: "view";
	title: string;
	component: ComponentType<IDockviewPanelProps>;
	closable: boolean;
	cascadeOnClose?: (string | RegExp)[];
	onClose?: () => void;
};

export type PanelDef = ToolPanelDef | ViewPanelDef;

export const PANELS: Record<PanelId, PanelDef> = {
	selections: {
		kind: "tool",
		icon: FilterCenterFocus,
		label: "Selections",
		component: StubSelections,
		default: { bar: "left", order: 0 },
	},
	modifiers: {
		kind: "tool",
		icon: Build,
		label: "Modifiers",
		component: StubModifiers,
		default: { bar: "left", order: 1 },
	},
	analysis: {
		kind: "tool",
		icon: Analytics,
		label: "Analysis",
		component: StubAnalysis,
		default: { bar: "left", order: 2 },
	},
	geometries: {
		kind: "tool",
		icon: Category,
		label: "Geometries",
		component: StubGeometries,
		default: { bar: "left", order: 3 },
	},
	"plots-browser": {
		kind: "tool",
		icon: ShowChart,
		label: "Plots",
		component: StubPlotsBrowser,
		default: { bar: "left", order: 4 },
	},
	rooms: {
		kind: "tool",
		icon: MeetingRoom,
		label: "Rooms",
		component: StubRooms,
		default: { bar: "left", order: 5 },
	},
	filesystem: {
		kind: "tool",
		icon: Folder,
		label: "Files",
		component: StubFilesystem,
		default: { bar: "left", order: 6 },
	},
	chat: {
		kind: "tool",
		icon: Chat,
		label: "Chat",
		component: StubChat,
		default: { bar: "right", order: 0 },
	},
	viewer: {
		kind: "view",
		title: "3D Viewer",
		component: StubViewerView,
		closable: true,
		cascadeOnClose: [/^plot-.+$/],
		// onClose is assigned in a later task (requires useLeaveRoom hook)
	},
};

export const TOOL_PANEL_IDS = (Object.entries(PANELS) as [PanelId, PanelDef][])
	.filter(([, def]) => def.kind === "tool")
	.map(([id]) => id);

export function getDefaultsForBar(bar: BarPosition): PanelId[] {
	return (Object.entries(PANELS) as [PanelId, PanelDef][])
		.filter(([, def]) => def.kind === "tool" && def.default.bar === bar)
		.sort(
			([, a], [, b]) =>
				(a.kind === "tool" ? (a.default.order ?? 0) : 0) -
				(b.kind === "tool" ? (b.default.order ?? 0) : 0),
		)
		.map(([id]) => id);
}
```

- [ ] **Step 3: Typecheck**

```bash
cd frontend && bunx tsc --noEmit
```

Expected: 0 errors. If `StubViewerView` prop typing fails, check the import path of `IDockviewPanelProps` against the installed version.

- [ ] **Step 4: Commit**

```bash
git add frontend/src/panels/registry.ts frontend/src/panels/stubs.tsx
git commit -m "feat(panels): add panel registry types and stubs"
```

---

## Task 2: Activity Bar Slice (Zustand)

**Files:**
- Create: `frontend/src/stores/slices/activityBarSlice.ts`
- Modify: `frontend/src/store.tsx`

- [ ] **Step 1: Write the slice**

Create `frontend/src/stores/slices/activityBarSlice.ts`:

```ts
import type { StateCreator } from "zustand";
import type { AppState } from "../../store";
import {
	type BarPosition,
	getDefaultsForBar,
	type PanelId,
	PANELS,
} from "../../panels/registry";

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

function initialState() {
	return {
		leftBarIcons: getDefaultsForBar("left"),
		rightBarIcons: getDefaultsForBar("right"),
		bottomBarIcons: getDefaultsForBar("bottom"),
		activeLeft: null as PanelId | null,
		activeRight: null as PanelId | null,
		activeBottom: null as PanelId | null,
	};
}

const BAR_KEY: Record<BarPosition, "leftBarIcons" | "rightBarIcons" | "bottomBarIcons"> = {
	left: "leftBarIcons",
	right: "rightBarIcons",
	bottom: "bottomBarIcons",
};

const ACTIVE_KEY: Record<BarPosition, "activeLeft" | "activeRight" | "activeBottom"> = {
	left: "activeLeft",
	right: "activeRight",
	bottom: "activeBottom",
};

export const createActivityBarSlice: StateCreator<
	AppState,
	[],
	[],
	ActivityBarSlice
> = (set) => ({
	...initialState(),

	moveIconToBar: (id, bar, index) =>
		set((state) => {
			const sourceBar: BarPosition | null =
				state.leftBarIcons.includes(id)
					? "left"
					: state.rightBarIcons.includes(id)
						? "right"
						: state.bottomBarIcons.includes(id)
							? "bottom"
							: null;
			if (!sourceBar) return {};
			if (PANELS[id].kind !== "tool") return {};

			const patch: Partial<ActivityBarSlice> = {};

			// Remove from source bar
			const sourceKey = BAR_KEY[sourceBar];
			const sourceList = state[sourceKey].filter((x) => x !== id);
			patch[sourceKey] = sourceList as PanelId[] & never;

			// Clear active in source bar if this icon was active
			if (state[ACTIVE_KEY[sourceBar]] === id) {
				patch[ACTIVE_KEY[sourceBar]] = null as never;
			}

			// Insert into target bar (same or different)
			const targetKey = BAR_KEY[bar];
			const currentTarget = sourceBar === bar ? sourceList : state[targetKey];
			const nextTarget = [...currentTarget];
			const insertAt =
				index === undefined || index < 0 || index > nextTarget.length
					? nextTarget.length
					: index;
			nextTarget.splice(insertAt, 0, id);
			patch[targetKey] = nextTarget as PanelId[] & never;

			return patch as ActivityBarSlice;
		}),

	toggleActive: (bar, id) =>
		set((state) => {
			const key = ACTIVE_KEY[bar];
			const current = state[key];
			return { [key]: current === id ? null : id } as Partial<ActivityBarSlice>;
		}),

	resetLayout: () => set(initialState()),
});
```

- [ ] **Step 2: Wire into root store**

Edit `frontend/src/store.tsx`:

```tsx
// Add import near other slice imports
import {
	type ActivityBarSlice,
	createActivityBarSlice,
} from "./stores/slices/activityBarSlice";

// Extend AppState type
export type AppState = ConnectionSlice &
	PlaybackSlice &
	SceneSlice &
	LockSlice &
	UISlice &
	ActivityBarSlice;

// In useAppStore create(...)
export const useAppStore = create<AppState>((...a) => ({
	...createConnectionSlice(...a),
	...createPlaybackSlice(...a),
	...createSceneSlice(...a),
	...createLockSlice(...a),
	...createUISlice(...a),
	...createActivityBarSlice(...a),
}));
```

- [ ] **Step 3: Typecheck**

```bash
cd frontend && bunx tsc --noEmit
```

Expected: 0 errors. Typescript will flag any collision between `resetLayout` in other slices — there should be none; if `set` typing gets noisy in the slice, cast the patch once (as shown above).

- [ ] **Step 4: Validate initial state from browser console**

```bash
playwright-cli -s=dockview open http://localhost:5173/rooms/dockview-test
playwright-cli -s=dockview eval "window.__ZUSTAND_DEBUG__ = true; 'ok'"   # no-op; just warms page
playwright-cli -s=dockview eval "JSON.stringify(Object.keys(document))"   # sanity check
playwright-cli -s=dockview console
```

Expected: no new console errors. Slice is not yet rendered — only state is wired.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/stores/slices/activityBarSlice.ts frontend/src/store.tsx
git commit -m "feat(stores): add activity bar slice"
```

---

## Task 3: ActivityBar Component

**Files:**
- Create: `frontend/src/panels/ActivityBar.tsx`

- [ ] **Step 1: Write the component**

Create `frontend/src/panels/ActivityBar.tsx`:

```tsx
import { Box, IconButton, Tooltip } from "@mui/material";
import { useCallback } from "react";
import { useShallow } from "zustand/react/shallow";
import { useAppStore } from "../store";
import { type BarPosition, type PanelId, PANELS } from "./registry";

const DRAG_MIME = "application/x-zndraw-panel-id";

const BAR_SX = {
	left: { flexDirection: "column", width: 48, borderRight: 1 },
	right: { flexDirection: "column", width: 48, borderLeft: 1 },
	bottom: { flexDirection: "row", height: 40, borderTop: 1, width: "100%" },
} as const;

const ACTIVE_INDICATOR: Record<
	BarPosition,
	(color: string) => Record<string, string>
> = {
	left: (c) => ({ borderLeft: `2px solid ${c}` }),
	right: (c) => ({ borderRight: `2px solid ${c}` }),
	bottom: (c) => ({ borderTop: `2px solid ${c}` }),
};

interface ActivityBarProps {
	position: BarPosition;
}

export function ActivityBar({ position }: ActivityBarProps) {
	const icons = useAppStore(
		useShallow((s) =>
			position === "left"
				? s.leftBarIcons
				: position === "right"
					? s.rightBarIcons
					: s.bottomBarIcons,
		),
	);
	const active = useAppStore((s) =>
		position === "left"
			? s.activeLeft
			: position === "right"
				? s.activeRight
				: s.activeBottom,
	);
	const toggleActive = useAppStore((s) => s.toggleActive);
	const moveIconToBar = useAppStore((s) => s.moveIconToBar);

	const onDragStart = useCallback(
		(e: React.DragEvent, id: PanelId) => {
			e.dataTransfer.setData(DRAG_MIME, id);
			e.dataTransfer.effectAllowed = "move";
		},
		[],
	);

	const onDragOver = useCallback((e: React.DragEvent) => {
		if (e.dataTransfer.types.includes(DRAG_MIME)) {
			e.preventDefault();
			e.dataTransfer.dropEffect = "move";
		}
	}, []);

	const onDropOnBar = useCallback(
		(e: React.DragEvent) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			if (!id) return;
			e.preventDefault();
			moveIconToBar(id as PanelId, position);
		},
		[moveIconToBar, position],
	);

	const onDropOnIcon = useCallback(
		(e: React.DragEvent, overIdx: number) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			if (!id) return;
			e.preventDefault();
			e.stopPropagation();
			moveIconToBar(id as PanelId, position, overIdx);
		},
		[moveIconToBar, position],
	);

	return (
		<Box
			data-testid={`activity-bar-${position}`}
			onDragOver={onDragOver}
			onDrop={onDropOnBar}
			sx={{
				display: "flex",
				borderColor: "divider",
				bgcolor: "background.paper",
				alignItems: "center",
				...BAR_SX[position],
			}}
		>
			{icons.map((id, idx) => {
				const def = PANELS[id];
				if (def.kind !== "tool") return null;
				const Icon = def.icon;
				const isActive = active === id;
				return (
					<Tooltip key={id} title={def.label} placement={position === "left" ? "right" : position === "right" ? "left" : "top"}>
						<IconButton
							data-testid={`activity-icon-${id}`}
							draggable
							onDragStart={(e) => onDragStart(e, id)}
							onDragOver={onDragOver}
							onDrop={(e) => onDropOnIcon(e, idx)}
							onClick={() => toggleActive(position, id)}
							color={isActive ? "primary" : "default"}
							sx={{
								borderRadius: 0,
								m: position === "bottom" ? "0 2px" : "2px 0",
								...(isActive
									? ACTIVE_INDICATOR[position]("var(--mui-palette-primary-main)")
									: {}),
							}}
						>
							<Icon />
						</IconButton>
					</Tooltip>
				);
			})}
		</Box>
	);
}
```

- [ ] **Step 2: Typecheck + lint**

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

Expected: 0 errors.

- [ ] **Step 3: Commit (no render yet; mounted in Task 12)**

```bash
git add frontend/src/panels/ActivityBar.tsx
git commit -m "feat(panels): activity bar with click toggle and drag-between-bars"
```

---

## Task 4: SidebarZone and BottomZone

**Files:**
- Create: `frontend/src/panels/SidebarZone.tsx`
- Create: `frontend/src/panels/BottomZone.tsx`

- [ ] **Step 1: Write `SidebarZone.tsx`**

Create `frontend/src/panels/SidebarZone.tsx`:

```tsx
import { Box } from "@mui/material";
import { useAppStore } from "../store";
import type { BarPosition } from "./registry";
import { PANELS } from "./registry";

interface SidebarZoneProps {
	position: Exclude<BarPosition, "bottom">;
}

const WIDTH = 320;

export function SidebarZone({ position }: SidebarZoneProps) {
	const active = useAppStore((s) =>
		position === "left" ? s.activeLeft : s.activeRight,
	);
	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	return (
		<Box
			data-testid={`sidebar-zone-${position}`}
			sx={{
				width: WIDTH,
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
		</Box>
	);
}
```

- [ ] **Step 2: Write `BottomZone.tsx`**

Create `frontend/src/panels/BottomZone.tsx`:

```tsx
import { Box } from "@mui/material";
import { useAppStore } from "../store";
import { PANELS } from "./registry";

const HEIGHT = 260;

export function BottomZone() {
	const active = useAppStore((s) => s.activeBottom);
	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	return (
		<Box
			data-testid="bottom-zone"
			sx={{
				height: HEIGHT,
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
		</Box>
	);
}
```

- [ ] **Step 3: Typecheck**

```bash
cd frontend && bunx tsc --noEmit
```

Expected: 0 errors.

- [ ] **Step 4: Commit**

```bash
git add frontend/src/panels/SidebarZone.tsx frontend/src/panels/BottomZone.tsx
git commit -m "feat(panels): sidebar and bottom zone containers"
```

---

## Task 5: DockviewLayout with Viewer-Only Initial State

**Files:**
- Create: `frontend/src/panels/DockviewLayout.tsx`
- Create: `frontend/src/panels/ViewerView.tsx`

- [ ] **Step 1: Write `ViewerView.tsx`**

Create `frontend/src/panels/ViewerView.tsx`:

```tsx
import { Box } from "@mui/material";
import type { IDockviewPanelProps } from "dockview-react";
import { useEffect, useRef, useState } from "react";
import { MyScene } from "../components/Canvas";

export function ViewerView(props: IDockviewPanelProps) {
	const { api } = props;
	const [remountKey, setRemountKey] = useState(0);
	const previousLocation = useRef(api.location.type);

	useEffect(() => {
		const disposable = api.onDidLocationChange((e) => {
			if (e.location.type !== previousLocation.current) {
				previousLocation.current = e.location.type;
				setRemountKey((k) => k + 1);
			}
		});
		return () => disposable.dispose();
	}, [api]);

	return (
		<Box
			key={remountKey}
			data-testid="viewer-view"
			sx={{ position: "relative", width: "100%", height: "100%", overflow: "hidden" }}
		>
			<MyScene />
		</Box>
	);
}
```

- [ ] **Step 2: Write `DockviewLayout.tsx`**

Create `frontend/src/panels/DockviewLayout.tsx`:

```tsx
import { Box, Typography } from "@mui/material";
import type { DockviewApi, DockviewReadyEvent } from "dockview-react";
import { DockviewReact } from "dockview-react";
import "dockview-react/dist/styles/dockview.css";
import { useCallback, useEffect, useRef, useState } from "react";
import { PANELS } from "./registry";
import { ViewerView } from "./ViewerView";

const components = {
	viewer: ViewerView,
	// plotView registered in Task 7
};

const DRAG_MIME_PLOT = "application/x-zndraw-plot-key";

function addViewerPanel(api: DockviewApi) {
	api.addPanel({
		id: "viewer",
		component: "viewer",
		title: PANELS.viewer.kind === "view" ? PANELS.viewer.title : "Viewer",
		params: { canClose: true },
	});
}

export function DockviewLayout() {
	const apiRef = useRef<DockviewApi | null>(null);
	const [isEmpty, setIsEmpty] = useState(false);

	const onReady = useCallback((event: DockviewReadyEvent) => {
		apiRef.current = event.api;
		addViewerPanel(event.api);
		setIsEmpty(event.api.panels.length === 0);

		event.api.onUnhandledDragOverEvent((e) => {
			// accept only panel-drag-from-plots-browser via onDidDrop path (Task 8)
			const dt = (e.nativeEvent as DragEvent)?.dataTransfer;
			if (dt?.types.includes(DRAG_MIME_PLOT)) e.accept();
		});

		event.api.onDidRemovePanel(() => {
			setIsEmpty(event.api.panels.length === 0);
		});
		event.api.onDidAddPanel(() => {
			setIsEmpty(event.api.panels.length === 0);
		});
	}, []);

	// Re-init viewer when returning from empty state (Task 6 wires room-join → re-add)
	useEffect(() => {
		const api = apiRef.current;
		if (!api) return;
		// No-op here; actual re-init triggered by useLeaveRoom / room switch.
	}, []);

	return (
		<Box sx={{ flexGrow: 1, position: "relative", minWidth: 0, minHeight: 0 }}>
			<DockviewReact
				className="dockview-theme-light"
				onReady={onReady}
				components={components}
				floatingGroupBounds="boundedWithinViewport"
			/>
			{isEmpty && (
				<Box
					data-testid="dockview-welcome"
					sx={{
						position: "absolute",
						inset: 0,
						display: "flex",
						alignItems: "center",
						justifyContent: "center",
						pointerEvents: "none",
					}}
				>
					<Typography variant="body1" color="text.secondary">
						No room selected. Open the Rooms panel to pick one.
					</Typography>
				</Box>
			)}
		</Box>
	);
}
```

- [ ] **Step 3: Export module for external use**

Create `frontend/src/panels/index.ts`:

```ts
export { ActivityBar } from "./ActivityBar";
export { BottomZone } from "./BottomZone";
export { DockviewLayout } from "./DockviewLayout";
export { SidebarZone } from "./SidebarZone";
export { PANELS } from "./registry";
export type { PanelId, BarPosition } from "./registry";
```

- [ ] **Step 4: Typecheck**

```bash
cd frontend && bunx tsc --noEmit
```

Expected: 0 errors. If dockview-react's CSS import path differs, check `node_modules/dockview-react/dist/styles/dockview.css` exists; otherwise use the path from its `package.json` `exports` field.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/panels/DockviewLayout.tsx frontend/src/panels/ViewerView.tsx frontend/src/panels/index.ts
git commit -m "feat(panels): dockview layout with viewer-only initial state"
```

---

## Task 6: useLeaveRoom + Viewer onClose Wiring

**Files:**
- Create: `frontend/src/hooks/useLeaveRoom.ts`
- Modify: `frontend/src/panels/ViewerView.tsx`
- Modify: `frontend/src/panels/DockviewLayout.tsx`

- [ ] **Step 1: Write the hook**

Create `frontend/src/hooks/useLeaveRoom.ts`:

```ts
import { useCallback } from "react";
import type { DockviewApi } from "dockview-react";
import { useAppStore } from "../store";

interface LeaveRoomOptions {
	api: DockviewApi | null;
	showConfirm?: (message: string) => Promise<boolean>;
}

export function useLeaveRoom({ api, showConfirm }: LeaveRoomOptions) {
	const showSnackbar = useAppStore((s) => s.showSnackbar);

	return useCallback(
		async (opts?: { skipConfirm?: boolean }) => {
			if (!api) return;
			const plotPanels = api.panels.filter((p) => p.id.startsWith("plot-"));
			const hasPlots = plotPanels.length > 0;

			if (hasPlots && !opts?.skipConfirm) {
				const msg = `Leave room? ${plotPanels.length} plot(s) will close.`;
				const ok = showConfirm ? await showConfirm(msg) : window.confirm(msg);
				if (!ok) return;
			}

			// Close plot panels first (cascade)
			for (const panel of plotPanels) {
				panel.api.close();
			}
			// Close viewer if present
			const viewer = api.getPanel("viewer");
			if (viewer) viewer.api.close();

			// URL + state cleanup
			window.history.pushState({}, "", "/");
			showSnackbar("Left room", "info");
			// Note: full disconnect / store reset is handled by the socket
			// manager watching the URL change (existing behavior).
		},
		[api, showConfirm, showSnackbar],
	);
}
```

- [ ] **Step 2: Expose api to ViewerView onClose**

Edit `frontend/src/panels/DockviewLayout.tsx` — expose api via module-level ref so ViewerView can trigger cascade:

Add near top of file:

```ts
let sharedApi: DockviewApi | null = null;
export function getDockviewApi(): DockviewApi | null {
	return sharedApi;
}
```

Inside `onReady`, replace `apiRef.current = event.api;` with:

```ts
apiRef.current = event.api;
sharedApi = event.api;
```

- [ ] **Step 3: Hook cascade into ViewerView close**

Edit `frontend/src/panels/ViewerView.tsx` — append:

```tsx
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import { getDockviewApi } from "./DockviewLayout";
```

And inside the component, add:

```tsx
const leaveRoom = useLeaveRoom({ api: getDockviewApi() });
useEffect(() => {
	const d = api.onWillFocus(() => {}); // no-op; ensures api is alive
	const onBefore = api.onDidClose?.(() => {});
	const closeHandler = api.onDidDispose?.(() => leaveRoom({ skipConfirm: false }));
	return () => {
		d.dispose();
		onBefore?.dispose();
		closeHandler?.dispose();
	};
}, [api, leaveRoom]);
```

If the dockview version you installed exposes only `onWillClose`, use it instead (cancellable); otherwise use `onDidDispose` as above. Confirm by running typecheck.

**Preferred implementation (cancellable):** if `onWillClose` exists, use it with `e.preventDefault()` on cancel, then trigger `leaveRoom` on confirm and call `api.close()` manually.

- [ ] **Step 4: Wire registry `onClose`**

Edit `frontend/src/panels/registry.ts` — the `viewer` entry already has `cascadeOnClose`. Leave the concrete close cascade to `ViewerView` (cleaner — it has access to hooks). No change to registry needed.

- [ ] **Step 5: Typecheck**

```bash
cd frontend && bunx tsc --noEmit
```

Expected: 0 errors. Resolve any unused symbols.

- [ ] **Step 6: Playwright validation (after Task 12 mounts the layout — skip here; validate this task's plumbing via typecheck only)**

This task has no user-visible change yet; the test is deferred to Task 12. A pre-wiring sanity check:

```bash
cd frontend && bun run lint
```

- [ ] **Step 7: Commit**

```bash
git add frontend/src/hooks/useLeaveRoom.ts frontend/src/panels/DockviewLayout.tsx frontend/src/panels/ViewerView.tsx
git commit -m "feat(panels): useLeaveRoom hook + viewer close cascade"
```

---

## Task 7: PlotView + Factory + Socket Handler Migration

**Files:**
- Create: `frontend/src/panels/PlotView.tsx`
- Create: `frontend/src/panels/plotViewFactory.ts`
- Modify: `frontend/src/panels/DockviewLayout.tsx`
- Modify: `frontend/src/hooks/socketHandlers/figureHandlers.ts`

- [ ] **Step 1: Write `PlotView.tsx`** — migrate FigureWindow rendering, drop Rnd wrapper

Create `frontend/src/panels/PlotView.tsx`. Keep the plotly initialization, marker tracks, and step/selection-marker updates from `FigureWindow.tsx` — strip everything related to `<Rnd>`, `WindowTemplate`, autocomplete header, zIndex, drag-handle. The entry-point is the IDockviewPanelProps `params.figureKey`.

```tsx
import ErrorIcon from "@mui/icons-material/Error";
import { Box, CircularProgress, Typography } from "@mui/material";
import type { IDockviewPanelProps } from "dockview-react";
import { useEffect, useMemo, useRef } from "react";
import { useFigure } from "../hooks/useFigures";
import { useAppStore } from "../store";

interface PlotViewParams {
	figureKey: string;
}

export function PlotView(props: IDockviewPanelProps<PlotViewParams>) {
	const { figureKey } = props.params;
	const plotContainer = useRef<HTMLDivElement>(null);

	const {
		data: figureResponse,
		isLoading,
		isError,
		error,
	} = useFigure(figureKey);

	const currentFrame = useAppStore((s) => s.currentFrame);
	const frame_selection = useAppStore((s) => s.frame_selection);
	const selections = useAppStore((s) => s.selections);

	const plotlyJson = useMemo(() => {
		if (figureResponse?.figure?.type === "plotly" && figureResponse.figure.data) {
			try {
				return JSON.parse(figureResponse.figure.data);
			} catch {
				return null;
			}
		}
		return null;
	}, [figureResponse]);

	// Plot init — copy the effect from FigureWindow.tsx lines ~640-985
	// but render into plotContainer.current without any Rnd bounds.
	useEffect(() => {
		let cancelled = false;
		if (!plotlyJson || !plotContainer.current) return;
		(async () => {
			const PlotlyJS = await import("plotly.js-dist-min");
			if (cancelled || !plotContainer.current) return;
			await PlotlyJS.newPlot(
				plotContainer.current,
				plotlyJson.data,
				{ ...plotlyJson.layout, autosize: true },
				{ responsive: true, displaylogo: false },
			);
		})();
		return () => {
			cancelled = true;
			if (plotContainer.current) {
				import("plotly.js-dist-min").then((PlotlyJS) =>
					plotContainer.current && PlotlyJS.purge(plotContainer.current),
				);
			}
		};
	}, [plotlyJson]);

	// TODO(task-7-followup): port step marker + selection marker effects
	// from FigureWindow.tsx 660-985 unchanged (they reference currentFrame,
	// frame_selection, selections). Insert them here verbatim, replacing
	// window refs with plotContainer.current.
	useEffect(() => {
		void currentFrame;
		void frame_selection;
		void selections;
	}, [currentFrame, frame_selection, selections]);

	if (isLoading) {
		return (
			<Box sx={{ display: "flex", alignItems: "center", justifyContent: "center", height: "100%" }}>
				<CircularProgress />
			</Box>
		);
	}
	if (isError) {
		return (
			<Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "center", height: "100%", gap: 1 }}>
				<ErrorIcon color="error" />
				<Typography>{error?.message ?? "Failed to load figure"}</Typography>
			</Box>
		);
	}
	return (
		<Box
			ref={plotContainer}
			data-testid={`plot-view-${figureKey}`}
			sx={{ width: "100%", height: "100%" }}
		/>
	);
}
```

**IMPORTANT:** replace the `TODO(task-7-followup)` block by porting the marker/track effects from `FigureWindow.tsx` verbatim — that code is battle-tested. Use `FigureWindow.tsx` lines 456–985 as the source. Remove only Rnd/autocomplete-specific state.

- [ ] **Step 2: Write the factory**

Create `frontend/src/panels/plotViewFactory.ts`:

```ts
import type { DockviewApi } from "dockview-react";

export function plotPanelId(figureKey: string): string {
	return `plot-${figureKey}`;
}

export function openPlotTab(
	api: DockviewApi,
	figureKey: string,
	position?: { referenceGroupId?: string; direction?: "right" | "below" | "within" },
): void {
	const id = plotPanelId(figureKey);
	const existing = api.getPanel(id);
	if (existing) {
		existing.api.setActive();
		return;
	}
	api.addPanel({
		id,
		component: "plotView",
		title: figureKey,
		params: { figureKey },
		position: position?.referenceGroupId
			? { referenceGroup: api.getGroup(position.referenceGroupId), direction: position.direction ?? "within" }
			: undefined,
	});
}

export function closePlotTab(api: DockviewApi, figureKey: string): void {
	const panel = api.getPanel(plotPanelId(figureKey));
	panel?.api.close();
}
```

- [ ] **Step 3: Register the component in `DockviewLayout.tsx`**

Edit `frontend/src/panels/DockviewLayout.tsx`:

```tsx
import { PlotView } from "./PlotView";
import { openPlotTab } from "./plotViewFactory";

const components = {
	viewer: ViewerView,
	plotView: PlotView,
};
```

And add an `onDidDrop` handler:

```tsx
function onDidDrop(event: DockviewDropEvent) {
	const dt = (event.nativeEvent as DragEvent)?.dataTransfer;
	const key = dt?.getData(DRAG_MIME_PLOT);
	if (!key) return;
	event.api.addPanel({
		id: `plot-${key}`,
		component: "plotView",
		title: key,
		params: { figureKey: key },
		position: {
			referenceGroup: event.group ?? undefined,
			direction: event.position ?? "within",
		},
	});
}

// In the JSX:
<DockviewReact
	/* ...existing... */
	onDidDrop={onDidDrop}
/>
```

- [ ] **Step 4: Migrate figure handlers**

Edit `frontend/src/hooks/socketHandlers/figureHandlers.ts`. Replace all `useWindowManagerStore.getState().openWindow/closeWindow/openWindows` with:

```ts
import { getDockviewApi } from "../../panels/DockviewLayout";
import { closePlotTab, openPlotTab, plotPanelId } from "../../panels/plotViewFactory";

// inside onFiguresInvalidate:
if (operation === "delete") {
	// ... invalidate queries ...
	const api = getDockviewApi();
	if (api) closePlotTab(api, data.key);
} else if (operation === "set") {
	// ... invalidate queries ...
	const api = getDockviewApi();
	if (api && !api.getPanel(plotPanelId(data.key))) {
		openPlotTab(api, data.key);
	}
}
```

Remove the `MAX_AUTO_OPEN_WINDOWS` cap (spec: "auto-open as a tab in the active editor group" — no cap).

- [ ] **Step 5: Validate via zndraw-cli + playwright** (only fully validates once Task 12 mounts layout; still do a typecheck here)

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

Expected: 0 errors.

- [ ] **Step 6: Commit**

```bash
git add frontend/src/panels/PlotView.tsx frontend/src/panels/plotViewFactory.ts frontend/src/panels/DockviewLayout.tsx frontend/src/hooks/socketHandlers/figureHandlers.ts
git commit -m "feat(panels): PlotView and figure handler migration to dockview"
```

---

## Task 8: PlotsBrowserPanel

**Files:**
- Create: `frontend/src/panels/PlotsBrowserPanel.tsx`
- Modify: `frontend/src/panels/registry.ts`

- [ ] **Step 1: Write the component**

Create `frontend/src/panels/PlotsBrowserPanel.tsx`:

```tsx
import CloseIcon from "@mui/icons-material/Close";
import {
	Box,
	Divider,
	IconButton,
	List,
	ListItem,
	ListItemButton,
	ListItemText,
	Typography,
} from "@mui/material";
import { useMemo } from "react";
import { useFigureList } from "../hooks/useFigures";
import { getDockviewApi } from "./DockviewLayout";
import { closePlotTab, openPlotTab, plotPanelId } from "./plotViewFactory";

const DRAG_MIME_PLOT = "application/x-zndraw-plot-key";

export function PlotsBrowserPanel() {
	const { data, isLoading } = useFigureList();
	const allKeys = useMemo(() => data?.items ?? [], [data]);

	const openKeys = useMemo(() => {
		const api = getDockviewApi();
		if (!api) return new Set<string>();
		return new Set(
			api.panels
				.filter((p) => p.id.startsWith("plot-"))
				.map((p) => p.id.slice("plot-".length)),
		);
	}, [data]);

	const onRowClick = (key: string) => {
		const api = getDockviewApi();
		if (!api) return;
		openPlotTab(api, key);
	};

	const onRowClose = (key: string) => {
		const api = getDockviewApi();
		if (!api) return;
		closePlotTab(api, key);
	};

	const onDragStart = (e: React.DragEvent, key: string) => {
		e.dataTransfer.setData(DRAG_MIME_PLOT, key);
		e.dataTransfer.effectAllowed = "copy";
	};

	if (isLoading) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography>Loading plots…</Typography>
			</Box>
		);
	}

	if (allKeys.length === 0) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography color="text.secondary">No plots available</Typography>
			</Box>
		);
	}

	const openKeysArr = allKeys.filter((k) => openKeys.has(k));

	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Typography variant="overline" sx={{ px: 2, pt: 1 }}>
				Available Plots
			</Typography>
			<List dense sx={{ flexGrow: 1, overflow: "auto" }}>
				{allKeys.map((k) => (
					<ListItem
						key={k}
						disablePadding
						draggable
						onDragStart={(e) => onDragStart(e, k)}
					>
						<ListItemButton onClick={() => onRowClick(k)}>
							<Box
								sx={{
									width: 8,
									height: 8,
									borderRadius: "50%",
									bgcolor: openKeys.has(k) ? "primary.main" : "transparent",
									border: 1,
									borderColor: "primary.main",
									mr: 1.5,
								}}
							/>
							<ListItemText primary={k} />
						</ListItemButton>
					</ListItem>
				))}
			</List>
			{openKeysArr.length > 0 && (
				<>
					<Divider />
					<Typography variant="overline" sx={{ px: 2, pt: 1 }}>
						Currently Open ({openKeysArr.length})
					</Typography>
					<List dense>
						{openKeysArr.map((k) => (
							<ListItem
								key={k}
								secondaryAction={
									<IconButton
										edge="end"
										size="small"
										onClick={() => onRowClose(k)}
									>
										<CloseIcon fontSize="small" />
									</IconButton>
								}
							>
								<ListItemText primary={k} />
							</ListItem>
						))}
					</List>
				</>
			)}
		</Box>
	);
}
```

- [ ] **Step 2: Swap stub for real component in registry**

Edit `frontend/src/panels/registry.ts`:

```tsx
import { PlotsBrowserPanel } from "./PlotsBrowserPanel";
// ...
"plots-browser": {
	kind: "tool",
	icon: ShowChart,
	label: "Plots",
	component: PlotsBrowserPanel,
	default: { bar: "left", order: 4 },
},
```

- [ ] **Step 3: Typecheck + lint**

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

- [ ] **Step 4: Commit**

```bash
git add frontend/src/panels/PlotsBrowserPanel.tsx frontend/src/panels/registry.ts
git commit -m "feat(panels): plots browser panel with click + drag"
```

---

## Task 9: RoomsPanel

**Files:**
- Create: `frontend/src/panels/RoomsPanel.tsx`
- Modify: `frontend/src/panels/registry.ts`

- [ ] **Step 1: Write the component**

Reuse `useRoomsStore` from `frontend/src/roomsStore.tsx` and the existing per-row actions (star, copy link, leave/delete) from `RoomManagementMenu`. Study that file first — port its row rendering verbatim, minus the menu-dropdown wrapper.

Create `frontend/src/panels/RoomsPanel.tsx`:

```tsx
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import {
	Box,
	IconButton,
	List,
	ListItem,
	ListItemButton,
	ListItemText,
	Typography,
} from "@mui/material";
import { useEffect } from "react";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import { getDockviewApi } from "./DockviewLayout";

export function RoomsPanel() {
	const rooms = useRoomsStore((s) => s.roomsArray);
	const loading = useRoomsStore((s) => s.loading);
	const fetchRooms = useRoomsStore((s) => s.fetchRooms);
	const currentRoomId = useAppStore((s) => s.roomId);
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const leaveRoom = useLeaveRoom({ api: getDockviewApi() });

	useEffect(() => {
		fetchRooms();
	}, [fetchRooms]);

	const switchToRoom = async (roomId: string) => {
		if (roomId === currentRoomId) return;
		await leaveRoom({ skipConfirm: true });
		window.history.pushState({}, "", `/rooms/${roomId}`);
		window.dispatchEvent(new PopStateEvent("popstate"));
	};

	const copyLink = (roomId: string) => {
		navigator.clipboard.writeText(`${window.location.origin}/rooms/${roomId}`);
		showSnackbar("Link copied", "success");
	};

	if (loading) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography>Loading rooms…</Typography>
			</Box>
		);
	}
	if (rooms.length === 0) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography color="text.secondary">
					No rooms available — create one via the Filesystem panel.
				</Typography>
			</Box>
		);
	}
	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Typography variant="overline" sx={{ px: 2, pt: 1 }}>
				Rooms
			</Typography>
			<List dense sx={{ flexGrow: 1, overflow: "auto" }}>
				{rooms.map((r) => (
					<ListItem
						key={r.id}
						secondaryAction={
							<IconButton size="small" onClick={() => copyLink(r.id)}>
								<ContentCopyIcon fontSize="small" />
							</IconButton>
						}
					>
						<ListItemButton
							selected={r.id === currentRoomId}
							onClick={() => switchToRoom(r.id)}
						>
							<ListItemText primary={r.name ?? r.id} secondary={r.id} />
						</ListItemButton>
					</ListItem>
				))}
			</List>
		</Box>
	);
}
```

**Note:** `window.dispatchEvent(new PopStateEvent('popstate'))` is the idiomatic way to trigger react-router-dom's navigation listener after a manual `history.pushState`. Verify the app's existing router (created with `createBrowserRouter` in `App.tsx`) reacts to `popstate`. If not, use `useNavigate()` from react-router-dom inside the component instead:

```tsx
import { useNavigate } from "react-router-dom";
// ...
const navigate = useNavigate();
navigate(`/rooms/${roomId}`);
```

- [ ] **Step 2: Swap stub for real component**

Edit `frontend/src/panels/registry.ts`:

```tsx
import { RoomsPanel } from "./RoomsPanel";
// ...
rooms: {
	kind: "tool",
	icon: MeetingRoom,
	label: "Rooms",
	component: RoomsPanel,
	default: { bar: "left", order: 5 },
},
```

- [ ] **Step 3: Typecheck + lint**

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

- [ ] **Step 4: Commit**

```bash
git add frontend/src/panels/RoomsPanel.tsx frontend/src/panels/registry.ts
git commit -m "feat(panels): rooms panel with in-place switching"
```

---

## Task 10: FilesystemPanel

**Files:**
- Create: `frontend/src/panels/FilesystemPanel.tsx`
- Modify: `frontend/src/panels/registry.ts`

- [ ] **Step 1: Read `pages/filesystemBrowser.tsx` carefully** — it holds all the logic (providers, file listing, LoadFileDialog, navigation). Extract the inner content (the `<Container><Paper>...</Paper></Container>` section) into the new panel, removing the outer `<AppBar>` and `<Container>` chrome.

- [ ] **Step 2: Write the panel**

Create `frontend/src/panels/FilesystemPanel.tsx`. Start by copying `pages/filesystemBrowser.tsx` and deleting the AppBar toolbar block (lines with `<AppBar position="static">...</AppBar>`) plus the `<Container maxWidth="lg">` wrapper. Wrap content in:

```tsx
<Box
	data-testid="filesystem-panel"
	sx={{ display: "flex", flexDirection: "column", height: "100%", overflow: "auto", p: 2 }}
>
	{/* inner selector + breadcrumbs + glob + file list + dialog */}
</Box>
```

Replace `roomId` sourcing: previously from `useParams` — now use:

```tsx
const roomId = useAppStore((s) => s.roomId);
```

Replace the post-load navigation — currently `navigate(\`/rooms/${newRoomId}\`)`. Now go through `useLeaveRoom` first so current plots cascade close:

```tsx
const leaveRoom = useLeaveRoom({ api: getDockviewApi() });
const navigate = useNavigate();

const handleLoadFile = async (params: LoadFileParams) => {
	// ... existing POST to submit task ...
	const { new_room_id } = response.data;
	await leaveRoom({ skipConfirm: true });
	navigate(`/rooms/${new_room_id}`);
};
```

- [ ] **Step 3: Swap stub for real component**

Edit `frontend/src/panels/registry.ts`:

```tsx
import { FilesystemPanel } from "./FilesystemPanel";
// ...
filesystem: {
	kind: "tool",
	icon: Folder,
	label: "Files",
	component: FilesystemPanel,
	default: { bar: "left", order: 6 },
},
```

- [ ] **Step 4: Typecheck + lint**

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

- [ ] **Step 5: Commit**

```bash
git add frontend/src/panels/FilesystemPanel.tsx frontend/src/panels/registry.ts
git commit -m "feat(panels): filesystem panel (unwrapped from page route)"
```

---

## Task 11: ChatPanel + Unread Badge Migration

**Files:**
- Create: `frontend/src/panels/ChatPanel.tsx`
- Modify: `frontend/src/panels/registry.ts`
- Modify: `frontend/src/hooks/socketHandlers/chatHandlers.ts`

- [ ] **Step 1: Read `src/components/ChatWindow.tsx`** — identify the inner content (message list + input + markdown pipeline). Everything outside `<Rnd>` gets kept; the `<Rnd>` and its drag-handle go away.

- [ ] **Step 2: Write the panel**

Create `frontend/src/panels/ChatPanel.tsx`. Start by copying `ChatWindow.tsx` and:
1. Remove the outer `<Rnd>` and its state (size, position, zIndex refs).
2. Remove `dragHandleClassName` and the drag-handle `<Box>`.
3. Replace the `<Paper>` wrapper with:

```tsx
<Box
	data-testid="chat-panel"
	sx={{ display: "flex", flexDirection: "column", height: "100%", overflow: "hidden" }}
>
	{/* existing message list + input */}
</Box>
```

4. When the panel mounts, call `resetChatUnread()` (it's now the active panel in its bar).

```tsx
useEffect(() => {
	useAppStore.getState().resetChatUnread();
}, []);
```

- [ ] **Step 3: Swap stub**

Edit `frontend/src/panels/registry.ts`:

```tsx
import { ChatPanel } from "./ChatPanel";
// ...
chat: {
	kind: "tool",
	icon: Chat,
	label: "Chat",
	component: ChatPanel,
	default: { bar: "right", order: 0 },
},
```

- [ ] **Step 4: Fix `chatHandlers.ts` unread gating**

Edit `frontend/src/hooks/socketHandlers/chatHandlers.ts`. Replace:

```ts
const { chatOpen, incrementChatUnread } = useAppStore.getState();
if (!chatOpen) {
	incrementChatUnread();
}
```

with:

```ts
const state = useAppStore.getState();
const chatIsActive =
	state.activeLeft === "chat" ||
	state.activeRight === "chat" ||
	state.activeBottom === "chat";
if (!chatIsActive) {
	state.incrementChatUnread();
}
```

- [ ] **Step 5: Remove `chatOpen` / `setChatOpen` from uiSlice**

Edit `frontend/src/stores/slices/uiSlice.ts`:
- Delete `chatOpen: boolean;`, `setChatOpen: (open: boolean) => void;` from the interface.
- Delete `chatOpen: false,` and `setChatOpen: (open) => set({ chatOpen: open }),` from the implementation.
- Keep `chatUnreadCount`, `incrementChatUnread`, `resetChatUnread`, `typingUsers`, `addTypingUser`, `removeTypingUser`.

- [ ] **Step 6: Display unread badge on activity bar icon**

Edit `frontend/src/panels/ActivityBar.tsx` — the chat icon needs a badge. Add inside the `icons.map` loop, before `<IconButton>`:

```tsx
import Badge from "@mui/material/Badge";
// ...
const unread = useAppStore((s) => s.chatUnreadCount);
// inside map:
const iconInner = (
	<Icon />
);
const iconMaybeWithBadge =
	id === "chat" && !isActive && unread > 0 ? (
		<Badge badgeContent={unread} color="error" max={99}>
			{iconInner}
		</Badge>
	) : (
		iconInner
	);
// replace <Icon /> with {iconMaybeWithBadge}
```

- [ ] **Step 7: Typecheck + lint**

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

Expected: will still fail if `landingPage.tsx` or another file still reads `chatOpen`. Let it fail for now — Task 12 fixes the consumer. If fixing here is trivial (one reference), fix in this task.

- [ ] **Step 8: Commit**

```bash
git add frontend/src/panels/ChatPanel.tsx frontend/src/panels/registry.ts frontend/src/panels/ActivityBar.tsx frontend/src/hooks/socketHandlers/chatHandlers.ts frontend/src/stores/slices/uiSlice.ts
git commit -m "feat(panels): chat panel with bar badge"
```

---

## Task 12: Replace landingPage Layout + AppBar Cleanup

**Files:**
- Modify: `frontend/src/pages/landingPage.tsx`

- [ ] **Step 1: Read current `landingPage.tsx` fully** — note all imports that need to stay (dialogs, hooks, socket manager).

- [ ] **Step 2: Swap stubs for real tool panels**

Edit `frontend/src/panels/registry.ts` imports and entries to use:

```tsx
import { GeometryPanel } from "../components/geometry/GeometryPanel";
import { SecondaryPanel } from "../components/SecondaryPanel";
import { SelectionsPanel } from "../components/SelectionsPanel";
// ...
selections: { kind: "tool", icon: FilterCenterFocus, label: "Selections", component: SelectionsPanel, default: { bar: "left", order: 0 } },
modifiers:  { kind: "tool", icon: Build, label: "Modifiers", component: () => <SecondaryPanel panelTitle="modifiers" />, default: { bar: "left", order: 1 } },
analysis:   { kind: "tool", icon: Analytics, label: "Analysis", component: () => <SecondaryPanel panelTitle="analysis" />, default: { bar: "left", order: 2 } },
geometries: { kind: "tool", icon: Category, label: "Geometries", component: GeometryPanel, default: { bar: "left", order: 3 } },
```

Note: `SecondaryPanel` and `SelectionsPanel` currently read `selectedCategory`/`selectedExtensions` from `useFormStore`. Move their `selectedExtensions` state to the existing `formStore` **only if still used** — otherwise (once all consumers migrate to activity bar), delete `formStore.ts` in Task 14. For now leave it; if the panels read `selectedCategory` to decide rendering, keep the getter returning the bar-derived active panel id.

Quickest path: drop `selectedCategory` reads entirely; these components are now only rendered when active, so they don't need to check. Inside `SecondaryPanel`/`SelectionsPanel`, delete any code of the form:

```tsx
const selectedCategory = useFormStore(s => s.selectedCategory);
if (selectedCategory !== "modifiers") return null;
```

Keep `selectedExtensions`/`setSelectedExtension` usage (per-extension-name form state).

- [ ] **Step 3: Replace `pages/landingPage.tsx` JSX**

Target structure:

```tsx
import {
	ActivityBar,
	BottomZone,
	DockviewLayout,
	SidebarZone,
} from "../panels";
// ... keep all other imports (socket manager, dialogs, app store, FrameProgressBar, etc.)

// Keep AppBar contents, but REMOVE:
//   - AddPlotButton import + usage
//   - chat IconButton + badge (moved to ActivityBar)
//   - RoomManagementMenu import + usage
//   - chatOpen / setChatOpen state reads
//   - ChatWindow lazy import and render
//   - SideBar, MyScene, WindowManager, DropOverlay, drag-boundary-container

// Add a "Reset layout" item to the "..." menu (existing):
//   onClick: () => {
//     useAppStore.getState().resetLayout();
//     const api = getDockviewApi();
//     if (api) {
//       api.panels.forEach(p => p.api.close());
//       api.addPanel({ id: "viewer", component: "viewer", title: "3D Viewer" });
//     }
//   }

return (
	<Box sx={{ display: "flex", flexDirection: "column", height: "100vh", width: "100vw", overflow: "hidden" }}>
		<CssBaseline />
		<AppBar position="static" sx={{ height: APPBAR_HEIGHT }}>
			<Toolbar>{/* kept toolbar items (theme, drawing, editing, screenshot, upload, profile, ...) */}</Toolbar>
		</AppBar>

		<Box sx={{ display: "flex", flexGrow: 1, minHeight: 0 }}>
			<ActivityBar position="left" />
			<SidebarZone position="left" />
			<DockviewLayout />
			<SidebarZone position="right" />
			<ActivityBar position="right" />
		</Box>

		<BottomZone />
		<ActivityBar position="bottom" />
		<FrameProgressBar />

		{/* Dialogs unchanged */}
		<ConnectionDialog />
		<LoginDialog />
		<RegisterDialog />
		<AdminPanel />
		<UserProfileDialog />
		<SiMGenTutorialDialog />
		<Snackbar />
		<ProgressNotifications />
	</Box>
);
```

- [ ] **Step 4: Typecheck + lint — must pass**

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

Expected: 0 errors. If `chatOpen` or `setChatOpen` references remain, remove them now.

- [ ] **Step 5: Playwright validation (REQUIRED)**

With `uv run zndraw` + `bun run dev` running:

```bash
# Default startup — left bar icons visible, viewer filling center, no sidebar open
playwright-cli -s=dockview open http://localhost:5173/rooms/dockview-test
playwright-cli -s=dockview snapshot --filename=task12-startup.yaml
playwright-cli -s=dockview screenshot --filename=task12-startup.png
playwright-cli -s=dockview console  # must be empty

# Open each left-bar panel
for panel in selections modifiers analysis geometries plots-browser rooms filesystem; do
  playwright-cli -s=dockview click "[data-testid=activity-icon-${panel}]"
  playwright-cli -s=dockview snapshot --filename=task12-panel-${panel}.yaml
  playwright-cli -s=dockview click "[data-testid=activity-icon-${panel}]"  # close
done

# Open chat (right bar)
playwright-cli -s=dockview click "[data-testid=activity-icon-chat]"
playwright-cli -s=dockview snapshot --filename=task12-chat.yaml

# Drag selections icon to right bar
playwright-cli -s=dockview drag "[data-testid=activity-icon-selections]" "[data-testid=activity-bar-right]"
playwright-cli -s=dockview snapshot --filename=task12-dragged.yaml

# Reset layout via "..." menu
# (click the ... IconButton, then "Reset layout")
playwright-cli -s=dockview click "[aria-label='more']"
playwright-cli -s=dockview click "text=Reset layout"
playwright-cli -s=dockview snapshot --filename=task12-reset.yaml

playwright-cli -s=dockview console
```

Verify each snapshot shows the expected structure (left bar has icons, opening an icon reveals the panel content in a `sidebar-zone-left` container, chat badge works, drag moves icon to new bar, reset restores defaults). No console errors.

- [ ] **Step 6: Commit**

```bash
git add frontend/src/pages/landingPage.tsx frontend/src/panels/registry.ts frontend/src/components/SelectionsPanel.tsx frontend/src/components/SecondaryPanel.tsx
git commit -m "feat(ui): replace landing page layout with dockview + activity bars"
```

---

## Task 13: Redirect /filesystem Route

**Files:**
- Modify: `frontend/src/App.tsx`

- [ ] **Step 1: Replace route**

Edit `frontend/src/App.tsx` — replace the `/rooms/:roomId/files` entry with a redirect component:

```tsx
import { Navigate, useParams, useSearchParams } from "react-router-dom";

function FilesystemRedirect() {
	const { roomId } = useParams<{ roomId: string }>();
	return <Navigate to={`/rooms/${roomId}?panel=filesystem`} replace />;
}

// in router:
{
	path: "/rooms/:roomId/files",
	element: <FilesystemRedirect />,
},
```

- [ ] **Step 2: Auto-open filesystem panel on `?panel=filesystem`**

Edit `frontend/src/pages/landingPage.tsx` — on mount, read `useSearchParams()` and, if `panel` is a known tool, set that panel active on its default bar:

```tsx
import { useSearchParams } from "react-router-dom";
import { PANELS, type PanelId } from "../panels";
// ...
const [params] = useSearchParams();
useEffect(() => {
	const panel = params.get("panel") as PanelId | null;
	if (panel && panel in PANELS && PANELS[panel].kind === "tool") {
		const def = PANELS[panel];
		if (def.kind === "tool" && def.default.bar !== "editor") {
			useAppStore.getState().toggleActive(def.default.bar, panel);
		}
	}
}, [params]);
```

- [ ] **Step 3: Validate**

```bash
playwright-cli -s=dockview goto http://localhost:5173/rooms/dockview-test/files
# Expected: redirect to /rooms/dockview-test?panel=filesystem and filesystem panel visible
playwright-cli -s=dockview snapshot --filename=task13-redirect.yaml
playwright-cli -s=dockview console
```

- [ ] **Step 4: Commit**

```bash
git add frontend/src/App.tsx frontend/src/pages/landingPage.tsx
git commit -m "feat(routes): redirect /files route to panel-based filesystem"
```

---

## Task 14: Delete Removed Files and Symbols

**Files:**
- Delete: `frontend/src/components/SideBar.tsx`
- Delete: `frontend/src/components/PrimaryDrawer.tsx`
- Delete: `frontend/src/components/WindowManager.tsx`
- Delete: `frontend/src/components/FigureWindow.tsx`
- Delete: `frontend/src/components/AddPlotButton.tsx`
- Delete: `frontend/src/components/ChatWindow.tsx`
- Delete: `frontend/src/stores/windowManagerStore.ts`
- Delete: `frontend/src/formStore.ts`
- Delete: `frontend/src/pages/filesystemBrowser.tsx`
- Modify: `frontend/package.json` (remove `react-rnd`)

- [ ] **Step 1: Run ripgrep sweep to find remaining imports**

```bash
cd /Users/fzills/tools/zndraw-fastapi/frontend
bunx rg -l --fixed-strings \
  -e 'SideBar' -e 'PrimaryDrawer' -e 'WindowManager' -e 'FigureWindow' \
  -e 'AddPlotButton' -e 'ChatWindow' -e 'windowManagerStore' -e 'formStore' \
  -e 'filesystemBrowser' -e 'react-rnd' \
  src/
```

Expected: only the files about to be deleted should reference these. Any other hits = fix before deletion.

- [ ] **Step 2: Delete files**

```bash
cd /Users/fzills/tools/zndraw-fastapi
git rm frontend/src/components/SideBar.tsx
git rm frontend/src/components/PrimaryDrawer.tsx
git rm frontend/src/components/WindowManager.tsx
git rm frontend/src/components/FigureWindow.tsx
git rm frontend/src/components/AddPlotButton.tsx
git rm frontend/src/components/ChatWindow.tsx
git rm frontend/src/stores/windowManagerStore.ts
git rm frontend/src/formStore.ts
git rm frontend/src/pages/filesystemBrowser.tsx
```

- [ ] **Step 3: Check that `SecondaryPanel` / `SelectionsPanel` don't need formStore**

```bash
cd frontend
bunx rg -n 'useFormStore|formStore' src/
```

If `selectedExtensions` is still used by these panels, inline that state into a small local slice or into `SceneSlice` / a new extensionFormSlice. If not used, this returns no matches.

If still needed, create `frontend/src/stores/slices/extensionFormSlice.ts`:

```ts
import type { StateCreator } from "zustand";
import type { AppState } from "../../store";

export interface ExtensionFormSlice {
	selectedExtensions: Record<string, string | null>;
	setSelectedExtension: (category: string, extension: string | null) => void;
}

export const createExtensionFormSlice: StateCreator<
	AppState,
	[],
	[],
	ExtensionFormSlice
> = (set) => ({
	selectedExtensions: {},
	setSelectedExtension: (category, extension) =>
		set((state) => ({
			selectedExtensions: { ...state.selectedExtensions, [category]: extension },
		})),
});
```

Wire into `store.tsx` and update the two panels' imports. This is cleaner than keeping the orphaned `formStore.ts`.

- [ ] **Step 4: Remove react-rnd**

```bash
cd frontend
bun remove react-rnd
```

- [ ] **Step 5: Verification sweep**

```bash
cd /Users/fzills/tools/zndraw-fastapi/frontend
bunx rg -i 'react-rnd|WindowManager|PrimaryDrawer|selectedCategory|openWindow|FigureWindow|AddPlotButton|SideBar\.tsx|windowManagerStore|formStore|filesystemBrowser\.tsx' src/
```

Expected: **zero matches**. If any remain, they are artifacts — delete.

- [ ] **Step 6: Typecheck + lint**

```bash
cd frontend && bunx tsc --noEmit && bun run lint
```

Expected: 0 errors.

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "chore(frontend): remove legacy sidebar/window-manager/chat-window and react-rnd"
```

---

## Task 15: E2E Playwright Test for Core Flows

**Files:**
- Create: `frontend/e2e/dockview-layout.spec.ts`

- [ ] **Step 1: Write the test**

Create `frontend/e2e/dockview-layout.spec.ts`:

```ts
import { expect, test } from "@playwright/test";

test.describe("dockview layout", () => {
	test("default startup shows activity bars and viewer", async ({ page }) => {
		await page.goto("/rooms/dockview-test");
		await expect(page.getByTestId("activity-bar-left")).toBeVisible();
		await expect(page.getByTestId("activity-bar-right")).toBeVisible();
		await expect(page.getByTestId("activity-bar-bottom")).toBeVisible();
		await expect(page.getByTestId("viewer-view")).toBeVisible();
		await expect(page.getByTestId("sidebar-zone-left")).toBeHidden();
	});

	test("clicking selections icon opens the left sidebar zone", async ({ page }) => {
		await page.goto("/rooms/dockview-test");
		await page.getByTestId("activity-icon-selections").click();
		await expect(page.getByTestId("sidebar-zone-left")).toBeVisible();
		await page.getByTestId("activity-icon-selections").click();
		await expect(page.getByTestId("sidebar-zone-left")).toBeHidden();
	});

	test("only one panel per bar", async ({ page }) => {
		await page.goto("/rooms/dockview-test");
		await page.getByTestId("activity-icon-selections").click();
		await page.getByTestId("activity-icon-modifiers").click();
		// sidebar-zone-left still exists but shows only the modifiers panel
		await expect(page.getByTestId("sidebar-zone-left")).toBeVisible();
	});

	test("opening plots browser and clicking a plot opens a tab", async ({ page }) => {
		await page.goto("/rooms/dockview-test");
		await page.getByTestId("activity-icon-plots-browser").click();
		// Assumes the test room has at least one figure seeded
		const firstRow = page.getByTestId("sidebar-zone-left").locator("li button").first();
		await expect(firstRow).toBeVisible();
		await firstRow.click();
		// Plot panel added; dockview renders it as [data-testid^=plot-view-]
		await expect(page.locator('[data-testid^="plot-view-"]')).toBeVisible();
	});

	test("closing viewer shows welcome state", async ({ page }) => {
		await page.goto("/rooms/dockview-test");
		await expect(page.getByTestId("viewer-view")).toBeVisible();
		// Confirm dialog is handled auto-dismiss
		page.on("dialog", (d) => d.accept());
		// Simulate clicking the viewer tab close via dockview UI
		const closeBtn = page.locator(".dv-tab .dv-default-tab-action").first();
		await closeBtn.click();
		await expect(page.getByTestId("dockview-welcome")).toBeVisible();
	});
});
```

- [ ] **Step 2: Seed test data (preflight)**

```bash
export ZNDRAW_ROOM=dockview-test
uv run zndraw-cli rooms create --id dockview-test 2>/dev/null || true
uv run zndraw-cli frames extend test_water.xyz
uv run python -c "
from zndraw import ZnDraw
import plotly.graph_objects as go
vis = ZnDraw(room='dockview-test')
vis.figures['test_plot'] = go.Figure(go.Scatter(x=[0,1,2], y=[0,1,4]))
"
```

- [ ] **Step 3: Build + run test**

```bash
cd frontend && bun run build
ZNDRAW_URL=http://localhost:8000 bunx playwright test e2e/dockview-layout.spec.ts
```

Expected: all 5 tests pass.

- [ ] **Step 4: Fix failures before proceeding.** Do NOT `git commit --no-verify` and do NOT disable the tests.

- [ ] **Step 5: Commit**

```bash
git add frontend/e2e/dockview-layout.spec.ts
git commit -m "test(e2e): dockview layout, panel toggling, plot opening, viewer close cascade"
```

---

## Task 16: Manual Smoke Test + Final Verification

**Files:** (none — validation only)

- [ ] **Step 1: Final rg sweep**

```bash
cd /Users/fzills/tools/zndraw-fastapi/frontend
bunx rg -i 'react-rnd|WindowManager|PrimaryDrawer|selectedCategory|openWindow|FigureWindow|AddPlotButton|SideBar\.tsx|windowManagerStore|formStore|filesystemBrowser\.tsx' src/ e2e/
```

Expected: **zero matches**. If any remain, go back and remove.

- [ ] **Step 2: Lint + typecheck**

```bash
cd frontend && bun run lint && bunx tsc --noEmit && bun run build
```

Expected: all three exit 0.

- [ ] **Step 3: Python tests unchanged**

```bash
cd /Users/fzills/tools/zndraw-fastapi
uv run pytest -x -q
```

Expected: all pass (this is a frontend-only refactor — backend shouldn't regress).

- [ ] **Step 4: Full E2E suite**

```bash
cd frontend
ZNDRAW_URL=http://localhost:8000 bunx playwright test
```

Expected: all pass. `chat-features.spec.ts`, `ui-panels-chat.spec.ts`, and similar may need selector updates — fix any breakage, do not skip.

- [ ] **Step 5: Manual smoke with playwright-cli**

```bash
playwright-cli -s=dockview open http://localhost:5173/rooms/dockview-test
playwright-cli -s=dockview snapshot --filename=smoke-startup.yaml

# Tools flow
playwright-cli -s=dockview click "[data-testid=activity-icon-selections]"
playwright-cli -s=dockview click "[data-testid=activity-icon-modifiers]"
playwright-cli -s=dockview click "[data-testid=activity-icon-geometries]"

# Plots flow — click row + drag row
playwright-cli -s=dockview click "[data-testid=activity-icon-plots-browser]"
playwright-cli -s=dockview click "text=test_plot"

# Auto-open: set a figure via CLI and verify a tab appears
uv run python -c "
from zndraw import ZnDraw
import plotly.graph_objects as go
vis = ZnDraw(room='dockview-test')
vis.figures['auto_opened'] = go.Figure(go.Scatter(x=[0,1], y=[1,2]))
"
playwright-cli -s=dockview snapshot --filename=smoke-auto-plot.yaml
# Expected: sidebar-zone-left shows 'auto_opened' and a plot-view-auto_opened panel is in the editor area

# Auto-close: delete the figure
uv run zndraw-cli figures delete auto_opened
playwright-cli -s=dockview snapshot --filename=smoke-auto-delete.yaml
# Expected: plot-view-auto_opened panel is gone

# Drag a plot tab to split (shift+drag the tab) — manual gesture, just confirm the editor area still works
playwright-cli -s=dockview screenshot --filename=smoke-final.png

# Room switching
playwright-cli -s=dockview click "[data-testid=activity-icon-rooms]"
uv run zndraw-cli rooms create --id dockview-smoke-2 2>/dev/null || true
playwright-cli -s=dockview click "text=dockview-smoke-2"
# Accept confirm if any plot tabs were open
playwright-cli -s=dockview dialog-accept
playwright-cli -s=dockview snapshot --filename=smoke-switched-room.yaml
playwright-cli -s=dockview eval "window.location.pathname"
# Expected output: "/rooms/dockview-smoke-2"

# Reset layout via ... menu
playwright-cli -s=dockview click "[aria-label='more']"
playwright-cli -s=dockview click "text=Reset layout"
playwright-cli -s=dockview snapshot --filename=smoke-reset.yaml

# Console must be empty
playwright-cli -s=dockview console
playwright-cli -s=dockview close
```

Save `/tmp/dockview-validation/*.yaml` + `*.png` for review.

- [ ] **Step 6: zndraw CLI end-to-end**

```bash
# Verify backend still works with the new frontend
uv run zndraw-cli rooms list
uv run zndraw-cli figures list
uv run zndraw-cli step set 5   # viewer should jump to frame 5
```

Watch the browser — each command should reflect live.

- [ ] **Step 7: No commit.** If all checks passed, the work is ready for PR. Use superpowers:requesting-code-review next.

---

## Self-Review

Spec coverage audit (each line = one spec requirement, one task reference):

- VS Code-like three bars + editor → Tasks 3, 4, 5, 12
- Drag icons between bars → Task 3
- One panel per bar → Task 2 (toggleActive), Task 4 (SidebarZone)
- Editor area splittable / floatable / popout → Task 5, 7 (default dockview behavior)
- Plots Browser click + drag → Task 8
- Rooms panel, in-place switching → Task 9
- Filesystem panel, `/filesystem` redirect → Task 10, 13
- Viewer close = leave room, cascade plots → Task 6
- Auto-open on figure_invalidate op=set → Task 7
- Auto-close on figure_invalidate op=delete → Task 7
- Chat panel (right bar, no Rnd, bar-icon badge) → Task 11
- `chatOpen` / `setChatOpen` removed, readers migrated → Task 11
- AppBar cleanup (AddPlot, chat icon, RoomManagementMenu removed; Reset layout added) → Task 12
- Single PANELS registry → Tasks 1, 8, 9, 10, 11, 12
- Full replacement / no artifacts → Tasks 14, 16 (rg sweep)
- Non-goals (persistence, multi-viewport, etc.) → not implemented, correct per spec
- Backend unchanged → Task 16 step 3

Placeholder scan: plan contains no "TBD", "later", "as appropriate". Marker-track porting in Task 7 step 1 is flagged as "port verbatim from FigureWindow.tsx lines 456–985" with source file + line numbers — the engineer has exact coordinates.

Type consistency: `PanelId` union (Task 1), `BarPosition` type (Task 1), `ActivityBarSlice` API (Task 2), `openPlotTab` / `closePlotTab` / `plotPanelId` (Task 7), `getDockviewApi` (Task 6 → used in 7, 8, 9, 10, 11, 12, 14). All match.

---

Plan complete. Two execution options:

**1. Subagent-Driven (recommended)** — Dispatch a fresh subagent per task, review between tasks, fast iteration. Uses superpowers:subagent-driven-development.

**2. Inline Execution** — Execute tasks in this session using superpowers:executing-plans, batch execution with checkpoints.

Which approach?
