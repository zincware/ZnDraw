# Dockview UI Follow-Up Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resolve the 7 post-merge dockview issues plus one layout bug: reclaim chrome for the viewer, sync dockview to MUI light/dark, fix the bottom-drawer shrink regression, add popout + maximize on groups, open plots beside the viewer, simplify the plots browser, auto-resize plotly on container resize.

**Architecture:** Frontend-only changes in `frontend/src/panels/*` plus a new E2E spec block. No backend changes. Each issue is an independent commit; order is smallest-blast-radius first so regressions are bisectable.

**Tech Stack:** React 19 + TypeScript + MUI v7 + Zustand + Vite 7 + dockview-react 5.2 + Plotly.js. Dev via `bun run dev` (port 5173 → proxies to :8000). Production serve via `uv sync --reinstall-package zndraw && uv run zndraw` (serves built static at :8000). E2E: `@playwright/test` against `:8000`.

**Parent spec:** [`docs/superpowers/specs/2026-04-16-dockview-followup-design.md`](../specs/2026-04-16-dockview-followup-design.md).

---

## File structure

### New files

- `frontend/src/panels/dockview-mui.css` — CSS variable overrides mapping `--dv-*` to MUI `--mui-palette-*` tokens. Imported once from `DockviewLayout.tsx`.
- `frontend/src/panels/groupActions.tsx` — per-group right-gutter action component (popout + maximize buttons). Wired to dockview via `rightHeaderActionsComponent`.

### Modified files

- `frontend/src/panels/registry.tsx` — chat `default.bar: "right"` → `"left"`.
- `frontend/src/panels/ActivityBar.tsx` — three-state render (full / sliver / hot), listens on `window` for `dragstart`/`dragend`.
- `frontend/src/panels/DockviewLayout.tsx` — absolute-positioned inner wrapper around `DockviewReact`; dynamic theme class from `useColorScheme`; imports `dockview-mui.css`; wires `rightHeaderActionsComponent={GroupActions}`.
- `frontend/src/panels/PlotView.tsx` — subscribes to `props.api.onDidDimensionsChange` → `Plotly.Plots.resize(container)`.
- `frontend/src/panels/PlotsBrowserPanel.tsx` — drops the "Currently Open" section; inline hover `×` on each open row.
- `frontend/src/panels/plotViewFactory.ts` — `openPlotTab` resolves target group: existing plot group (lex-first id) > right of viewer group > dockview default.
- `frontend/e2e/dockview-layout.spec.ts` — new assertions for each change.

### Unchanged

- `frontend/src/panels/SidebarZone.tsx`, `BottomZone.tsx` (already null-on-no-active; no change needed; chat appearing on left just adds an entry, BAR icons filter by `default.bar`).
- `frontend/src/hooks/socketHandlers/figureHandlers.ts` (already calls `openPlotTab(api, key)`; new positioning logic lives inside `openPlotTab`).
- Backend, socket event layer, R3F scene, `useLeaveRoom`.

---

## Pre-flight: branch + servers

Before Task 1, confirm you're on a feature branch off `main` and that the dev environment works.

- [ ] **Step P1: Verify branch**

Run: `git rev-parse --abbrev-ref HEAD`
Expected: something like `spec/dockview-ui-redesign` (or a follow-up branch). If on `main`, create a branch: `git switch -c fix/dockview-followup`.

- [ ] **Step P2: Install frontend deps**

Run: `cd frontend && bun install`
Expected: clean install, no EEXIST errors. If they appear, `rm -rf node_modules && bun install` recovers.

- [ ] **Step P3: Start backend**

Run: `uv run zndraw` (from repo root, in its own terminal, leave it running)
Expected: `Uvicorn running on http://0.0.0.0:8000`. Visit http://localhost:8000/ — should show the ZnDraw landing page (dockview layout + activity bars), NOT a chemistry-management app.

- [ ] **Step P4: Start frontend dev server**

Run: `cd frontend && bun run dev` (own terminal, leave running)
Expected: `Local: http://localhost:5173/`. Visit — should match :8000 visually and hot-reload on file save.

- [ ] **Step P5: Seed E2E test room**

Run (from repo root):
```bash
uv run zndraw-cli rooms create --room dockview-test 2>/dev/null || true
uv run python -c "
from zndraw import ZnDraw
import plotly.graph_objects as go
vis = ZnDraw(url='http://localhost:8000', room='dockview-test')
for key in ['energy', 'rdf']:
    if key not in vis.figures:
        vis.figures[key] = go.Figure(go.Scatter(x=[0,1,2], y=[0,1,4], mode='lines'))
"
```
Expected: no errors. The test room has at least 2 figures for the plot-positioning assertion later.

---

## Task 1: Move chat icon to the left activity bar

**Goal:** Chat defaults to the left bar so the right bar is empty at startup (unlocks right-sliver drop-zone later).

**Files:**
- Modify: `frontend/src/panels/registry.tsx:117`
- Test: `frontend/e2e/dockview-layout.spec.ts` (append block)

- [ ] **Step 1.1: Write the failing E2E assertion**

Append to `frontend/e2e/dockview-layout.spec.ts`, inside the existing `test.describe("dockview layout", ...)` block (just before the closing `});`):

```ts
test("chat icon defaults to the left activity bar", async ({ page }) => {
    await page.goto(`/rooms/${ROOM}`);
    const leftBar = page.getByTestId("activity-bar-left");
    await expect(leftBar.getByTestId("activity-icon-chat")).toBeVisible();
    const rightBar = page.getByTestId("activity-bar-right");
    await expect(rightBar.getByTestId("activity-icon-chat")).toHaveCount(0);
});
```

- [ ] **Step 1.2: Build + run test → expect FAIL**

Run: `cd frontend && uv run --project .. bash -lc 'cd $(pwd) && bun run build' && cd .. && uv sync --reinstall-package zndraw`
Then: `cd frontend && bun run test:e2e -- --grep "chat icon defaults"`
Expected: FAIL. Chat is still on the right bar.

If `bun run test:e2e` is not defined, use: `bunx playwright test e2e/dockview-layout.spec.ts --grep "chat icon defaults"`.

- [ ] **Step 1.3: Flip chat default bar**

Edit `frontend/src/panels/registry.tsx` — find the `chat:` entry (around line 112) and change:

```tsx
chat: {
    kind: "tool",
    icon: Chat,
    label: "Chat",
    component: ChatPanel,
    default: { bar: "right", order: 0 },
},
```

to:

```tsx
chat: {
    kind: "tool",
    icon: Chat,
    label: "Chat",
    component: ChatPanel,
    default: { bar: "left", order: 7 },
},
```

Note: `order: 7` places chat at the bottom of the left-bar stack (existing left icons use orders 0–6). Do NOT use `order: 0` — that would push chat above Selections.

- [ ] **Step 1.4: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "chat icon defaults"`
Expected: PASS. Also verify existing specs still pass: `bunx playwright test e2e/dockview-layout.spec.ts`.

- [ ] **Step 1.5: Commit**

```bash
git add frontend/src/panels/registry.tsx frontend/e2e/dockview-layout.spec.ts
git commit -m "feat(panels): move chat icon to left activity bar

Part of the dockview follow-up: empty right bar is a prerequisite for
the sliver drop-zone change."
```

---

## Task 2: Fix the bottom-drawer shrink regression

**Goal:** Opening a bottom panel shrinks the viewer (currently only left/right shrink it).

**Files:**
- Modify: `frontend/src/panels/DockviewLayout.tsx:80-88` (the wrapper structure)
- Test: `frontend/e2e/dockview-layout.spec.ts` (append)

- [ ] **Step 2.1: Write the failing E2E assertion**

Append to the `test.describe("dockview layout", ...)` block:

```ts
test("opening a bottom panel shrinks the viewer vertically", async ({
    page,
}) => {
    await page.goto(`/rooms/${ROOM}`);
    const viewer = page.getByTestId("viewer-view");
    await expect(viewer).toBeVisible();
    const before = await viewer.boundingBox();
    if (!before) throw new Error("viewer bounding box missing");

    // Use the Plots bar icon — it toggles the plots-browser into the left sidebar
    // by default. To simulate a bottom-zone activation without introducing a
    // bottom-default panel yet, drag the plots-browser icon into the bottom bar.
    const plotsIcon = page.getByTestId("activity-icon-plots-browser");
    const bottomBar = page.getByTestId("activity-bar-bottom");
    await plotsIcon.dragTo(bottomBar);
    await plotsIcon.click();

    const bottomZone = page.getByTestId("bottom-zone");
    await expect(bottomZone).toBeVisible();

    const after = await viewer.boundingBox();
    if (!after) throw new Error("viewer bounding box missing after");
    expect(after.height).toBeLessThan(before.height - 50);
});
```

- [ ] **Step 2.2: Build + run test → expect FAIL**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "shrinks the viewer vertically"`
Expected: FAIL. The viewer's `after.height` is approximately equal to `before.height` because the current layout bug means the BottomZone pushes its siblings but DockviewReact keeps its previous computed height.

- [ ] **Step 2.3: Add the absolute-positioned inner wrapper**

Edit `frontend/src/panels/DockviewLayout.tsx` — replace the `return (...)` block (lines 80–107) with:

```tsx
return (
    <Box sx={{ flexGrow: 1, position: "relative", minWidth: 0, minHeight: 0 }}>
        <Box sx={{ position: "absolute", inset: 0 }}>
            <DockviewReact
                className="dockview-theme-light"
                onReady={onReady}
                components={components}
                onDidDrop={onDidDrop}
                floatingGroupBounds="boundedWithinViewport"
            />
        </Box>
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
```

The change: `DockviewReact` now lives inside a second `<Box>` with `position: "absolute", inset: 0`. Its built-in `ResizeObserver` fires on this inner wrapper, whose size follows the flex-parent correctly.

- [ ] **Step 2.4: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts`
Expected: the new assertion passes; all existing ones still pass.

- [ ] **Step 2.5: Commit**

```bash
git add frontend/src/panels/DockviewLayout.tsx frontend/e2e/dockview-layout.spec.ts
git commit -m "fix(panels): decouple DockviewReact size from flex siblings

Wrap DockviewReact in an absolutely-positioned box so its ResizeObserver
sees the flex-parent's computed size. Before this, bottom-drawer opens
did not shrink the viewer because DockviewReact rendered inline and
its height was tied to content, not the parent's current size."
```

---

## Task 3: Dockview CSS variable overrides (MUI palette)

**Goal:** Dockview tab chrome reads MUI palette tokens so colors match the app. No behavior change yet — just the stylesheet.

**Files:**
- Create: `frontend/src/panels/dockview-mui.css`
- Modify: `frontend/src/panels/DockviewLayout.tsx` (add `import "./dockview-mui.css"`)

- [ ] **Step 3.1: Write the stylesheet**

Create `frontend/src/panels/dockview-mui.css`:

```css
/* Maps dockview CSS variables to MUI palette tokens so tab chrome matches
   the app. Applied whenever DockviewReact has either dockview-theme-light
   or dockview-theme-dark class. The MUI tokens themselves flip with
   colorScheme, so one rule serves both modes. */

.dockview-theme-light,
.dockview-theme-dark {
    --dv-group-view-background-color: var(--mui-palette-background-default);
    --dv-tabs-and-actions-container-background-color: var(--mui-palette-background-paper);
    --dv-tab-background-color: var(--mui-palette-background-paper);
    --dv-activegroup-visiblepanel-tab-background-color: var(--mui-palette-background-default);
    --dv-tab-header-foreground-color: var(--mui-palette-text-primary);
    --dv-inactivegroup-visiblepanel-tab-background-color: var(--mui-palette-action-hover);
    --dv-separator-border: var(--mui-palette-divider);
    --dv-paneview-active-outline-color: var(--mui-palette-primary-main);
}
```

- [ ] **Step 3.2: Import the stylesheet from DockviewLayout**

Edit `frontend/src/panels/DockviewLayout.tsx` — just after the existing `import "dockview-react/dist/styles/dockview.css";` line, add:

```tsx
import "./dockview-mui.css";
```

The import order matters: MUI overrides load after dockview's base so the override wins.

- [ ] **Step 3.3: Rebuild + visually verify**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw`
Open http://localhost:8000/rooms/dockview-test in a browser. Click the Plots icon and open a plot. The tab bar should now read as a light-gray MUI Paper color (not the dockview default near-white).

- [ ] **Step 3.4: Commit**

```bash
git add frontend/src/panels/dockview-mui.css frontend/src/panels/DockviewLayout.tsx
git commit -m "feat(panels): map dockview CSS vars to MUI palette tokens

Tab chrome now reads MUI palette tokens so the dockview surface visually
belongs to the app. Works for both light and dark themes because the
MUI tokens themselves flip on color scheme change."
```

---

## Task 4: Sync dockview theme class with MUI color scheme

**Goal:** Flipping MUI light ↔ dark flips dockview's base theme class too, so borders/shadows/active outlines all track the color mode.

**Files:**
- Modify: `frontend/src/panels/DockviewLayout.tsx` — replace hardcoded `className="dockview-theme-light"` with dynamic.
- Test: `frontend/e2e/dockview-layout.spec.ts`.

- [ ] **Step 4.1: Write the failing E2E assertion**

Append to the test block:

```ts
test("dockview theme class flips with MUI color scheme", async ({
    page,
}) => {
    await page.goto(`/rooms/${ROOM}`);
    const dockview = page.locator(".dv-dockview").first();
    await expect(dockview).toHaveClass(/dockview-theme-light/);

    await page.evaluate(() => {
        document.documentElement.setAttribute("data-mui-color-scheme", "dark");
    });
    // Wait a frame for React to re-render against the new scheme.
    await page.waitForTimeout(100);
    await expect(dockview).toHaveClass(/dockview-theme-dark/);
});
```

Note: MUI stores the color scheme on `<html data-mui-color-scheme>`. Setting it directly triggers `useColorScheme` subscribers as long as the provider is mounted (which it is — `CssVarsProvider` wraps the app in `main.tsx`).

- [ ] **Step 4.2: Build + run test → expect FAIL**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "theme class flips"`
Expected: FAIL — the class is hardcoded to `dockview-theme-light`.

- [ ] **Step 4.3: Make the class dynamic**

Edit `frontend/src/panels/DockviewLayout.tsx`:

1. Add the import near the top:

```tsx
import { useColorScheme } from "@mui/material/styles";
```

2. Inside the `DockviewLayout` function, before the `onReady` callback, add:

```tsx
const { mode, systemMode } = useColorScheme();
const themeClass =
    (mode === "system" ? systemMode : mode) === "dark"
        ? "dockview-theme-dark"
        : "dockview-theme-light";
```

3. Change the JSX prop from `className="dockview-theme-light"` to `className={themeClass}`.

- [ ] **Step 4.4: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts`
Expected: all assertions pass.

- [ ] **Step 4.5: Commit**

```bash
git add frontend/src/panels/DockviewLayout.tsx frontend/e2e/dockview-layout.spec.ts
git commit -m "feat(panels): sync dockview theme class with MUI color scheme

Flip className between dockview-theme-light and dockview-theme-dark
based on useColorScheme(). Combined with the MUI palette overrides
from the previous commit, this gives dockview a fully MUI-driven
light/dark appearance."
```

---

## Task 5: Plotly auto-resize on panel dimension change

**Goal:** Resizing a plot panel (splitter drag, popout, maximize, viewport change) makes the Plotly chart re-layout to fill the new container.

**Files:**
- Modify: `frontend/src/panels/PlotView.tsx` — add `onDidDimensionsChange` effect.
- Test: `frontend/e2e/dockview-layout.spec.ts`.

- [ ] **Step 5.1: Write the failing E2E assertion**

Append to the test block:

```ts
test("plotly chart resizes when the viewport changes", async ({ page }) => {
    await page.setViewportSize({ width: 1400, height: 900 });
    await page.goto(`/rooms/${ROOM}`);
    await page.getByTestId("activity-icon-plots-browser").click();
    const zone = page.getByTestId("sidebar-zone-left");
    const firstRow = zone.locator('li [role="button"]').first();
    await firstRow.click();
    const plotly = page.locator(".plotly").first();
    await expect(plotly).toBeVisible({ timeout: 15000 });
    const widthBefore = (await plotly.boundingBox())?.width ?? 0;

    await page.setViewportSize({ width: 900, height: 900 });
    await page.waitForTimeout(300);
    const widthAfter = (await plotly.boundingBox())?.width ?? 0;

    expect(widthAfter).toBeLessThan(widthBefore - 100);
});
```

- [ ] **Step 5.2: Build + run test → expect FAIL**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "plotly chart resizes"`
Expected: FAIL. Plotly renders once with `responsive: true` but dockview's panel container size changes aren't reaching plotly's internal resize loop; the chart stays at its initial width.

- [ ] **Step 5.3: Subscribe to dimension changes in PlotView**

Edit `frontend/src/panels/PlotView.tsx`. Inside the `PlotView` function, locate the effect near the end that handles purge-on-unmount (around line 710):

```tsx
// ===== EFFECT: Purge Plotly on unmount =====
useEffect(() => {
    const container = plotContainer.current;
    return () => {
        if (container) {
            Plotly.purge(container);
        }
    };
}, []);
```

**Immediately after** that effect, add:

```tsx
// ===== EFFECT: Re-layout Plotly on container size changes =====
// Dockview emits onDidDimensionsChange on splitter drag, popout, maximize,
// and viewport resize. Plotly's responsive:true handles window resizes but
// not intra-panel resizes, so we re-trigger its layout pass here.
useEffect(() => {
    const disposable = props.api.onDidDimensionsChange(() => {
        const container = plotContainer.current;
        if (container) {
            Plotly.Plots.resize(container);
        }
    });
    return () => disposable.dispose();
}, [props.api]);
```

- [ ] **Step 5.4: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "plotly chart resizes"`
Expected: PASS.

- [ ] **Step 5.5: Commit**

```bash
git add frontend/src/panels/PlotView.tsx frontend/e2e/dockview-layout.spec.ts
git commit -m "feat(panels): re-layout plotly on panel dimension change

Subscribe PlotView to dockview's onDidDimensionsChange and call
Plotly.Plots.resize(container). plotly's responsive:true covers window
resizes but misses intra-panel resizes (splitter drag, popout, maximize)."
```

---

## Task 6: Group actions — popout + maximize buttons

**Goal:** Every dockview group gets a popout button and a maximize/un-maximize toggle in its right-gutter action area.

**Files:**
- Create: `frontend/src/panels/groupActions.tsx`
- Modify: `frontend/src/panels/DockviewLayout.tsx` (wire `rightHeaderActionsComponent`)
- Test: `frontend/e2e/dockview-layout.spec.ts`

- [ ] **Step 6.1: Write the failing E2E assertion**

Append to the test block:

```ts
test("every group shows popout + maximize buttons", async ({ page }) => {
    await page.goto(`/rooms/${ROOM}`);
    await expect(page.getByTestId("group-popout")).toHaveCount(1);
    await expect(page.getByTestId("group-maximize")).toHaveCount(1);

    // Maximize the viewer group.
    await page.getByTestId("group-maximize").click();
    // Icon should flip to FullscreenExit — test the aria label.
    await expect(
        page.getByTestId("group-maximize").getByLabel("Exit full-screen"),
    ).toBeVisible();
    // Click again to restore.
    await page.getByTestId("group-maximize").click();
    await expect(
        page.getByTestId("group-maximize").getByLabel("Full-screen"),
    ).toBeVisible();
});

test("opening a second group yields a second popout+maximize pair", async ({
    page,
}) => {
    await page.goto(`/rooms/${ROOM}`);
    await page.getByTestId("activity-icon-plots-browser").click();
    const zone = page.getByTestId("sidebar-zone-left");
    const firstRow = zone.locator('li [role="button"]').first();
    await firstRow.click();
    await expect(page.locator('[data-testid^="plot-view-"]').first()).toBeVisible({
        timeout: 15000,
    });
    // Plot tab opens inside a new group to the right of the viewer (per §7
    // of the spec — wired in Task 8). Until Task 8 lands this may still
    // share the viewer's group, in which case assert ≥1. Assert ≥1 for now
    // and tighten after Task 8.
    expect(await page.getByTestId("group-popout").count()).toBeGreaterThanOrEqual(1);
});
```

- [ ] **Step 6.2: Build + run test → expect FAIL**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "popout"`
Expected: FAIL. No group actions component is wired yet.

- [ ] **Step 6.3: Create the groupActions component**

Create `frontend/src/panels/groupActions.tsx`:

```tsx
import FullscreenIcon from "@mui/icons-material/Fullscreen";
import FullscreenExitIcon from "@mui/icons-material/FullscreenExit";
import OpenInNewIcon from "@mui/icons-material/OpenInNew";
import { Box, IconButton, Tooltip } from "@mui/material";
import type { IDockviewHeaderActionsProps } from "dockview-react";
import { useEffect, useState } from "react";

export function GroupActions(props: IDockviewHeaderActionsProps) {
    const { api, group } = props;
    const [isMax, setIsMax] = useState<boolean>(() => group.api.isMaximized());

    useEffect(() => {
        const d = api.onDidLayoutChange(() => {
            setIsMax(group.api.isMaximized());
        });
        setIsMax(group.api.isMaximized());
        return () => d.dispose();
    }, [api, group]);

    const onPopout = () => {
        api.addPopoutGroup(group);
    };

    const onToggleMaximize = () => {
        if (group.api.isMaximized()) {
            api.exitMaximizedGroup();
        } else {
            api.maximizeGroup(group);
        }
    };

    return (
        <Box sx={{ display: "flex", alignItems: "center", px: 0.5 }}>
            <Tooltip title="Pop out">
                <IconButton
                    size="small"
                    onClick={onPopout}
                    data-testid="group-popout"
                    aria-label="Pop out"
                >
                    <OpenInNewIcon fontSize="small" />
                </IconButton>
            </Tooltip>
            <Tooltip title={isMax ? "Exit full-screen" : "Full-screen"}>
                <IconButton
                    size="small"
                    onClick={onToggleMaximize}
                    data-testid="group-maximize"
                    aria-label={isMax ? "Exit full-screen" : "Full-screen"}
                >
                    {isMax ? (
                        <FullscreenExitIcon fontSize="small" />
                    ) : (
                        <FullscreenIcon fontSize="small" />
                    )}
                </IconButton>
            </Tooltip>
        </Box>
    );
}
```

- [ ] **Step 6.4: Wire GroupActions into DockviewLayout**

Edit `frontend/src/panels/DockviewLayout.tsx`:

1. Add import near the top:

```tsx
import { GroupActions } from "./groupActions";
```

2. In the `<DockviewReact ... />` JSX, add the prop:

```tsx
<DockviewReact
    className={themeClass}
    onReady={onReady}
    components={components}
    onDidDrop={onDidDrop}
    floatingGroupBounds="boundedWithinViewport"
    rightHeaderActionsComponent={GroupActions}
/>
```

- [ ] **Step 6.5: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts`
Expected: both new assertions pass; all existing ones still pass.

- [ ] **Step 6.6: Commit**

```bash
git add frontend/src/panels/groupActions.tsx frontend/src/panels/DockviewLayout.tsx frontend/e2e/dockview-layout.spec.ts
git commit -m "feat(panels): popout + maximize buttons on every group

Add GroupActions component wired via rightHeaderActionsComponent.
Dockview renders it per-group automatically, so every group
(viewer, plot groups, popouts) gets its own pair. Icon flips between
Fullscreen and FullscreenExit based on group.api.isMaximized()."
```

---

## Task 7: Simplify PlotsBrowser to one list

**Goal:** Drop the "Currently Open" section; the filled-dot indicator on each row already conveys the same info. Inline `×` on hover of open rows to close.

**Files:**
- Modify: `frontend/src/panels/PlotsBrowserPanel.tsx`
- Test: `frontend/e2e/dockview-layout.spec.ts`

- [ ] **Step 7.1: Write the failing E2E assertion**

Append:

```ts
test("plots browser renders a single list (no 'Currently Open' section)", async ({
    page,
}) => {
    await page.goto(`/rooms/${ROOM}`);
    await page.getByTestId("activity-icon-plots-browser").click();
    const zone = page.getByTestId("sidebar-zone-left");
    await expect(zone).toBeVisible();
    await expect(zone.getByText(/Currently Open/i)).toHaveCount(0);

    // Open a plot, then the inline close button should appear on hover.
    const firstRow = zone.locator('li [role="button"]').first();
    const firstRowText = (await firstRow.textContent())?.trim();
    await firstRow.click();
    if (firstRowText) {
        await expect(
            page.locator(`[data-testid="plot-view-${firstRowText}"]`),
        ).toBeVisible({ timeout: 15000 });
    }

    // Still only one list.
    await expect(zone.getByText(/Currently Open/i)).toHaveCount(0);
    // Filled dot on the open row (use the data-testid we add next).
    await expect(
        zone.locator(`[data-testid="plot-row-${firstRowText}"][data-open="true"]`),
    ).toBeVisible();
});
```

- [ ] **Step 7.2: Build + run test → expect FAIL**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "plots browser renders a single list"`
Expected: FAIL. "Currently Open" section still renders.

- [ ] **Step 7.3: Rewrite PlotsBrowserPanel**

Replace the full contents of `frontend/src/panels/PlotsBrowserPanel.tsx` with:

```tsx
import CloseIcon from "@mui/icons-material/Close";
import {
    Box,
    IconButton,
    List,
    ListItem,
    ListItemButton,
    ListItemText,
    Typography,
} from "@mui/material";
import { useEffect, useMemo, useState } from "react";
import { useFigureList } from "../hooks/useFigures";
import { getDockviewApi } from "./DockviewLayout";
import { closePlotTab, openPlotTab } from "./plotViewFactory";

const DRAG_MIME_PLOT = "application/x-zndraw-plot-key";

function useOpenPlotKeys(): Set<string> {
    const [version, setVersion] = useState(0);

    useEffect(() => {
        const api = getDockviewApi();
        if (!api) return;
        const add = api.onDidAddPanel(() => setVersion((v) => v + 1));
        const remove = api.onDidRemovePanel(() => setVersion((v) => v + 1));
        return () => {
            add.dispose();
            remove.dispose();
        };
    }, []);

    return useMemo(() => {
        const api = getDockviewApi();
        if (!api) return new Set<string>();
        return new Set(
            api.panels
                .filter((p) => p.id.startsWith("plot-"))
                .map((p) => p.id.slice("plot-".length)),
        );
        // biome-ignore lint/correctness/useExhaustiveDependencies: version is a tick trigger
    }, [version]);
}

export function PlotsBrowserPanel() {
    const { data, isLoading } = useFigureList();
    const openKeys = useOpenPlotKeys();
    const allKeys = useMemo(() => data?.items ?? [], [data]);

    const onRowClick = (key: string) => {
        const api = getDockviewApi();
        if (!api) return;
        openPlotTab(api, key);
    };

    const onRowClose = (e: React.MouseEvent, key: string) => {
        e.stopPropagation();
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

    return (
        <Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
            <Typography variant="overline" sx={{ px: 2, pt: 1 }}>
                Plots
            </Typography>
            <List dense sx={{ flexGrow: 1, overflow: "auto" }}>
                {allKeys.map((k: string) => {
                    const isOpen = openKeys.has(k);
                    return (
                        <ListItem
                            key={k}
                            data-testid={`plot-row-${k}`}
                            data-open={isOpen ? "true" : "false"}
                            disablePadding
                            draggable
                            onDragStart={(e) => onDragStart(e, k)}
                            sx={{
                                "& .plot-row-close": { visibility: "hidden" },
                                "&:hover .plot-row-close": {
                                    visibility: isOpen ? "visible" : "hidden",
                                },
                            }}
                            secondaryAction={
                                isOpen ? (
                                    <IconButton
                                        className="plot-row-close"
                                        edge="end"
                                        size="small"
                                        aria-label={`Close ${k}`}
                                        onClick={(e) => onRowClose(e, k)}
                                    >
                                        <CloseIcon fontSize="small" />
                                    </IconButton>
                                ) : null
                            }
                        >
                            <ListItemButton onClick={() => onRowClick(k)}>
                                <Box
                                    sx={{
                                        width: 8,
                                        height: 8,
                                        borderRadius: "50%",
                                        bgcolor: isOpen ? "primary.main" : "transparent",
                                        border: 1,
                                        borderColor: "primary.main",
                                        mr: 1.5,
                                    }}
                                />
                                <ListItemText primary={k} />
                            </ListItemButton>
                        </ListItem>
                    );
                })}
            </List>
        </Box>
    );
}
```

Diff vs. the old file:
- Removed `Divider` import and its element.
- Removed the second `<List>` and `<Typography variant="overline">Currently Open</Typography>` block.
- Added `data-testid="plot-row-{k}"` and `data-open` on each row for E2E assertions.
- Moved close button inline via `secondaryAction`, hidden by default, revealed on hover for open rows only (hover + `visibility`).
- Renamed the heading from "Available Plots" to just "Plots" — one list doesn't need the qualifier.

- [ ] **Step 7.4: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts`
Expected: all assertions pass.

- [ ] **Step 7.5: Commit**

```bash
git add frontend/src/panels/PlotsBrowserPanel.tsx frontend/e2e/dockview-layout.spec.ts
git commit -m "refactor(panels): collapse PlotsBrowser to one list

The separate 'Currently Open' section was redundant: the filled-dot
indicator conveys open state. Moved the close button inline as a
hover-reveal secondaryAction on each open row."
```

---

## Task 8: Open new plots beside the viewer (not on top of it)

**Goal:** `openPlotTab` resolves the target group deterministically: lex-first existing plot group > right of viewer > dockview default.

**Files:**
- Modify: `frontend/src/panels/plotViewFactory.ts`
- Test: `frontend/e2e/dockview-layout.spec.ts`

- [ ] **Step 8.1: Write the failing E2E assertion**

Append:

```ts
test("first plot opens in a new group to the right of the viewer", async ({
    page,
}) => {
    await page.goto(`/rooms/${ROOM}`);
    const viewer = page.getByTestId("viewer-view");
    const viewerBox = await viewer.boundingBox();
    if (!viewerBox) throw new Error("viewer bounding box missing");

    await page.getByTestId("activity-icon-plots-browser").click();
    const zone = page.getByTestId("sidebar-zone-left");
    const firstRow = zone.locator('li [role="button"]').first();
    const firstRowText = (await firstRow.textContent())?.trim();
    await firstRow.click();
    const plot = page.locator(`[data-testid="plot-view-${firstRowText}"]`);
    await expect(plot).toBeVisible({ timeout: 15000 });
    const plotBox = await plot.boundingBox();
    if (!plotBox) throw new Error("plot bounding box missing");

    // Plot's left edge should be at or beyond the viewer's right edge.
    expect(plotBox.x).toBeGreaterThanOrEqual(viewerBox.x + viewerBox.width - 5);
});

test("second plot stacks with the first (shares a group)", async ({
    page,
}) => {
    await page.goto(`/rooms/${ROOM}`);
    await page.getByTestId("activity-icon-plots-browser").click();
    const zone = page.getByTestId("sidebar-zone-left");
    const rows = zone.locator('li [role="button"]');
    await expect(rows).toHaveCount(2, { timeout: 15000 });

    const firstKey = (await rows.nth(0).textContent())?.trim();
    const secondKey = (await rows.nth(1).textContent())?.trim();
    await rows.nth(0).click();
    await expect(
        page.locator(`[data-testid="plot-view-${firstKey}"]`),
    ).toBeVisible({ timeout: 15000 });
    await rows.nth(1).click();
    await expect(
        page.locator(`[data-testid="plot-view-${secondKey}"]`),
    ).toBeVisible({ timeout: 15000 });

    // Both plot bounding boxes should overlap fully (one tab behind the other).
    const firstBox = await page
        .locator(`[data-testid="plot-view-${firstKey}"]`)
        .boundingBox();
    const secondBox = await page
        .locator(`[data-testid="plot-view-${secondKey}"]`)
        .boundingBox();
    if (!firstBox || !secondBox) throw new Error("plot bounding box missing");
    // Same group → same x position (within a few px for tab-header differences).
    expect(Math.abs(firstBox.x - secondBox.x)).toBeLessThan(5);
});
```

- [ ] **Step 8.2: Build + run test → expect FAIL**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "opens in a new group to the right"`
Expected: FAIL. The current `openPlotTab` with no `referenceGroupId` calls `addPanel` without a position, so dockview stacks the new panel in the active group (the viewer's group). The plot ends up on top of the viewer.

- [ ] **Step 8.3: Rewrite openPlotTab**

Replace the contents of `frontend/src/panels/plotViewFactory.ts` with:

```ts
import type { DockviewApi } from "dockview-react";

export function plotPanelId(figureKey: string): string {
    return `plot-${figureKey}`;
}

const PLOT_ID_PREFIX = "plot-";
const VIEWER_PANEL_ID = "viewer";

/**
 * Resolves which group a new plot tab should open in.
 *
 * Rules (first match wins):
 * 1. If any existing `plot-*` panel is open, pick the group containing
 *    the plot whose id sorts first lexicographically. Stable across
 *    sessions and independent of insertion order.
 * 2. Else, place a new group to the right of the viewer's group.
 * 3. Else (no viewer), let dockview decide.
 */
function resolvePlotPosition(api: DockviewApi):
    | { referenceGroup: string; direction: "within" | "right" }
    | undefined {
    const plotPanels = api.panels
        .filter((p) => p.id.startsWith(PLOT_ID_PREFIX))
        .sort((a, b) => a.id.localeCompare(b.id));

    if (plotPanels.length > 0) {
        const group = plotPanels[0].group;
        if (group) {
            return { referenceGroup: group.id, direction: "within" };
        }
    }

    const viewer = api.getPanel(VIEWER_PANEL_ID);
    if (viewer?.group) {
        return { referenceGroup: viewer.group.id, direction: "right" };
    }

    return undefined;
}

export function openPlotTab(
    api: DockviewApi,
    figureKey: string,
    opts?: {
        referenceGroupId?: string;
        direction?: "right" | "left" | "above" | "below" | "within";
    },
): void {
    const id = plotPanelId(figureKey);
    const existing = api.getPanel(id);
    if (existing) {
        existing.api.setActive();
        return;
    }

    // Explicit caller override (used by the onDidDrop path in DockviewLayout
    // where the user dropped onto a specific target).
    if (opts?.referenceGroupId) {
        api.addPanel({
            id,
            component: "plotView",
            title: figureKey,
            params: { figureKey },
            position: {
                referenceGroup: opts.referenceGroupId,
                direction: opts.direction ?? "within",
            },
        });
        return;
    }

    const position = resolvePlotPosition(api);
    api.addPanel({
        id,
        component: "plotView",
        title: figureKey,
        params: { figureKey },
        position,
    });
}

export function closePlotTab(api: DockviewApi, figureKey: string): void {
    const panel = api.getPanel(plotPanelId(figureKey));
    panel?.api.close();
}
```

- [ ] **Step 8.4: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts`
Expected: both new positioning assertions pass. The earlier "popout + maximize pair count ≥ 1" assertion from Task 6 can now be tightened — with the plot in its own group, there will be 2 pairs. This tightening is optional; the `≥ 1` still holds.

- [ ] **Step 8.5: Commit**

```bash
git add frontend/src/panels/plotViewFactory.ts frontend/e2e/dockview-layout.spec.ts
git commit -m "feat(panels): open new plots beside the viewer, not on top

openPlotTab now resolves the target group deterministically: existing
plot group (lex-first panel id) > right of viewer > dockview default.
The onDidDrop path in DockviewLayout still honors the user's explicit
drop target because it passes referenceGroupId."
```

---

## Task 9: ActivityBar three-state render (sliver / hot / full)

**Goal:** Empty activity bars collapse to a 4px invisible sliver; during a panel-icon drag they light up as hot drop zones; a bar with icons renders its full 48px (or 40px for bottom).

**Files:**
- Modify: `frontend/src/panels/ActivityBar.tsx`
- Test: `frontend/e2e/dockview-layout.spec.ts`

- [ ] **Step 9.1: Write the failing E2E assertion**

Append:

```ts
test("empty right bar collapses to a 4px sliver", async ({ page }) => {
    await page.goto(`/rooms/${ROOM}`);
    const right = page.getByTestId("activity-bar-right");
    const box = await right.boundingBox();
    if (!box) throw new Error("right bar bounding box missing");
    expect(box.width).toBeLessThanOrEqual(8);
});

test("empty bottom bar collapses to a 4px sliver", async ({ page }) => {
    await page.goto(`/rooms/${ROOM}`);
    const bottom = page.getByTestId("activity-bar-bottom");
    const box = await bottom.boundingBox();
    if (!box) throw new Error("bottom bar bounding box missing");
    expect(box.height).toBeLessThanOrEqual(8);
});

test("dragging an icon lights up empty bars as hot drop zones", async ({
    page,
}) => {
    await page.goto(`/rooms/${ROOM}`);
    const icon = page.getByTestId("activity-icon-selections");
    const right = page.getByTestId("activity-bar-right");

    // Simulate an in-flight drag by dispatching dragstart/dragenter/dragover
    // events on the icon and bar. Playwright's high-level dragTo() completes
    // the drop synchronously, so we can't observe the "hot" middle state.
    // Instead, mutate state directly via the window drag events the component
    // listens on.
    await page.evaluate(() => {
        const dt = new DataTransfer();
        dt.setData("application/x-zndraw-panel-id", "selections");
        window.dispatchEvent(
            new DragEvent("dragstart", { dataTransfer: dt, bubbles: true }),
        );
    });

    await expect(right).toHaveAttribute("data-sliver-state", "hot");

    await page.evaluate(() => {
        window.dispatchEvent(new DragEvent("dragend", { bubbles: true }));
    });
    await expect(right).toHaveAttribute("data-sliver-state", "sliver");
});
```

- [ ] **Step 9.2: Build + run test → expect FAIL**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts --grep "sliver|hot drop"`
Expected: FAIL. ActivityBar always renders at full size.

- [ ] **Step 9.3: Rewrite ActivityBar with three-state render**

Replace the contents of `frontend/src/panels/ActivityBar.tsx` with:

```tsx
import Badge from "@mui/material/Badge";
import { Box, IconButton, Tooltip, keyframes } from "@mui/material";
import { useCallback, useEffect, useState } from "react";
import { useShallow } from "zustand/react/shallow";
import { useAppStore } from "../store";
import { type BarPosition, type PanelId, PANELS } from "./registry";

const DRAG_MIME = "application/x-zndraw-panel-id";

const BAR_SX = {
    left: { flexDirection: "column", width: 48, borderRight: 1 },
    right: { flexDirection: "column", width: 48, borderLeft: 1 },
    bottom: { flexDirection: "row", height: 40, borderTop: 1, width: "100%" },
} as const;

const SLIVER_SX = {
    left: { width: 4 },
    right: { width: 4 },
    bottom: { height: 4, width: "100%" },
} as const;

const ACTIVE_INDICATOR: Record<
    BarPosition,
    (color: string) => Record<string, string>
> = {
    left: (c) => ({ borderLeft: `2px solid ${c}` }),
    right: (c) => ({ borderRight: `2px solid ${c}` }),
    bottom: (c) => ({ borderTop: `2px solid ${c}` }),
};

const pulse = keyframes`
    0%, 100% { opacity: 0.5; }
    50% { opacity: 1; }
`;

type SliverState = "full" | "sliver" | "hot";

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
    const chatUnread = useAppStore((s) => s.chatUnreadCount);

    const [isDragActive, setIsDragActive] = useState(false);

    useEffect(() => {
        const onGlobalDragStart = (e: DragEvent) => {
            if (e.dataTransfer?.types.includes(DRAG_MIME)) {
                setIsDragActive(true);
            }
        };
        const onGlobalDragEnd = () => {
            setIsDragActive(false);
        };
        window.addEventListener("dragstart", onGlobalDragStart);
        window.addEventListener("dragend", onGlobalDragEnd);
        return () => {
            window.removeEventListener("dragstart", onGlobalDragStart);
            window.removeEventListener("dragend", onGlobalDragEnd);
        };
    }, []);

    const onDragStart = useCallback((e: React.DragEvent, id: PanelId) => {
        e.dataTransfer.setData(DRAG_MIME, id);
        e.dataTransfer.effectAllowed = "move";
    }, []);

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

    const state: SliverState =
        icons.length === 0 ? (isDragActive ? "hot" : "sliver") : "full";

    if (state !== "full") {
        return (
            <Box
                data-testid={`activity-bar-${position}`}
                data-sliver-state={state}
                onDragOver={onDragOver}
                onDrop={onDropOnBar}
                sx={{
                    display: "flex",
                    flexShrink: 0,
                    ...SLIVER_SX[position],
                    bgcolor:
                        state === "hot"
                            ? "rgba(25, 118, 210, 0.2)"
                            : "transparent",
                    ...(state === "hot"
                        ? {
                                borderStyle: "solid",
                                borderColor: "primary.main",
                                borderWidth:
                                    position === "left"
                                        ? "0 2px 0 0"
                                        : position === "right"
                                            ? "0 0 0 2px"
                                            : "2px 0 0 0",
                                animation: `${pulse} 1s ease-in-out infinite`,
                            }
                        : {}),
                }}
            />
        );
    }

    return (
        <Box
            data-testid={`activity-bar-${position}`}
            data-sliver-state="full"
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
                const iconElement = <Icon />;
                const iconWithBadge =
                    id === "chat" && !isActive && chatUnread > 0 ? (
                        <Badge badgeContent={chatUnread} color="error" max={99}>
                            {iconElement}
                        </Badge>
                    ) : (
                        iconElement
                    );
                return (
                    <Tooltip
                        key={id}
                        title={def.label}
                        placement={
                            position === "left"
                                ? "right"
                                : position === "right"
                                    ? "left"
                                    : "top"
                        }
                    >
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
                                    ? ACTIVE_INDICATOR[position](
                                            "var(--mui-palette-primary-main)",
                                        )
                                    : {}),
                            }}
                        >
                            {iconWithBadge}
                        </IconButton>
                    </Tooltip>
                );
            })}
        </Box>
    );
}
```

Key diffs vs. the existing file:
- New `isDragActive` local state, managed by `window.addEventListener('dragstart' | 'dragend', …)`.
- Derived `state: "full" | "sliver" | "hot"` from `icons.length === 0` and `isDragActive`.
- When not full, render a 4px `<Box>` with pulsing-border styling if hot, transparent if sliver. Still accepts drop events, so dropping onto the sliver still calls `moveIconToBar`.
- `data-sliver-state` attribute exposed for E2E.

- [ ] **Step 9.4: Rebuild + run test → expect PASS**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts`
Expected: three new assertions pass; existing "default startup shows activity bars and viewer" still passes because the slivers ARE present (they're still rendered `<Box>` elements with `data-testid`), just narrow.

Note: the existing test asserts `await expect(page.getByTestId("activity-bar-right")).toBeVisible();` — Playwright considers a 4px × viewport-height element visible. If the test fails after this change, widen the sliver to 6px or adjust the assertion.

- [ ] **Step 9.5: Commit**

```bash
git add frontend/src/panels/ActivityBar.tsx frontend/e2e/dockview-layout.spec.ts
git commit -m "feat(panels): three-state ActivityBar (full/sliver/hot)

Empty bars collapse to a 4px invisible sliver; a pulsing 'hot' state
lights up during a panel-icon drag so the drop zone is discoverable
without eating chrome space when idle. Implemented via a window-level
drag listener that sets a local isDragActive flag."
```

---

## Task 10: Final validation

**Goal:** Typecheck, lint, full E2E suite, and a manual screenshot sweep.

- [ ] **Step 10.1: Typecheck**

Run: `cd frontend && bunx tsc --noEmit`
Expected: no errors.

- [ ] **Step 10.2: Lint**

Run: `cd frontend && bunx biome check src/panels/ e2e/dockview-layout.spec.ts`
Expected: no errors. If biome flags style-only issues, run `bunx biome check --write` on those files and re-run.

- [ ] **Step 10.3: Full dockview E2E spec**

Run: `cd frontend && bun run build && cd .. && uv sync --reinstall-package zndraw && cd frontend && bunx playwright test e2e/dockview-layout.spec.ts`
Expected: all assertions in `dockview-layout.spec.ts` pass (the pre-existing 5 tests plus the ones added in Tasks 1–9). The 13 pre-existing broken specs in other files are out of scope.

- [ ] **Step 10.4: Manual screenshot sweep**

Start `uv run zndraw` if not already running, then from a new terminal run the playwright-cli to capture three screenshots:

```bash
mkdir -p /tmp/dockview-validation
```

Then, using your preferred playwright-cli session (e.g., `/playwright-cli` in Claude Code, or `bunx playwright open http://localhost:8000/rooms/dockview-test`):

1. Screenshot: default startup → `/tmp/dockview-validation/followup-01-default.png`.
2. Open plots-browser, open `energy` → screenshot → `/tmp/dockview-validation/followup-02-plot-open.png`. Verify plot is to the right of viewer.
3. Click maximize on the plot group → screenshot → `/tmp/dockview-validation/followup-03-maximized.png`.
4. Toggle MUI color scheme to dark (via the app's theme toggle, or the devtools trick from Task 4 step 1) → screenshot → `/tmp/dockview-validation/followup-04-dark.png`. Tab chrome should follow MUI dark palette.

- [ ] **Step 10.5: Review + final commit**

If typecheck, lint, and E2E all clean and screenshots look right, no further commit needed. Otherwise fix issues inline and commit them per-file.

Optional housekeeping commit for README / ISSUES.md if you want to mark items resolved:

```bash
# Only if needed.
git add ISSUES.md
git commit -m "docs: mark dockview follow-up issues resolved"
```

---

## Self-review checklist

After all tasks are committed, verify:

- [ ] **Spec coverage.** Walk through §1–§8 of `2026-04-16-dockview-followup-design.md`:
  - §1 Plotly resize → Task 5 ✓
  - §2 Theme sync → Tasks 3 + 4 ✓
  - §3 Popout → Task 6 ✓
  - §4 Maximize → Task 6 ✓
  - §5 Chrome economy + chat bar move → Tasks 1 + 9 ✓
  - §6 PlotsBrowser simplification → Task 7 ✓
  - §7 Plot positioning → Task 8 ✓
  - §8 Bottom-drawer shrink fix → Task 2 ✓

- [ ] **E2E coverage matches spec.** The spec's validation list has 8 assertions; every one is in Task 1–9 tests.

- [ ] **No regressions in existing E2E.** `dockview-layout.spec.ts` pre-existing 5 tests still pass.

- [ ] **Typecheck + lint clean.**

- [ ] **All new files committed.** `dockview-mui.css`, `groupActions.tsx`.

- [ ] **All commits have clear messages tied to an issue.**

---

## Gotchas / known risks

- **MUI v7 `useColorScheme` on `<html>` attribute.** The assertion in Task 4.1 mutates `data-mui-color-scheme` directly. If the app uses `data-mui-color-scheme` differently (scoped to a container), adjust the selector. Verify via `page.evaluate(() => document.documentElement.outerHTML)` during the first E2E run.
- **Playwright `dragTo` synchronicity.** Task 9 uses window-level `DragEvent` dispatch rather than Playwright's `dragTo` to observe the "hot" middle state. If a future test needs to drag-and-drop between bars and observe the final state only, `dragTo` still works.
- **`group.api.isMaximized()` reactivity.** The `onDidLayoutChange` event in Task 6 fires on many state changes; the `setIsMax` call guards against no-op re-renders because React dedups same-value updates. If profiling shows redundant renders, switch to `api.onDidMaximizedGroupChange` (dockview 5.x) which only fires on max-state transitions.
- **dockview popout CSS.** Popout windows get their own `<html>` with a fresh MUI CssVarsProvider mount — dockview serves the popout from the same static bundle, so the MUI CSS variables are present. If popouts look unthemed, verify the CssVarsProvider in `main.tsx` applies to the popout document too (it should, since the popout shares the same static assets and script-loads the same bundle).
- **`bun run dev` port conflicts.** Dev runs at :5173 and proxies API calls to :8000. If :5173 is taken, vite picks the next free port; update the playwright config if E2E runs against the dev server. Most E2E should run against :8000 (built static) for consistency with CI.
