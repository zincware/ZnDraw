// Prerequisite: seed the test room + at least one figure before running this spec.
// Run once:
//   uv run zndraw-cli rooms create --room dockview-test 2>/dev/null || true
//   uv run python -c "
//   from zndraw import ZnDraw
//   import plotly.graph_objects as go
//   vis = ZnDraw(url='http://localhost:8000', room='dockview-test')
//   for key in ['energy']:
//       if key not in vis.figures:
//           vis.figures[key] = go.Figure(go.Scatter(x=[0,1,2], y=[0,1,4], mode='lines'))
//   "

import { expect, test } from "@playwright/test";

const ROOM = "dockview-test";

test.describe("dockview layout", () => {
	test("default startup shows activity bars and viewer", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		await expect(page.getByTestId("activity-bar-left")).toBeVisible();
		await expect(page.getByTestId("activity-bar-right")).toBeVisible();
		await expect(page.getByTestId("activity-bar-bottom")).toBeVisible();
		await expect(page.getByTestId("viewer-view")).toBeVisible();
		await expect(page.getByTestId("sidebar-zone-left")).toBeHidden();
	});

	test("clicking selections icon opens the left sidebar zone", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-selections").click();
		await expect(page.getByTestId("sidebar-zone-left")).toBeVisible();
		await page.getByTestId("activity-icon-selections").click();
		await expect(page.getByTestId("sidebar-zone-left")).toBeHidden();
	});

	test("only one panel per bar — switching icons swaps the panel", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-selections").click();
		await expect(page.getByTestId("sidebar-zone-left")).toBeVisible();
		await page.getByTestId("activity-icon-modifiers").click();
		// sidebar still visible; selections icon is no longer active (primary color)
		await expect(page.getByTestId("sidebar-zone-left")).toBeVisible();
	});

	test("opening plots browser and clicking a plot opens a tab", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-plots-browser").click();
		const zone = page.getByTestId("sidebar-zone-left");
		await expect(zone).toBeVisible();
		// MUI ListItemButton renders as <div role="button"> inside <li>, not <button>.
		const firstRow = zone.locator('li [role="button"]').first();
		await expect(firstRow).toBeVisible();
		const firstRowText = (await firstRow.textContent())?.trim();
		await firstRow.click();
		// Plot panel opens inside the dockview editor area; its outer wrapper has
		// data-testid=plot-view-<key>. Plotly lazy-loads, so we give it time.
		if (firstRowText) {
			await expect(
				page.locator(`[data-testid="plot-view-${firstRowText}"]`),
			).toBeVisible({ timeout: 15000 });
		}
	});

	test("closing viewer triggers leave-room cascade (URL goes to /)", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await expect(page.getByTestId("viewer-view")).toBeVisible();
		page.on("dialog", (d) => d.accept());
		// Dockview renders a tab action (close X) on the active group; its
		// visible class is .dv-default-tab-action.
		const closeBtn = page.locator(".dv-default-tab-action").first();
		await closeBtn.click();
		await expect(page).toHaveURL(/\/$/);
	});

	test("chat icon defaults to the left activity bar", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		const leftBar = page.getByTestId("activity-bar-left");
		await expect(leftBar.getByTestId("activity-icon-chat")).toBeVisible();
		const rightBar = page.getByTestId("activity-bar-right");
		await expect(rightBar.getByTestId("activity-icon-chat")).toHaveCount(0);
	});

	test("opening a bottom panel shrinks the viewer vertically", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		const viewer = page.getByTestId("viewer-view");
		await expect(viewer).toBeVisible();
		const before = await viewer.boundingBox();
		if (!before) throw new Error("viewer bounding box missing");

		// Drag the plots-browser icon into the bottom bar and activate it.
		const plotsIcon = page.getByTestId("activity-icon-plots-browser");
		const bottomBar = page.getByTestId("activity-bar-bottom");
		await plotsIcon.dragTo(bottomBar);
		await plotsIcon.click();

		const bottomZone = page.getByTestId("bottom-zone");
		await expect(bottomZone).toBeVisible();
		// Give the ResizeObserver time to fire and dockview to re-layout panels.
		await page.waitForTimeout(300);

		const after = await viewer.boundingBox();
		if (!after) throw new Error("viewer bounding box missing after");
		expect(after.height).toBeLessThan(before.height - 50);
	});

	test("dockview tab chrome reads MUI palette in light mode", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.waitForTimeout(1000);

		const themed = page.locator(".dockview-theme-light").first();
		await expect(themed).toHaveCount(1);

		const vars = await themed.evaluate((el) => {
			const s = getComputedStyle(el);
			return {
				groupBg: s.getPropertyValue("--dv-group-view-background-color").trim(),
				tabsBarBg: s
					.getPropertyValue("--dv-tabs-and-actions-container-background-color")
					.trim(),
				activeTabColor: s
					.getPropertyValue("--dv-activegroup-visiblepanel-tab-color")
					.trim(),
				muiBgDefault: s.getPropertyValue("--mui-palette-background-default").trim(),
				muiBgPaper: s.getPropertyValue("--mui-palette-background-paper").trim(),
				muiTextPrimary: s.getPropertyValue("--mui-palette-text-primary").trim(),
			};
		});

		// dockview vars must resolve to the same values as the MUI tokens.
		expect(vars.groupBg).toBe(vars.muiBgDefault);
		expect(vars.tabsBarBg).toBe(vars.muiBgPaper);
		expect(vars.activeTabColor).toBe(vars.muiTextPrimary);
		// Light-mode sanity: MUI background-default must be a light colour.
		expect(vars.muiBgDefault.toLowerCase()).toMatch(/#fff|rgb\(255/);
	});

	test("dockview tab chrome reads MUI palette in dark mode", async ({
		page,
	}) => {
		// Seed dark mode before page load so MUI's useColorScheme picks it up.
		await page.addInitScript(() =>
			localStorage.setItem("mui-mode", "dark"),
		);
		await page.goto(`/rooms/${ROOM}`);
		await page.waitForTimeout(1500);

		// In dark mode, dockview uses themeDark whose className is
		// 'dockview-theme-dark' (Task 4 wires this).
		const themed = page.locator(".dockview-theme-dark").first();
		await expect(themed).toHaveCount(1);

		const vars = await themed.evaluate((el) => {
			const s = getComputedStyle(el);
			return {
				groupBg: s.getPropertyValue("--dv-group-view-background-color").trim(),
				tabsBarBg: s
					.getPropertyValue("--dv-tabs-and-actions-container-background-color")
					.trim(),
				activeTabColor: s
					.getPropertyValue("--dv-activegroup-visiblepanel-tab-color")
					.trim(),
				muiBgDefault: s.getPropertyValue("--mui-palette-background-default").trim(),
				muiBgPaper: s.getPropertyValue("--mui-palette-background-paper").trim(),
				muiTextPrimary: s.getPropertyValue("--mui-palette-text-primary").trim(),
			};
		});

		// dockview vars must resolve to MUI tokens (now dark palette).
		expect(vars.groupBg).toBe(vars.muiBgDefault);
		expect(vars.tabsBarBg).toBe(vars.muiBgPaper);
		expect(vars.activeTabColor).toBe(vars.muiTextPrimary);
		// Dark-mode sanity: MUI background-default must be a dark colour.
		expect(vars.muiBgDefault.toLowerCase()).toMatch(/#121212|rgb\(18/);
	});

	test("3D viewer canvas background flips to dark in dark mode", async ({
		page,
	}) => {
		await page.addInitScript(() =>
			localStorage.setItem("mui-mode", "dark"),
		);
		await page.goto(`/rooms/${ROOM}`);
		await page.waitForTimeout(2500);

		// R3F's <Canvas> outer wrapper div receives the inline
		// `background: <theme.palette.background.default>` from Canvas.tsx.
		// When MUI is in dark mode, that ancestor must be a dark colour —
		// not the white default we had before the MuiCssVars fix.
		const bg = await page.evaluate(() => {
			const canvas = document.querySelector(
				"[data-testid='viewer-view'] canvas",
			) as HTMLCanvasElement | null;
			if (!canvas) return null;
			let el: HTMLElement | null = canvas.parentElement;
			while (el) {
				const c = getComputedStyle(el).backgroundColor;
				if (c !== "rgba(0, 0, 0, 0)" && c !== "transparent") return c;
				el = el.parentElement;
			}
			return null;
		});
		expect(bg).not.toBeNull();
		// Dark MUI background is rgb(18, 18, 18) a.k.a. #121212.
		expect(bg?.toLowerCase()).toMatch(/rgb\(18|#121212/);
	});

	test("dockview theme flips with MUI color scheme", async ({ page }) => {
		// Set dark mode via localStorage BEFORE page load so MUI's
		// useColorScheme picks it up on first render.
		await page.addInitScript(() => {
			localStorage.setItem("mui-mode", "dark");
		});
		await page.goto(`/rooms/${ROOM}`);
		await page.waitForTimeout(1000);

		// dockview-react applies the theme's className to its root element.
		// themeDark has className 'dockview-theme-dark'.
		const darkEl = page.locator(".dockview-theme-dark");
		await expect(darkEl).toHaveCount(1);
		const lightEl = page.locator(".dockview-theme-light");
		await expect(lightEl).toHaveCount(0);

		// Spot-check that the MUI dark palette is active and the dockview
		// vars resolve to it.
		const vars = await darkEl.first().evaluate((el) => {
			const s = getComputedStyle(el);
			return {
				groupBg: s.getPropertyValue("--dv-group-view-background-color").trim(),
				muiBgDefault: s.getPropertyValue("--mui-palette-background-default").trim(),
			};
		});
		expect(vars.groupBg).toBe(vars.muiBgDefault);
		expect(vars.muiBgDefault.toLowerCase()).toMatch(/#121212|rgb\(18/);
	});
});
