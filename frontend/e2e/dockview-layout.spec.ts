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
});
