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
		await page.waitForTimeout(500);
		const widthAfter = (await plotly.boundingBox())?.width ?? 0;

		expect(widthAfter).toBeLessThan(widthBefore - 100);
	});

	test("every group shows popout + maximize buttons", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		await expect(page.getByTestId("group-popout")).toHaveCount(1);
		await expect(page.getByTestId("group-maximize")).toHaveCount(1);

		// Maximize the viewer group and verify the icon flips.
		await page.getByTestId("group-maximize").click();
		await expect(page.getByTestId("group-maximize")).toHaveAttribute(
			"aria-label",
			"Exit full-screen",
		);
		// Click again to restore.
		await page.getByTestId("group-maximize").click();
		await expect(page.getByTestId("group-maximize")).toHaveAttribute(
			"aria-label",
			"Full-screen",
		);
	});

	test("plots browser renders a single list (no 'Currently Open' section)", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-plots-browser").click();
		const zone = page.getByTestId("sidebar-zone-left");
		await expect(zone).toBeVisible();
		await expect(zone.getByText(/Currently Open/i)).toHaveCount(0);

		const firstRow = zone.locator('li [role="button"]').first();
		const firstRowText = (await firstRow.textContent())?.trim();
		await firstRow.click();
		if (firstRowText) {
			await expect(
				page.locator(`[data-testid="plot-view-${firstRowText}"]`),
			).toBeVisible({ timeout: 15000 });
		}

		// Still only one list — no "Currently Open" after opening.
		await expect(zone.getByText(/Currently Open/i)).toHaveCount(0);
		// Filled dot on the open row (via data-open="true" on the row).
		await expect(
			zone.locator(`[data-testid="plot-row-${firstRowText}"][data-open="true"]`),
		).toBeVisible();
	});

	test("first plot opens in a new group to the right of the viewer", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		const viewer = page.getByTestId("viewer-view");

		await page.getByTestId("activity-icon-plots-browser").click();
		const zone = page.getByTestId("sidebar-zone-left");
		const firstRow = zone.locator('li [role="button"]').first();
		const firstRowText = (await firstRow.textContent())?.trim();
		await firstRow.click();
		const plot = page.locator(`[data-testid="plot-view-${firstRowText}"]`);
		await expect(plot).toBeVisible({ timeout: 15000 });

		// Re-measure both after the plot opens so the layout is stable.
		const viewerBox = await viewer.boundingBox();
		const plotBox = await plot.boundingBox();
		if (!viewerBox) throw new Error("viewer bounding box missing");
		if (!plotBox) throw new Error("plot bounding box missing");

		// Plot's left edge should be at/beyond the viewer's right edge.
		expect(plotBox.x).toBeGreaterThanOrEqual(viewerBox.x + viewerBox.width - 5);
	});

	test("second plot stacks with the first (shares a group)", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-plots-browser").click();
		const zone = page.getByTestId("sidebar-zone-left");
		const rows = zone.locator('li [role="button"]');
		await expect(rows.first()).toBeVisible({ timeout: 15000 });
		const count = await rows.count();
		expect(count).toBeGreaterThanOrEqual(2);

		const firstKey = (await rows.nth(0).textContent())?.trim();
		const secondKey = (await rows.nth(1).textContent())?.trim();
		await rows.nth(0).click();
		await expect(
			page.locator(`[data-testid="plot-view-${firstKey}"]`),
		).toBeVisible({ timeout: 15000 });
		// Capture first plot x while it is the active (visible) tab.
		const firstBox = await page
			.locator(`[data-testid="plot-view-${firstKey}"]`)
			.boundingBox();
		if (!firstBox) throw new Error("first plot bounding box missing");

		await rows.nth(1).click();
		await expect(
			page.locator(`[data-testid="plot-view-${secondKey}"]`),
		).toBeVisible({ timeout: 15000 });
		const secondBox = await page
			.locator(`[data-testid="plot-view-${secondKey}"]`)
			.boundingBox();
		if (!secondBox) throw new Error("second plot bounding box missing");

		// Same group → panels share x-position (within a few px for tab edges).
		expect(Math.abs(firstBox.x - secondBox.x)).toBeLessThan(5);
	});

	test("empty right bar collapses to a sliver (≤8px wide)", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.waitForTimeout(500);
		const right = page.getByTestId("activity-bar-right");
		const box = await right.boundingBox();
		if (!box) throw new Error("right bar bounding box missing");
		expect(box.width).toBeLessThanOrEqual(8);
	});

	test("empty bottom bar collapses to a sliver (≤8px tall)", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.waitForTimeout(500);
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

		// Simulate an in-flight drag by dispatching a window dragstart
		// event with the panel-id MIME in dataTransfer.types. We can't
		// use Playwright's dragTo because it completes the drop synchronously.
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

	/**
	 * Drag a resize handle by `deltaX` pixels using direct PointerEvent
	 * dispatch. Playwright's `page.mouse.*` does not drive handlers that
	 * rely on `setPointerCapture` + window pointer listeners, which is
	 * exactly what the resize handles use.
	 */
	async function dragResizeHandle(
		page: import("@playwright/test").Page,
		testId: string,
		deltaX: number,
		deltaY = 0,
	) {
		await page.evaluate(
			({ testId, dx, dy }) => {
				const handle = document.querySelector(
					`[data-testid="${testId}"]`,
				) as HTMLElement | null;
				if (!handle) throw new Error(`no handle ${testId}`);
				const rect = handle.getBoundingClientRect();
				const startX = rect.left + rect.width / 2;
				const startY = rect.top + rect.height / 2;
				const down = new PointerEvent("pointerdown", {
					bubbles: true,
					cancelable: true,
					clientX: startX,
					clientY: startY,
					pointerId: 1,
					pointerType: "mouse",
				});
				handle.dispatchEvent(down);
				const move = new PointerEvent("pointermove", {
					bubbles: true,
					cancelable: true,
					clientX: startX + dx,
					clientY: startY + dy,
					pointerId: 1,
					pointerType: "mouse",
				});
				window.dispatchEvent(move);
				const up = new PointerEvent("pointerup", {
					bubbles: true,
					cancelable: true,
					clientX: startX + dx,
					clientY: startY + dy,
					pointerId: 1,
					pointerType: "mouse",
				});
				window.dispatchEvent(up);
			},
			{ testId, dx: deltaX, dy: deltaY },
		);
	}

	test("sidebar resize: dragging the inner edge widens the zone", async ({
		page,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-geometries").click();
		const zone = page.getByTestId("sidebar-zone-left");
		const before = await zone.boundingBox();
		if (!before) throw new Error("zone bounding box missing");

		await dragResizeHandle(page, "sidebar-resize-left", 100);

		const after = await zone.boundingBox();
		if (!after) throw new Error("zone bounding box missing after");
		expect(after.width).toBeGreaterThan(before.width + 50);
	});

	test("sidebar resize clamps at SIDEBAR_MAX_PX", async ({ page }) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-geometries").click();
		const zone = page.getByTestId("sidebar-zone-left");

		await dragResizeHandle(page, "sidebar-resize-left", 2000);

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
		await dragResizeHandle(page, "sidebar-resize-left", 150);
		const widthAfterResize = (
			await page.getByTestId("sidebar-zone-left").boundingBox()
		)?.width;
		if (widthAfterResize === undefined) throw new Error("no width");
		expect(widthAfterResize).toBeGreaterThan(460); // sanity: drag actually moved it

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
		await dragResizeHandle(page, "sidebar-resize-left", 150);
		const widthAfterResize = (
			await page.getByTestId("sidebar-zone-left").boundingBox()
		)?.width;
		expect(widthAfterResize ?? 0).toBeGreaterThan(460); // sanity

		await page.reload();
		await page.getByTestId("activity-icon-geometries").click();
		const box = await page.getByTestId("sidebar-zone-left").boundingBox();
		expect(box?.width).toBeGreaterThanOrEqual(318);
		expect(box?.width).toBeLessThanOrEqual(322);
	});

	test("dragging widens empty bars to 56 px", async ({ page }) => {
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

		await page.getByTestId("rooms-search").locator("input").fill("zz-no-such-room-zz");
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
		const deleteItem = page.getByRole("menuitem", { name: /Delete/i });
		await expect(deleteItem).toHaveAttribute("aria-disabled", "true");
	});
});
