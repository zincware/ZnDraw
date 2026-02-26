import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-geometry";
const ROOM_DEL = "test-geometry-del";

/**
 * Open the geometry panel by clicking the sidebar "Manage geometries" button.
 * The button is a ListItemButton inside a Tooltip with title "Manage geometries".
 */
async function openGeometryPanel(page: import("@playwright/test").Page) {
	const manageBtn = page.getByRole("button", { name: "Manage geometries" });
	await manageBtn.click();
	// Wait for the DataGrid to appear inside the panel
	await page.waitForSelector('[role="grid"]', {
		state: "visible",
		timeout: 10000,
	});
}

function setupRoom(room: string) {
	CLI(`rooms create --room-id ${room}`);
	PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms('H4', positions=[(0,0,0),(2,0,0),(0,2,0),(2,2,0)]))
from zndraw.geometries import Box, Sphere, Curve, Arrow, Floor
vis.geometries['floor'] = Floor(active=True, height=-2.0, color='#808080')
vis.geometries['box'] = Box(position=[(0,2,0)], size=[(4,4,4)], color=['#e74c3c'], material='MeshToonMaterial')
vis.geometries['sphere'] = Sphere(position=[(8,2,0)], radius=[2.0], color=['#3498db'], material='MeshPhysicalMaterial_glass')
vis.geometries['curve'] = Curve(position=[(-6,0,-6),(-3,4,-3),(0,0,0),(3,4,3),(6,0,6)], color='#2ecc71')
vis.geometries['arrow'] = Arrow(position=[(12,0,0)], direction=[(0,5,0)], color=['#f39c12'])
`);
}

test.describe("Geometry & Drawing", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		setupRoom(ROOM);
	});

	test("geometries render and panel shows all entries", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		await page.screenshot({ path: "e2e/screenshots/geometry-scene.png" });

		await openGeometryPanel(page);
		await page.waitForTimeout(500);

		await page.screenshot({ path: "e2e/screenshots/geometry-panel.png" });

		const grid = page.locator('[role="grid"]');
		await expect(grid).toBeVisible();

		for (const name of ["box", "sphere", "curve", "arrow", "floor"]) {
			await expect(grid.getByText(name, { exact: true }).first()).toBeVisible({
				timeout: 5000,
			});
		}
	});

	test("editing mode toggles with e key", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Click the canvas to ensure it has focus before pressing 'e'
		const canvas = page.locator("canvas");
		await canvas.click();
		await page.waitForTimeout(500);

		// Press 'e' to enter editing mode
		await page.keyboard.press("e");
		await page.waitForTimeout(1000);

		const editChip = page.getByText("Translate");
		await expect(editChip).toBeVisible({ timeout: 10000 });

		await page.screenshot({
			path: "e2e/screenshots/geometry-edit-mode.png",
		});

		// Click canvas again to ensure focus before pressing 'e'
		await canvas.click();
		await page.waitForTimeout(500);

		// Press 'e' again to exit editing mode
		await page.keyboard.press("e");

		await expect(editChip).not.toBeVisible({ timeout: 5000 });

		await page.screenshot({
			path: "e2e/screenshots/geometry-exit-edit-mode.png",
		});
	});
});

test.describe("Geometry Deletion", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		setupRoom(ROOM_DEL);
	});

	test("geometry deletion via CLI reflects in UI", async ({ page }) => {
		test.setTimeout(120_000);
		await page.goto(`${BASE_URL}/rooms/${ROOM_DEL}`);
		await waitForScene(page);

		await openGeometryPanel(page);
		const grid = page.locator('[role="grid"]');
		await expect(grid.getByText("arrow", { exact: true }).first()).toBeVisible({
			timeout: 5000,
		});

		await page.screenshot({
			path: "e2e/screenshots/geometry-before-delete.png",
		});

		// Delete the 'arrow' geometry via CLI
		CLI(`geometries delete ${ROOM_DEL} arrow`);
		await page.waitForTimeout(1000);

		// Reload and verify the arrow is gone from the panel
		await page.reload();
		await waitForScene(page);
		await openGeometryPanel(page);

		const gridAfter = page.locator('[role="grid"]');
		await expect(gridAfter).toBeVisible();

		await expect(gridAfter.getByText("arrow", { exact: true })).not.toBeVisible(
			{ timeout: 5000 },
		);

		for (const name of ["box", "sphere", "curve", "floor"]) {
			await expect(
				gridAfter.getByText(name, { exact: true }).first(),
			).toBeVisible({ timeout: 5000 });
		}

		await page.screenshot({
			path: "e2e/screenshots/geometry-after-delete.png",
		});
	});
});
