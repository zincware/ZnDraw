import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-frame-invalidation";

test.describe("Frame Invalidation", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		CLI(`rooms create --room-id ${ROOM}`);
		PY(`
from zndraw import ZnDraw
import ase
from ase.build import molecule
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[:]
# Start with ethanol (9 atoms: C2H5OH)
vis.append(molecule('CH3CH2OH'))
`);
	});

	test("in-place frame update changes displayed atoms", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Verify initial state: 1 frame
		await expect(page.getByText("/ 1")).toBeVisible({ timeout: 15000 });

		await page.screenshot({
			path: "e2e/screenshots/invalidation-before-update.png",
		});

		// Replace frame 0 with benzene (12 atoms) via Python
		PY(`
from zndraw import ZnDraw
from ase.build import molecule
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
vis.atoms = molecule('C6H6')
`);

		// The frontend should refetch the frame via frames_invalidate socket event.
		// Wait for a visual change — take screenshot for manual inspection.
		await page.waitForTimeout(3000);

		await page.screenshot({
			path: "e2e/screenshots/invalidation-after-update.png",
		});

		// Verify frame count is still 1 (in-place update, not append)
		await expect(page.getByText("/ 1")).toBeVisible();
	});

	test("appending frames updates frame count", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		const beforeJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const beforeCount = beforeJson.total_frames;

		// Append 5 more frames
		PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
for i in range(5):
    vis.append(ase.Atoms('H2', positions=[(0,0,0),(0.74,0,0)]))
`);

		// Frontend should update frame count via frames_invalidate
		const expectedTotal = beforeCount + 5;
		await expect(page.getByText(`/ ${expectedTotal}`)).toBeVisible({
			timeout: 15000,
		});
	});

	test("deleting frames updates frame count", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		const beforeJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const beforeCount = beforeJson.total_frames;

		// Delete first frame
		PY(`
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[0]
`);

		// Frontend should update frame count
		const expectedTotal = beforeCount - 1;
		await expect(page.getByText(`/ ${expectedTotal}`)).toBeVisible({
			timeout: 15000,
		});
	});
});
