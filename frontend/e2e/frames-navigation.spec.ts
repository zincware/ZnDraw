import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-frames";

/**
 * Wait for the frame counter to show the expected total frame count,
 * confirming that all frames have been indexed by the frontend.
 */
async function waitForFrameCount(
	page: import("@playwright/test").Page,
	total: number,
) {
	await expect(page.getByText(`/ ${total}`)).toBeVisible({ timeout: 15000 });
}

test.describe("Frames & Navigation", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		CLI(`rooms create --room-id ${ROOM}`);
		// Clear any existing data from previous runs and set up fresh
		PY(`
from zndraw import ZnDraw
import ase
import numpy as np
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
# Clear existing frames
del vis[:]
for i in range(20):
    atoms = ase.Atoms('OH2', positions=[(0,0,0),(0.96,0,0),(0.59,0.77,0)])
    atoms.positions += np.random.default_rng(i).normal(0, 0.05, (3,3))
    atoms.info['energy'] = -10.0 + np.sin(i/3.0)
    vis.append(atoms)
vis.bookmarks[5] = 'Interesting'
vis.bookmarks[15] = 'Peak energy'
vis.selection = [0, 1]
vis.selection_groups['oxygen'] = {'particles': [0]}
vis.selection_groups['hydrogen'] = {'particles': [1, 2]}
vis.step = 0
`);
	});

	test("room loads with correct frame count", async ({ page }) => {
		// Get the actual frame count from CLI first
		const countJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const totalFrames = countJson.total_frames;

		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);
		await waitForFrameCount(page, totalFrames);

		// The frame counter shows "1 / N" (1-indexed display, starting at frame 0)
		const frameCounter = page.getByText(`1 / ${totalFrames}`);
		await expect(frameCounter).toBeVisible();

		// The slider should be present and at value 0
		const slider = page.getByRole("slider", { name: "Frame Progress" });
		await expect(slider).toBeVisible();

		await page.screenshot({ path: "e2e/screenshots/frames-initial.png" });
	});

	test("step navigation updates view", async ({ page }) => {
		const countJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const totalFrames = countJson.total_frames;

		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);
		await waitForFrameCount(page, totalFrames);

		// Set step to 10 via CLI
		CLI(`step set ${ROOM} 10`);

		// Wait for the UI to reflect the change via socket.io sync.
		// The frame counter should update to "11 / N" (1-indexed).
		await expect(page.getByText(`11 / ${totalFrames}`)).toBeVisible({
			timeout: 10000,
		});

		// The slider value should be 10
		const slider = page.getByRole("slider", { name: "Frame Progress" });
		await expect(slider).toHaveAttribute("aria-valuenow", "10");

		await page.screenshot({ path: "e2e/screenshots/frames-step10.png" });

		// Navigate to a different step
		CLI(`step set ${ROOM} 0`);
		await expect(page.getByText(`1 / ${totalFrames}`)).toBeVisible({
			timeout: 10000,
		});

		await page.screenshot({ path: "e2e/screenshots/frames-step0.png" });
	});

	test("bookmarks visible on timeline", async ({ page }) => {
		const countJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const totalFrames = countJson.total_frames;

		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);
		await waitForFrameCount(page, totalFrames);

		// Bookmark buttons should be visible on the timeline
		const interestingBookmark = page.getByRole("button", {
			name: /Bookmark: Interesting/,
		});
		const peakBookmark = page.getByRole("button", {
			name: /Bookmark: Peak energy/,
		});

		await expect(interestingBookmark).toBeVisible({ timeout: 10000 });
		await expect(peakBookmark).toBeVisible({ timeout: 10000 });

		// Verify the bookmark labels include the frame numbers
		await expect(interestingBookmark).toHaveAttribute(
			"aria-label",
			"Bookmark: Interesting (Frame 5)",
		);
		await expect(peakBookmark).toHaveAttribute(
			"aria-label",
			"Bookmark: Peak energy (Frame 15)",
		);

		// Navigate to the bookmark frame and verify the frame counter
		CLI(`step set ${ROOM} 5`);
		await expect(page.getByText(`6 / ${totalFrames}`)).toBeVisible({
			timeout: 10000,
		});

		await page.screenshot({ path: "e2e/screenshots/frames-bookmark.png" });

		// Verify bookmarks are still visible from the CLI
		const bookmarksOutput = CLI(`bookmarks list ${ROOM}`);
		expect(bookmarksOutput).toContain('"Interesting"');
		expect(bookmarksOutput).toContain('"Peak energy"');
	});

	test("selection panel opens and shows groups", async ({ page }) => {
		const countJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const totalFrames = countJson.total_frames;

		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);
		await waitForFrameCount(page, totalFrames);

		// Click the "Selection tools and groups" button in the sidebar
		const selectionBtn = page.getByRole("button", {
			name: "Selection tools and groups",
		});
		await expect(selectionBtn).toBeVisible();
		await selectionBtn.click();

		// The panel should open with the "Selections" heading
		const heading = page.getByRole("heading", { name: "Selections" });
		await expect(heading).toBeVisible({ timeout: 10000 });

		// "Selection Tools" accordion should be present
		const selectionToolsHeading = page.getByRole("heading", {
			name: "Selection Tools",
		});
		await expect(selectionToolsHeading.first()).toBeVisible();

		// "Selection Groups" accordion should be present
		const selectionGroupsHeading = page.getByRole("heading", {
			name: "Selection Groups",
		});
		await expect(selectionGroupsHeading.first()).toBeVisible();

		// The DataGrid should contain the named groups "hydrogen" and "oxygen"
		const grid = page.locator('[role="grid"]');
		await expect(grid).toBeVisible({ timeout: 10000 });

		// Verify group names appear in the grid
		await expect(grid.getByText("hydrogen")).toBeVisible({ timeout: 5000 });
		await expect(grid.getByText("oxygen")).toBeVisible({ timeout: 5000 });

		// The "particles" column header should be visible
		const particlesHeader = grid.getByRole("columnheader", {
			name: "particles",
		});
		await expect(particlesHeader).toBeVisible();

		await page.screenshot({
			path: "e2e/screenshots/frames-selection-panel.png",
		});

		// Close the panel by clicking the button again
		await selectionBtn.click();
		await expect(heading).not.toBeVisible({ timeout: 5000 });

		await page.screenshot({
			path: "e2e/screenshots/frames-selection-panel-closed.png",
		});
	});
});
