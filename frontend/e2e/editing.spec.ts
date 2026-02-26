import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-editing";
const ROOM_LOCK = "test-editing-lock";

function setupRoom(room: string) {
	CLI(`rooms create --room-id ${room}`);
	PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms('H4', positions=[(0,0,0),(3,0,0),(0,3,0),(3,3,0)]))
vis.step = 0
`);
}

test.describe("Editing Mode", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		setupRoom(ROOM);
	});

	test("select atom and enter editing mode", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		await page.screenshot({
			path: "e2e/screenshots/editing-before-select.png",
		});

		// Click near the center of the canvas to attempt atom selection
		const canvas = page.locator("canvas");
		const box = await canvas.boundingBox();
		if (box) {
			await canvas.click({
				position: { x: box.width / 2, y: box.height / 2 },
			});
		}
		await page.waitForTimeout(1000);

		await page.screenshot({
			path: "e2e/screenshots/editing-after-select.png",
		});

		// Press 'e' to enter editing mode
		await page.keyboard.press("e");
		await page.waitForTimeout(1000);

		// Verify editing mode: Translate chip visible
		const editChip = page.getByText("Translate");
		await expect(editChip).toBeVisible({ timeout: 5000 });

		await page.screenshot({
			path: "e2e/screenshots/editing-mode-entered.png",
		});

		// Exit editing mode
		await page.keyboard.press("e");

		await expect(editChip).not.toBeVisible({ timeout: 5000 });
	});

	test("keyboard translate moves atom", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Click canvas center to select an atom
		const canvas = page.locator("canvas");
		const box = await canvas.boundingBox();
		if (box) {
			await canvas.click({
				position: { x: box.width / 2, y: box.height / 2 },
			});
		}
		await page.waitForTimeout(1000);

		// Enter editing mode
		await page.keyboard.press("e");
		await page.waitForTimeout(1000);
		await expect(page.getByText("Translate")).toBeVisible({
			timeout: 5000,
		});

		await page.screenshot({
			path: "e2e/screenshots/editing-translate-before.png",
		});

		// Press 'y' to select Y axis, then ArrowUp to translate
		await page.keyboard.press("y");
		await page.waitForTimeout(500);
		for (let i = 0; i < 5; i++) {
			await page.keyboard.press("ArrowUp");
			await page.waitForTimeout(200);
		}

		await page.waitForTimeout(1000);

		await page.screenshot({
			path: "e2e/screenshots/editing-translate-after.png",
		});

		// Exit editing mode
		await page.keyboard.press("e");
		await page.waitForTimeout(500);

		// Manual inspection: compare before/after screenshots
		// The atom should have moved upward in the viewport
	});
});

test.describe("Edit Lock — Two Tabs Same User", () => {
	test.beforeAll(() => {
		setupRoom(ROOM_LOCK);
	});

	test("second tab of same user cannot acquire lock held by first tab", async ({
		context,
	}) => {
		const tab1 = await context.newPage();
		const tab2 = await context.newPage();

		await tab1.goto(`${BASE_URL}/rooms/${ROOM_LOCK}`);
		await waitForScene(tab1);
		await tab2.goto(`${BASE_URL}/rooms/${ROOM_LOCK}`);
		await waitForScene(tab2);

		// Tab 1: acquire the edit lock
		const canvas1 = tab1.locator("canvas");
		await canvas1.click();
		await tab1.waitForTimeout(500);
		await tab1.keyboard.press("e");
		await expect(tab1.getByText("Translate")).toBeVisible({ timeout: 10000 });

		await tab1.screenshot({
			path: "e2e/screenshots/lock-tab1-editing.png",
		});

		// Tab 2: try to enter editing mode — should be blocked (423)
		await tab2.bringToFront();
		const canvas2 = tab2.locator("canvas");
		await canvas2.click();
		await tab2.waitForTimeout(500);
		await tab2.keyboard.press("e");
		await tab2.waitForTimeout(2000);

		// Tab 2 should NOT be in editing mode
		await expect(tab2.getByText("Translate")).not.toBeVisible({
			timeout: 5000,
		});

		await tab2.screenshot({
			path: "e2e/screenshots/lock-tab2-blocked.png",
		});

		// Tab 1: exit editing mode to release the lock
		await tab1.bringToFront();
		await canvas1.click();
		await tab1.waitForTimeout(500);
		await tab1.keyboard.press("e");
		await expect(tab1.getByText("Translate")).not.toBeVisible({
			timeout: 5000,
		});

		// Wait for lock release to propagate via socket
		await tab2.waitForTimeout(2000);

		// Tab 2: now should be able to enter editing mode
		// Retry up to 3 times — headed Playwright may need re-focus between tabs
		await tab2.bringToFront();
		await canvas2.click();
		await tab2.waitForTimeout(500);

		for (let attempt = 0; attempt < 3; attempt++) {
			await tab2.keyboard.press("e");
			try {
				await expect(tab2.getByText("Translate")).toBeVisible({
					timeout: 5000,
				});
				break;
			} catch {
				if (attempt === 2)
					throw new Error("Tab2 failed to enter editing mode after 3 attempts");
				await canvas2.click();
				await tab2.waitForTimeout(1000);
			}
		}

		await tab2.screenshot({
			path: "e2e/screenshots/lock-tab2-acquired.png",
		});

		// Cleanup
		await canvas2.click();
		await tab2.waitForTimeout(500);
		await tab2.keyboard.press("e");
		await tab2.waitForTimeout(500);

		await tab1.close();
		await tab2.close();
	});
});
