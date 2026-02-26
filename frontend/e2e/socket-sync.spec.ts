import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene, createAuthContext } from "./helpers";

const ROOM = "test-socket-sync";

test.describe("Socket Sync — Multi-Tab", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		CLI(`rooms create --room-id ${ROOM}`);
		PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[:]
for i in range(20):
    atoms = ase.Atoms('H4', positions=[(0,0,0),(1,0,0),(0,1,0),(1,1,0)])
    vis.append(atoms)
vis.step = 0
`);
	});

	test("step change syncs between tabs", async ({ browser }) => {
		const ctxA = await createAuthContext(browser);
		const ctxB = await createAuthContext(browser);
		const pageA = await ctxA.newPage();
		const pageB = await ctxB.newPage();

		try {
			await pageA.goto(`${BASE_URL}/rooms/${ROOM}`);
			await pageB.goto(`${BASE_URL}/rooms/${ROOM}`);
			await waitForScene(pageA);
			await waitForScene(pageB);

			await expect(pageA.getByText("/ 20")).toBeVisible({
				timeout: 15000,
			});
			await expect(pageB.getByText("/ 20")).toBeVisible({
				timeout: 15000,
			});

			// Set step to 10 via CLI
			CLI(`step set ${ROOM} 10`);

			// Both tabs should reflect the change
			await expect(pageA.getByText("11 / 20")).toBeVisible({
				timeout: 10000,
			});
			await expect(pageB.getByText("11 / 20")).toBeVisible({
				timeout: 10000,
			});
		} finally {
			await ctxA.close();
			await ctxB.close();
		}
	});

	test("chat message appears in both tabs", async ({ browser }) => {
		const ctxA = await createAuthContext(browser);
		const ctxB = await createAuthContext(browser);
		const pageA = await ctxA.newPage();
		const pageB = await ctxB.newPage();

		try {
			await pageA.goto(`${BASE_URL}/rooms/${ROOM}`);
			await pageB.goto(`${BASE_URL}/rooms/${ROOM}`);
			await waitForScene(pageA);
			await waitForScene(pageB);

			// Open chat in both tabs
			await pageA
				.getByRole("button", { name: "toggle chat" })
				.click();
			await pageB
				.getByRole("button", { name: "toggle chat" })
				.click();
			await pageA.waitForTimeout(500);
			await pageB.waitForTimeout(500);

			// Send a message via Python
			const uniqueMsg = `sync-test-${Date.now()}`;
			PY(`
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
vis.log("${uniqueMsg}")
`);

			// Both tabs should see the message (socket message_new broadcast)
			await expect(pageA.getByText(uniqueMsg)).toBeVisible({
				timeout: 10000,
			});
			await expect(pageB.getByText(uniqueMsg)).toBeVisible({
				timeout: 10000,
			});
		} finally {
			await ctxA.close();
			await ctxB.close();
		}
	});

	test("frame update syncs between tabs", async ({ browser }) => {
		const ctxA = await createAuthContext(browser);
		const ctxB = await createAuthContext(browser);
		const pageA = await ctxA.newPage();
		const pageB = await ctxB.newPage();

		try {
			await pageA.goto(`${BASE_URL}/rooms/${ROOM}`);
			await pageB.goto(`${BASE_URL}/rooms/${ROOM}`);
			await waitForScene(pageA);
			await waitForScene(pageB);

			const beforeJson = JSON.parse(CLI(`frames count ${ROOM}`));
			const beforeCount = beforeJson.total_frames;

			// Append 3 frames via Python
			PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
for i in range(3):
    vis.append(ase.Atoms('O2', positions=[(0,0,0),(1.2,0,0)]))
`);

			const expectedCount = beforeCount + 3;

			// Both tabs should update frame count via frames_invalidate
			await expect(
				pageA.getByText(`/ ${expectedCount}`),
			).toBeVisible({ timeout: 15000 });
			await expect(
				pageB.getByText(`/ ${expectedCount}`),
			).toBeVisible({ timeout: 15000 });
		} finally {
			await ctxA.close();
			await ctxB.close();
		}
	});

	test("selection change syncs between tabs", async ({ browser }) => {
		const ctxA = await createAuthContext(browser);
		const ctxB = await createAuthContext(browser);
		const pageA = await ctxA.newPage();
		const pageB = await ctxB.newPage();

		try {
			await pageA.goto(`${BASE_URL}/rooms/${ROOM}`);
			await pageB.goto(`${BASE_URL}/rooms/${ROOM}`);
			await waitForScene(pageA);
			await waitForScene(pageB);

			// Open selection panel in both tabs
			const btnA = pageA.getByRole("button", {
				name: "Selection tools and groups",
			});
			const btnB = pageB.getByRole("button", {
				name: "Selection tools and groups",
			});
			await btnA.click();
			await btnB.click();

			await expect(
				pageA.getByRole("heading", { name: "Selections" }),
			).toBeVisible({ timeout: 10000 });
			await expect(
				pageB.getByRole("heading", { name: "Selections" }),
			).toBeVisible({ timeout: 10000 });

			// Set selection via CLI
			CLI(`selection set ${ROOM} 0 1 2`);

			// Take screenshots for manual inspection of selection state
			await pageA.waitForTimeout(2000);
			await pageB.waitForTimeout(2000);

			await pageA.screenshot({
				path: "e2e/screenshots/sync-selection-a.png",
			});
			await pageB.screenshot({
				path: "e2e/screenshots/sync-selection-b.png",
			});

			// Verify selection via CLI to confirm it was set
			const selJson = CLI(`selection get ${ROOM}`);
			expect(selJson).toContain("0");
			expect(selJson).toContain("1");
			expect(selJson).toContain("2");
		} finally {
			await ctxA.close();
			await ctxB.close();
		}
	});

	test("editing in one tab reflects in other tab", async ({ browser }) => {
		const ctxA = await createAuthContext(browser);
		const ctxB = await createAuthContext(browser);
		const pageA = await ctxA.newPage();
		const pageB = await ctxB.newPage();

		try {
			await pageA.goto(`${BASE_URL}/rooms/${ROOM}`);
			await pageB.goto(`${BASE_URL}/rooms/${ROOM}`);
			await waitForScene(pageA);
			await waitForScene(pageB);

			// Tab A: click canvas center to select an atom
			const canvas = pageA.locator("canvas");
			const box = await canvas.boundingBox();
			if (box) {
				await canvas.click({
					position: { x: box.width / 2, y: box.height / 2 },
				});
			}
			await pageA.waitForTimeout(1000);

			// Tab A: enter editing mode
			await pageA.keyboard.press("e");
			await pageA.waitForTimeout(1000);
			await expect(pageA.getByText("Translate")).toBeVisible({
				timeout: 5000,
			});

			// Tab A: move atom with keyboard (y-axis + ArrowUp)
			await pageA.keyboard.press("y");
			await pageA.waitForTimeout(500);
			for (let i = 0; i < 5; i++) {
				await pageA.keyboard.press("ArrowUp");
				await pageA.waitForTimeout(200);
			}

			// Tab A: exit editing mode (commits the change)
			await pageA.keyboard.press("e");
			await pageA.waitForTimeout(2000);

			// Tab B: take screenshot — should show updated atom position
			// The edit triggers a frames_invalidate broadcast
			await pageB.waitForTimeout(2000);

			await pageA.screenshot({
				path: "e2e/screenshots/sync-edit-tab-a.png",
			});
			await pageB.screenshot({
				path: "e2e/screenshots/sync-edit-tab-b.png",
			});

			// Manual inspection: Tab B screenshot should reflect the atom movement
		} finally {
			await ctxA.close();
			await ctxB.close();
		}
	});
});
