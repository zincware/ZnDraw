import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-chat-features";

test.describe("Chat Features", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		CLI(`rooms create --room-id ${ROOM}`);
		PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[:]
for i in range(10):
    vis.append(ase.Atoms('H2', positions=[(0,0,0),(0.74,0,0)]))
vis.step = 0
`);
	});

	test("SMILES code block renders molecule preview", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Send a message with a SMILES code block via Python
		PY(`
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
vis.log("\`\`\`smiles\\nCCO\\n\`\`\`")
`);

		// Open the chat panel
		await page.getByRole("button", { name: "toggle chat" }).click();
		await page.waitForTimeout(500);

		await expect(
			page.getByRole("heading", { name: "Chat" }),
		).toBeVisible();

		// MoleculePreview renders as <img alt="Molecule structure">
		const moleculeImg = page.getByAltText("Molecule structure");
		await expect(moleculeImg.first()).toBeVisible({ timeout: 15000 });

		await page.screenshot({
			path: "e2e/screenshots/chat-smiles-preview.png",
		});
	});

	test("@frame reference renders clickable chip that jumps to frame", async ({
		page,
	}) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Verify starting state
		await expect(page.getByText("/ 10")).toBeVisible({ timeout: 15000 });

		// Send a message with @5 frame reference
		// The remark plugin parses @5 (1-based user input) → frame index 4 (0-based).
		// FrameReference displays frame+1 = 5 as the chip label.
		PY(`
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
vis.log("See @5 for details")
`);

		// Open chat
		await page.getByRole("button", { name: "toggle chat" }).click();
		await page.waitForTimeout(500);

		// FrameReference renders a clickable MUI Chip with label "5"
		const frameChip = page
			.locator(".MuiChip-root")
			.filter({ hasText: "5" });
		await expect(frameChip.first()).toBeVisible({ timeout: 10000 });

		await page.screenshot({
			path: "e2e/screenshots/chat-frame-reference.png",
		});

		// Click the chip to jump to frame 4 (0-based) → display "5 / 10"
		await frameChip.first().click();

		// Frame counter should update to "5 / 10"
		await expect(page.getByText("5 / 10")).toBeVisible({ timeout: 10000 });

		await page.screenshot({
			path: "e2e/screenshots/chat-frame-reference-jumped.png",
		});
	});
});
