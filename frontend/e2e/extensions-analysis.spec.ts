import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-extensions";

/**
 * Open the modifier panel by clicking the sidebar "Modifier tools" button
 * and wait for the panel heading to appear.
 */
async function openModifierPanel(page: import("@playwright/test").Page) {
	await page.getByRole("button", { name: "Modifier tools" }).click();
	await page
		.getByRole("heading", { name: "modifiers" })
		.waitFor({ state: "visible", timeout: 10000 });
}

/**
 * Open the analysis panel by clicking the sidebar "Analysis tools" button
 * and wait for the panel heading to appear.
 */
async function openAnalysisPanel(page: import("@playwright/test").Page) {
	await page.getByRole("button", { name: "Analysis tools" }).click();
	await page
		.getByRole("heading", { name: "analysis" })
		.waitFor({ state: "visible", timeout: 10000 });
}

/**
 * Select an option from a MUI Select dropdown identified by its label.
 * MUI Select renders as a combobox — click to open the listbox, then
 * click the matching option.
 */
async function selectFromDropdown(
	page: import("@playwright/test").Page,
	label: string,
	optionText: string,
) {
	await page.getByRole("combobox", { name: label }).click();
	// Wait for MUI menu animation to settle
	await page.waitForTimeout(500);
	// MUI Select opens a listbox with role="option" items
	await page
		.getByRole("option", { name: optionText })
		.click({ timeout: 15000 });
}

test.describe("Extensions & Analysis", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		CLI(`rooms create --room-id ${ROOM}`);
		PY(`
from zndraw import ZnDraw
import ase
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[:]
for i in range(30):
    atoms = ase.Atoms('H4', positions=[(0,0,0),(1,0,0),(0,1,0),(1,1,0)])
    calc = SinglePointCalculator(atoms, energy=-10 + np.sin(i/5.0))
    atoms.calc = calc
    vis.append(atoms)
`);
	});

	test("AddFromSMILES modifier adds a frame", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Capture initial frame count via CLI
		const beforeJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const beforeCount = beforeJson.total_frames;

		await page.screenshot({ path: "e2e/screenshots/ext-initial.png" });

		// Open the modifiers panel
		await openModifierPanel(page);
		await page.waitForTimeout(500);

		await page.screenshot({
			path: "e2e/screenshots/ext-modifier-panel.png",
		});

		// Select AddFromSMILES from the modifiers Method dropdown
		await selectFromDropdown(page, "modifiers Method", "AddFromSMILES");
		await page.waitForTimeout(500);

		await page.screenshot({ path: "e2e/screenshots/ext-add-smiles.png" });

		// Fill the SMILES input with "CCO" (ethanol)
		const smilesInput = page.getByPlaceholder("Enter SMILES notation");
		await smilesInput.waitFor({ state: "visible", timeout: 10000 });
		await smilesInput.fill("CCO");
		await page.waitForTimeout(500);

		await page.screenshot({
			path: "e2e/screenshots/ext-smiles-filled.png",
		});

		// Click "Run Extension"
		await page.getByRole("button", { name: "Run Extension" }).click();

		// Wait for the internal executor to process the task
		await page.waitForTimeout(5000);

		await page.screenshot({
			path: "e2e/screenshots/ext-smiles-result.png",
		});

		// Verify the frame count increased via CLI
		const afterJson = JSON.parse(CLI(`frames count ${ROOM}`));
		const afterCount = afterJson.total_frames;
		expect(afterCount).toBeGreaterThan(beforeCount);
	});

	test("Properties1D analysis creates a figure", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Open the analysis panel
		await openAnalysisPanel(page);
		await page.waitForTimeout(500);

		await page.screenshot({
			path: "e2e/screenshots/ext-analysis-panel.png",
		});

		// Select Properties1D from the analysis Method dropdown
		await selectFromDropdown(page, "analysis Method", "Properties1D");
		// Wait for dynamic properties to load from server
		await page.waitForTimeout(2000);

		await page.screenshot({
			path: "e2e/screenshots/ext-analysis-selected.png",
		});

		// The Value field is an MUI Autocomplete. Click to open, then select.
		const valueInput = page.getByRole("combobox", { name: "Value" });
		await valueInput.waitFor({ state: "visible", timeout: 10000 });
		await valueInput.click();
		// Wait for async options to load
		await page.waitForTimeout(2000);

		await page.screenshot({
			path: "e2e/screenshots/ext-analysis-dropdown.png",
		});

		// Try to find an option containing "energy" — the key name may vary
		const energyOption = page
			.getByRole("option")
			.filter({ hasText: /energy/i });
		const optionCount = await energyOption.count();

		if (optionCount > 0) {
			await energyOption.first().click();
		} else {
			// Fall back: select the first available option
			const firstOption = page.getByRole("option").first();
			await firstOption.click({ timeout: 5000 });
		}
		await page.waitForTimeout(500);

		await page.screenshot({
			path: "e2e/screenshots/ext-analysis-configured.png",
		});

		// Click "Run Extension"
		await page.getByRole("button", { name: "Run Extension" }).click();

		// Wait for the internal executor to process the task
		await page.waitForTimeout(5000);

		await page.screenshot({
			path: "e2e/screenshots/ext-analysis-result.png",
		});

		// Verify a figure was created via CLI — items are strings (figure keys)
		const figuresJson = JSON.parse(CLI(`figures list ${ROOM}`));
		const figureKeys: string[] = figuresJson.items;
		expect(figureKeys.length).toBeGreaterThan(0);
	});
});
