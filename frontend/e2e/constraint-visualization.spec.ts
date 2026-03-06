import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-constraints";

/**
 * Open the geometry panel by clicking the sidebar "Manage geometries" button.
 */
async function openGeometryPanel(page: import("@playwright/test").Page) {
	const manageBtn = page.getByRole("button", { name: "Manage geometries" });
	await manageBtn.click();
	await page.waitForSelector('[role="grid"]', {
		state: "visible",
		timeout: 10000,
	});
}

function setupConstraintRoom() {
	CLI(`rooms create --room-id ${ROOM}`);
	PY(`
from zndraw import ZnDraw
from molify import smiles2conformers
from ase.constraints import FixAtoms

vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[:]

# Butyric acid: constrain the carbon chain
atoms = smiles2conformers('CCCC(=O)O', numConfs=1)[0]
carbon_indices = [i for i, s in enumerate(atoms.symbols) if s == 'C']
atoms.set_constraint(FixAtoms(indices=carbon_indices))
vis.append(atoms)
`);
}

test.describe("Constraint Visualization", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		setupConstraintRoom();
	});

	test("constraints-fixed-atoms geometry appears in panel", async ({
		page,
	}) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		await openGeometryPanel(page);
		const grid = page.locator('[role="grid"]');
		await expect(
			grid.getByText("constraints-fixed-atoms", { exact: true }).first(),
		).toBeVisible({ timeout: 5000 });

		await page.screenshot({
			path: "e2e/screenshots/constraint-geometry-panel.png",
		});
	});

	test("constraint geometry renders in scene", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Wait for the scene to fully render with constraints
		await page.waitForTimeout(3000);

		await page.screenshot({
			path: "e2e/screenshots/constraint-visualization.png",
		});
	});
});
