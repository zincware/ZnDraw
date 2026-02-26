import { test } from "@playwright/test";
import * as path from "path";
import { fileURLToPath } from "url";
import { BASE_URL, PY, waitForScene, spawnPY, waitForBgReady } from "./helpers";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const SCREENSHOTS_DIR = path.resolve(
	__dirname,
	"../../docs/source/_static/screenshots",
);

for (const theme of ["light", "dark"] as const) {
	const outDir = path.join(SCREENSHOTS_DIR, `${theme}mode`);
	const prefix = `docs-${theme}`;

	test.describe(`Documentation Screenshots (${theme})`, () => {
		test.describe.configure({ mode: "serial" });

		test.beforeAll(() => {
			// Compute bmim_bf4 into a theme-specific source room
			PY(`
import molify
from zndraw import ZnDraw
bf4 = molify.smiles2conformers("[B-](F)(F)(F)F", numConfs=100)
bmim = molify.smiles2conformers("CCCCN1C=C[N+](=C1)C", numConfs=100)
mols = [molify.pack([[bf4[i]], [bmim[i]]], counts=[1, 1], density=1200) for i in range(32)]
bmim_bf4 = molify.pack([mols], counts=[32], density=1200)
vis = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
del vis[:]
vis.append(bmim_bf4)
`);
			// Compute bmim_bf4_e_f (100 frames with energies/forces)
			PY(`
import numpy as np
import molify
import ase
from ase.calculators.singlepoint import SinglePointCalculator
from zndraw import ZnDraw

bf4 = molify.smiles2conformers("[B-](F)(F)(F)F", numConfs=100)
bmim = molify.smiles2conformers("CCCCN1C=C[N+](=C1)C", numConfs=100)
mols = [molify.pack([[bf4[i]], [bmim[i]]], counts=[1, 1], density=1200) for i in range(32)]
bmim_bf4 = molify.pack([mols], counts=[32], density=1200)

rng = np.random.default_rng(42)
n_frames = 100
t = np.linspace(0, 4 * np.pi, n_frames)
base_energy = 10 * np.sin(t)
noise_energy = rng.normal(0, 1.5, n_frames)
noise_dft = rng.normal(0, 1.5, n_frames)

frames = []
for i in range(n_frames):
    atoms = bmim_bf4.copy()
    atoms.calc = SinglePointCalculator(
        atoms,
        energy=base_energy[i] + noise_energy[i],
        forces=rng.uniform(-0.1, 0.1, size=(len(atoms), 3)),
    )
    atoms.calc.results["dft_energy"] = base_energy[i] + noise_dft[i]
    atoms.calc.results["dft_forces"] = rng.uniform(-0.1, 0.1, size=(len(atoms), 3))
    frames.append(atoms)

vis = ZnDraw(url='${BASE_URL}', room='${prefix}-ef-source')
del vis[:]
vis.extend(frames)
`);
		});

		test.beforeEach(async ({ page }) => {
			await page.addInitScript((mode) => {
				localStorage.setItem("mui-mode", mode);
			}, theme);
		});

		// Helper: clone bmim_bf4 from source room to a new room
		function cloneBmim(room: string) {
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
dst = ZnDraw(url='${BASE_URL}', room='${room}')
del dst[:]
dst.append(src[0])
`);
		}

		// Helper: clone bmim_bf4 and extend N copies
		function cloneBmimExtend(room: string, n: number) {
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
dst = ZnDraw(url='${BASE_URL}', room='${room}')
del dst[:]
dst.extend([src[0] for _ in range(${n})])
`);
		}

		// Helper: clone bmim_bf4_e_f frames to a new room
		function cloneEF(room: string, n?: number) {
			const sliceExpr = n !== undefined ? `src[:${n}]` : "src[:]";
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-ef-source')
dst = ZnDraw(url='${BASE_URL}', room='${room}')
del dst[:]
dst.extend(${sliceExpr})
`);
		}

		test("overview", async ({ page }) => {
			const room = `${prefix}-overview`;
			cloneBmim(room);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await waitForScene(page);
			await page.waitForTimeout(2000);
			await page.screenshot({ path: path.join(outDir, "overview.png") });
		});

		test("selection", async ({ page }) => {
			const room = `${prefix}-selection`;
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
dst = ZnDraw(url='${BASE_URL}', room='${room}')
del dst[:]
dst.append(src[0])
dst.selection = list(range(100))
dst.selection_groups["bf4"] = {"particles": list(range(32 * 5))}
dst.selection_groups["bmim"] = {"particles": list(range(32 * 25))}
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page
				.getByRole("button", { name: "Selection tools and groups" })
				.click();
			await page.mouse.click(0, 0);
			await page.waitForTimeout(2500);
			await page.screenshot({ path: path.join(outDir, "selection.png") });
		});

		test("drawing_mode", async ({ page }) => {
			const room = `${prefix}-drawing-mode`;
			PY(`
import ase
from zndraw import ZnDraw
from zndraw.geometries import Floor, Box, Curve
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms())
vis.geometries["floor"] = Floor(active=True, height=-5.0, color="#A0A0A0")
vis.geometries["box"] = Box(position=[(0, 0, 0)], size=[(10, 10, 10)])
vis.geometries["curve"] = Curve(
    position=[
        (-4.789177573713647, 5.474554893697157, 5.646802136489912),
        (-4.272226349849859, -4.421872431387516, 5.1448314636945724),
        (5.55918672834876, -4.537280536828369, 5.1448314636945724),
        (5.58483963937306, 5.279062292631962, 5.1448314636945724),
        (5.651912331581116, 5.290737777511481, -4.622419531622475),
        (5.651912331581116, -4.54817732764105, -4.742356696752667),
        (12.155146865890176, 1.703348269955838, 1.9575992621569651),
    ]
)
vis.selections["curve"] = [0]
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(2500);
			await page.keyboard.press("e");
			await page.screenshot({ path: path.join(outDir, "drawing_mode.png") });
		});

		test("geometries", async ({ page }) => {
			const room = `${prefix}-geometries`;
			PY(`
import ase
from zndraw import ZnDraw
from zndraw.geometries import Floor, Box, Sphere, Curve, Arrow
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms())
for key in list(vis.geometries.keys()):
    del vis.geometries[key]
vis.geometries["floor"] = Floor(active=True, height=-2.0, color="#808080")
vis.geometries["box"] = Box(
    position=[(0, 2, 0)], size=[(4, 4, 4)], color=["#e74c3c"],
    material="MeshToonMaterial",
)
vis.geometries["sphere"] = Sphere(
    position=[(8, 2, 0)], radius=[2.0], color=["#3498db"],
    material="MeshPhysicalMaterial_glass",
)
vis.geometries["curve"] = Curve(
    position=[(-6, 0, -6), (-3, 4, -3), (0, 0, 0), (3, 4, 3), (6, 0, 6)],
    color="#2ecc71",
)
vis.geometries["arrow"] = Arrow(
    position=[(12, 0, 0)], direction=[(0, 5, 0)], color=["#f39c12"],
    material="MeshPhysicalMaterial_shiny",
)
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(2000);
			await page.screenshot({ path: path.join(outDir, "geometries.png") });
		});

		test("geometry_viewer", async ({ page }) => {
			const room = `${prefix}-geometry-viewer`;
			PY(`
from zndraw import ZnDraw
from zndraw.geometries import Floor, Curve
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(src[0])
vis.geometries["floor"] = Floor(active=True, height=-5.0, color="#808080")
vis.geometries["my_curve"] = Curve(
    position=[(-10, 0, -10), (0, 10, 0), (10, 0, 10)],
    color="#2ecc71",
)
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(2500);
			await page.getByRole("button", { name: "Manage geometries" }).click();
			await page.screenshot({ path: path.join(outDir, "geometry_viewer.png") });
		});

		test("bookmarks", async ({ page }) => {
			const room = `${prefix}-bookmarks`;
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.extend([src[0] for _ in range(20)])
vis.bookmarks[5] = "Frame 5"
vis.bookmarks[17] = "Frame 17"
vis.step = 5
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page
				.getByRole("button", { name: "Bookmark: Frame 5 (Frame 5)" })
				.hover();
			await page.waitForTimeout(400);
			await page.screenshot({ path: path.join(outDir, "bookmarks.png") });
		});

		test("timeline", async ({ page }) => {
			const room = `${prefix}-timeline`;
			cloneBmimExtend(room, 50);
			PY(`
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${room}')
vis.step = 25
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1000);
			await page.screenshot({ path: path.join(outDir, "timeline.png") });
		});

		test("analysis_1d", async ({ page }) => {
			const room = `${prefix}-analysis-1d`;
			cloneEF(room);
			PY(`
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${room}')
vis.step = 16
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1000);
			await page.getByRole("button", { name: "Analysis tools" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("analysis Method").click();
			await page.getByRole("option", { name: "Properties1D" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("Value").click();
			await page.getByRole("option", { name: "calc.energy" }).click();
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Run Extension" }).click();
			await page.waitForTimeout(2000);
			await page.screenshot({ path: path.join(outDir, "analysis_1d.png") });
		});

		test("analysis_2d", async ({ page }) => {
			const room = `${prefix}-analysis-2d`;
			cloneEF(room);
			PY(`
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${room}')
vis.step = 16
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1000);
			// Run Properties1D first
			await page.getByRole("button", { name: "Analysis tools" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("analysis Method").click();
			await page.getByRole("option", { name: "Properties1D" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("Value").click();
			await page.getByRole("option", { name: "calc.energy" }).click();
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Run Extension" }).click();
			await page.waitForTimeout(2000);
			// Now run Properties2D
			await page.getByLabel("analysis Method").click();
			await page.getByRole("option", { name: "Properties2D" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("X data").click();
			await page.getByRole("option", { name: "calc.energy" }).click();
			await page.waitForTimeout(300);
			await page.getByLabel("Y data").click();
			await page.getByRole("option", { name: "calc.dft_energy" }).click();
			await page.waitForTimeout(300);
			await page.getByLabel("Color").click();
			await page.getByRole("option", { name: "calc.energy" }).click();
			await page.waitForTimeout(300);
			await page.getByRole("button", { name: "Run Extension" }).click();
			await page.waitForTimeout(2000);
			await page.screenshot({ path: path.join(outDir, "analysis_2d.png") });
		});

		test("molecule_builder", async ({ page }) => {
			const room = `${prefix}-molecule-builder`;
			PY(`
import ase
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms())
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Modifier tools" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("modifiers Method").click();
			await page.getByRole("option", { name: "AddFromSMILES" }).click();
			await page.waitForTimeout(500);
			await page
				.getByPlaceholder("Enter SMILES notation")
				.fill("NC(Cc1ccccc1)C(=O)O");
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Run Extension" }).click();
			await page.waitForTimeout(2000);
			await page.screenshot({
				path: path.join(outDir, "molecule_builder.png"),
			});
		});

		test("molecule_builder_editor", async ({ page }) => {
			const room = `${prefix}-molecule-builder-editor`;
			PY(`
import ase
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms())
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Modifier tools" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("modifiers Method").click();
			await page.getByRole("option", { name: "AddFromSMILES" }).click();
			await page.waitForTimeout(500);
			await page
				.getByPlaceholder("Enter SMILES notation")
				.fill("NC(Cc1ccccc1)C(=O)O");
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Run Extension" }).click();
			await page.waitForTimeout(2000);
			await page.getByPlaceholder("Enter SMILES notation").fill("CC(N)C(=O)O");
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Draw", exact: true }).click();
			await page
				.getByRole("dialog", { name: "Molecular Structure Editor" })
				.waitFor({ state: "visible" });
			await page.waitForTimeout(1000);
			await page.screenshot({
				path: path.join(outDir, "molecule_builder_editor.png"),
			});
		});

		test("rooms", async ({ page }) => {
			const room1 = `${prefix}-demo-simulation`;
			const room2 = `${prefix}-analysis-results`;
			const room3 = `${prefix}-molecule-library`;
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
v1 = ZnDraw(url='${BASE_URL}', room='${room1}')
del v1[:]
v1.extend([src[0] for _ in range(50)])
v2 = ZnDraw(url='${BASE_URL}', room='${room2}')
del v2[:]
v2.extend([src[0] for _ in range(25)])
v3 = ZnDraw(url='${BASE_URL}', room='${room3}')
del v3[:]
v3.extend([src[0] for _ in range(10)])
`);
			await page.goto(`${BASE_URL}/rooms`);
			await page.waitForTimeout(1500);
			await page.screenshot({ path: path.join(outDir, "rooms.png") });
		});

		test("file_browser", async ({ page }) => {
			const room = `${prefix}-file-browser`;
			cloneBmim(room);
			const bg = spawnPY(`
import time
import fsspec
from zndraw import ZnDraw

vis = ZnDraw(url='${BASE_URL}', room='${room}')
fs = fsspec.filesystem('memory')
fs.mkdir('/simulations')
fs.mkdir('/simulations/bmim_bf4')
with fs.open('/simulations/bmim_bf4/trajectory.xyz', 'w') as f:
    f.write('5\\nframe 0\\nH 0 0 0\\nH 1 0 0\\nH 0 1 0\\nH 1 1 0\\nH 0 0 1\\n')
with fs.open('/simulations/bmim_bf4/energies.csv', 'w') as f:
    f.write('step,energy\\n1,-100.5\\n')
fs.mkdir('/structures')
with fs.open('/structures/water.xyz', 'w') as f:
    f.write('3\\nwater\\nO 0 0 0\\nH 0.96 0 0\\nH -0.24 0.93 0\\n')
with fs.open('/structures/ethanol.xyz', 'w') as f:
    f.write('9\\nethanol\\nC 0 0 0\\nC 1.5 0 0\\nO 2.3 0 0\\nH 0 1 0\\nH 0 -1 0\\nH 0 0 1\\nH 1.5 1 0\\nH 1.5 -1 0\\nH 3.2 0 0\\n')
with fs.open('/molecule.xyz', 'w') as f:
    f.write('3\\nbenzene subset\\nC 0 0 0\\nC 1.4 0 0\\nC 0.7 1.2 0\\n')
vis.register_fs(fs, name='project-data')
time.sleep(60)
`);
			try {
				await waitForBgReady(10000);
				for (let attempt = 0; attempt < 3; attempt++) {
					await page.goto(`${BASE_URL}/rooms/${room}/files`);
					await page.waitForTimeout(2000);
					const noFs = page.getByText(
						"No filesystems are currently registered",
					);
					if (!(await noFs.isVisible().catch(() => false))) break;
					if (attempt < 2) await waitForBgReady(3000);
				}
				// Navigate into structures/ to show files alongside dirs
				await page.getByRole("button", { name: /structures/i }).click();
				await page.waitForTimeout(1500);
				await page.screenshot({ path: path.join(outDir, "file_browser.png") });
			} finally {
				bg.kill();
			}
		});

		test("camera", async ({ page }) => {
			const room = `${prefix}-camera`;
			PY(`
from zndraw import ZnDraw
from zndraw.geometries import Curve, Camera
from zndraw.transformations import CurveAttachment
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(src[0])
vis.geometries["cam_pos"] = Curve(
    position=[(25, 15, 25), (30, 20, 15), (25, 15, 5)],
    color="#3498db",
)
vis.geometries["cam_target"] = Curve(
    position=[(10, 10, 10), (12, 10, 10), (10, 10, 10)],
    color="#e74c3c",
)
vis.geometries["camera"] = Camera(
    position=CurveAttachment(geometry_key="cam_pos", progress=0.5),
    target=CurveAttachment(geometry_key="cam_target", progress=0.5),
    fov=60,
    helper_visible=True,
    helper_color="#00ff00",
)
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1500);
			await page.getByRole("button", { name: "Manage geometries" }).click();
			await page.waitForTimeout(1000);
			await page.screenshot({ path: path.join(outDir, "camera.png") });
		});

		test("python_connection", async ({ page }) => {
			const room = `${prefix}-python-connection`;
			cloneBmim(room);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await waitForScene(page);
			await page.getByRole("button", { name: "show connection info" }).click();
			await page.waitForTimeout(500);
			await page.screenshot({
				path: path.join(outDir, "python_connection.png"),
			});
		});

		test("chat", async ({ page }) => {
			const room = `${prefix}-chat`;
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(src[0])
vis.log("Welcome to ZnDraw! This chat supports **markdown** formatting.")
vis.log('Here\\'s a code example:\\n\\n\`\`\`py\\nprint("Hello ZnDraw!")\\n\`\`\`')
vis.log(
    "The potential energy is given by:\\n\\n"
    "$$E = \\\\sum_{i<j} 4\\\\varepsilon \\\\left[ "
    "\\\\left(\\\\frac{\\\\sigma}{r_{ij}}\\\\right)^{12} - "
    "\\\\left(\\\\frac{\\\\sigma}{r_{ij}}\\\\right)^{6} \\\\right]$$"
)
vis.log("You can also use inline math like $E = mc^2$ in your messages.")
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1000);
			await page.getByRole("button", { name: "toggle chat" }).click();
			await page.waitForTimeout(500);
			await page.screenshot({ path: path.join(outDir, "chat.png") });
		});

		test("custom_modifier", async ({ page }) => {
			const room = `${prefix}-custom-modifier`;
			const bg = spawnPY(`
import time
from pydantic import Field
from zndraw import ZnDraw
from zndraw.extensions import Extension, Category
import typing as t

class ScaleAtoms(Extension):
    """Scale atom positions by a factor."""
    category: t.ClassVar[Category] = Category.MODIFIER
    factor: float = Field(
        1.5, ge=0.1, le=5.0, description="Scale factor",
        json_schema_extra={"format": "range"},
    )
    center_first: bool = Field(
        True, description="Center atoms before scaling",
        json_schema_extra={"format": "checkbox"},
    )
    def run(self, vis, **kwargs):
        pass

src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(src[0])
vis.register_job(ScaleAtoms)
time.sleep(60)
`);
			try {
				await waitForBgReady();
				await page.goto(`${BASE_URL}/rooms/${room}`);
				await page.waitForTimeout(1000);
				await page.getByRole("button", { name: "Modifier tools" }).click();
				await page.waitForTimeout(500);
				await page.getByLabel("modifiers Method").click();
				await page.getByRole("option", { name: "ScaleAtoms" }).click();
				await page.waitForTimeout(500);
				await page.screenshot({
					path: path.join(outDir, "custom_modifier.png"),
				});
			} finally {
				bg.kill();
			}
		});

		test("progress_tracker", async ({ page }) => {
			const room = `${prefix}-progress-tracker`;
			cloneBmim(room);
			const bg = spawnPY(`
import time
from zndraw import ZnDraw
from zndraw.tqdm import ZnDrawTqdm
vis = ZnDraw(url='${BASE_URL}', room='${room}')
pbar = ZnDrawTqdm(total=100, vis=vis, description="Geometry Optimization", unit="steps")
pbar.update(42)
time.sleep(60)
pbar.close()
`);
			try {
				await waitForBgReady(5000);
				await page.goto(`${BASE_URL}/rooms/${room}`);
				await page.waitForTimeout(2000);
				await page.screenshot({
					path: path.join(outDir, "progress_tracker.png"),
				});
			} finally {
				bg.kill();
			}
		});

		test("locked_room", async ({ page }) => {
			const room = `${prefix}-locked-room`;
			cloneBmim(room);
			const bg = spawnPY(`
import time
from zndraw import ZnDraw
vis = ZnDraw(url='${BASE_URL}', room='${room}')
with vis.get_lock(msg="Uploading trajectory data..."):
    time.sleep(60)
`);
			try {
				await waitForBgReady(5000);
				await page.goto(`${BASE_URL}/rooms/${room}`);
				await page.waitForTimeout(2000);
				await page.screenshot({ path: path.join(outDir, "locked_room.png") });
			} finally {
				bg.kill();
			}
		});

		test("dynamic_properties_dropdown", async ({ page }) => {
			const room = `${prefix}-dynamic-props`;
			PY(`
import numpy as np
import ase
from zndraw import ZnDraw
from zndraw.geometries import Arrow
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
atoms = ase.Atoms("H4", positions=[(0, 0, 0), (2, 0, 0), (0, 2, 0), (2, 2, 0)])
atoms.arrays["forces"] = np.array(
    [[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0]], dtype=float,
)
vis.append(atoms)
vis.geometries["force_arrows"] = Arrow(
    position="arrays.positions", direction="arrays.forces",
    color=["#ff6600"], radius=0.1,
)
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(2000);
			await page.getByRole("button", { name: "Manage geometries" }).click();
			await page.waitForTimeout(500);
			await page.getByRole("row", { name: "force_arrows" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("Position").click();
			await page.waitForTimeout(500);
			await page.screenshot({
				path: path.join(outDir, "dynamic_properties_dropdown.png"),
			});
		});

		test("property_inspector", async ({ page }) => {
			const room = `${prefix}-property-inspector`;
			cloneEF(room, 10);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(2000);
			await page.keyboard.press("i");
			await page.waitForTimeout(500);
			await page.getByRole("button", { name: "Application settings" }).click();
			await page.waitForTimeout(500);
			await page.getByLabel("Settings Category").click();
			await page.waitForTimeout(300);
			await page.getByRole("option", { name: "Property Inspector" }).click();
			await page.waitForTimeout(500);
			await page.getByText("arrays.positions").click();
			await page.waitForTimeout(300);
			await page.getByText("calc.energy").click();
			await page.waitForTimeout(500);
			await page.mouse.move(950, 280);
			await page.waitForTimeout(500);
			await page.screenshot({
				path: path.join(outDir, "property_inspector.png"),
			});
		});

		test("chat_frame_reference", async ({ page }) => {
			const room = `${prefix}-chat-frame-ref`;
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.extend([src[0] for _ in range(20)])
vis.log("Initial structure loaded at @0")
vis.log("Check the transition at @10 and compare with @15")
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1000);
			await page.getByRole("button", { name: "toggle chat" }).click();
			await page.waitForTimeout(500);
			await page.screenshot({
				path: path.join(outDir, "chat_frame_reference.png"),
			});
		});

		test("chat_progress", async ({ page }) => {
			const room = `${prefix}-chat-progress`;
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(src[0])
vis.log(
    "Processing simulation data:\\n\\n"
    "\`\`\`progress\\n"
    "description: Analyzing frames\\n"
    "value: 75\\n"
    "max: 100\\n"
    "color: success\\n"
    "\`\`\`"
)
vis.log(
    "Waiting for calculation:\\n\\n"
    "\`\`\`progress\\n"
    "description: Optimizing geometry...\\n"
    "color: primary\\n"
    "\`\`\`"
)
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1000);
			await page.getByRole("button", { name: "toggle chat" }).click();
			await page.waitForTimeout(500);
			await page.screenshot({ path: path.join(outDir, "chat_progress.png") });
		});

		test("editing_mode", async ({ page }) => {
			const room = `${prefix}-editing-mode`;
			PY(`
import ase
from zndraw import ZnDraw
from zndraw.geometries import Box, Sphere
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms())
vis.geometries["box"] = Box(
    position=[(0, 2, 0)], size=[(4, 4, 4)], color=["#3498db"],
)
vis.geometries["sphere"] = Sphere(
    position=[(8, 2, 0)], radius=[2.0], color=["#e74c3c"],
)
vis.selections["box"] = [0]
vis.selections["sphere"] = [0]
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1500);
			await page.keyboard.press("e");
			await page.waitForTimeout(500);
			await page.screenshot({ path: path.join(outDir, "editing_mode.png") });
		});

		test("editing_axis_constraint", async ({ page }) => {
			const room = `${prefix}-editing-axis`;
			PY(`
import ase
from zndraw import ZnDraw
from zndraw.geometries import Box
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms())
vis.geometries["box"] = Box(
    position=[(0, 2, 0)], size=[(4, 4, 4)], color=["#3498db"],
)
vis.selections["box"] = [0]
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1500);
			await page.keyboard.press("e");
			await page.waitForTimeout(500);
			await page.keyboard.down("x");
			await page.waitForTimeout(300);
			await page.screenshot({
				path: path.join(outDir, "editing_axis_constraint.png"),
			});
			await page.keyboard.up("x");
		});

		test("curve_editing", async ({ page }) => {
			const room = `${prefix}-curve-editing`;
			PY(`
import ase
from zndraw import ZnDraw
from zndraw.geometries import Curve
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(ase.Atoms())
vis.geometries["curve"] = Curve(
    position=[(-6, 0, -6), (-3, 4, -3), (0, 0, 0), (3, 4, 3), (6, 0, 6)],
    color="#2ecc71",
    marker={"enabled": True, "size": 0.15, "opacity": 1.0},
    virtual_marker={"enabled": True, "size": 0.1, "opacity": 0.6},
)
vis.selections["curve"] = [1, 3]
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1500);
			await page.keyboard.press("e");
			await page.waitForTimeout(500);
			await page.screenshot({ path: path.join(outDir, "curve_editing.png") });
		});

		test("chat_smiles", async ({ page }) => {
			const room = `${prefix}-chat-smiles`;
			PY(`
from zndraw import ZnDraw
src = ZnDraw(url='${BASE_URL}', room='${prefix}-bmim-source')
vis = ZnDraw(url='${BASE_URL}', room='${room}')
del vis[:]
vis.append(src[0])
vis.log("Here's ethanol:\\n\\n\`\`\`smiles\\nCCO\\n\`\`\`")
vis.log("And here's caffeine:\\n\\n\`\`\`smiles\\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C\\n\`\`\`")
`);
			await page.goto(`${BASE_URL}/rooms/${room}`);
			await page.waitForTimeout(1500);
			await page.getByRole("button", { name: "toggle chat" }).click();
			await page.waitForTimeout(1000);
			await page.screenshot({ path: path.join(outDir, "chat_smiles.png") });
		});
	});
}
