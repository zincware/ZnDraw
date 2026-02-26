import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene, spawnPY, waitForBgReady } from "./helpers";

const ROOM_EXT = "test-registration-ext";
const ROOM_FS = "test-registration-fs";
const ROOM_MOUNT = "test-registration-mount";

test.describe("Registration", () => {
	test.describe.configure({ mode: "serial" });

	test("register_job makes extension appear in modifier panel", async ({
		page,
	}) => {
		CLI(`rooms create --room-id ${ROOM_EXT}`);
		PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${ROOM_EXT}')
del vis[:]
vis.append(ase.Atoms('H2', positions=[(0,0,0),(0.74,0,0)]))
`);

		// Spawn a background Python process that registers a custom extension
		// and stays alive (the extension is only available while the process runs)
		const bg = spawnPY(`
import time
from zndraw import ZnDraw
from zndraw.extensions import Extension, Category
from pydantic import Field
import typing as t

class E2ETestModifier(Extension):
    """A test modifier registered by E2E tests."""
    category: t.ClassVar[Category] = Category.MODIFIER
    greeting: str = Field(default="hello", title="Greeting")

    def run(self, vis):
        pass

vis = ZnDraw(url='${BASE_URL}', room='${ROOM_EXT}')
vis.register_job(E2ETestModifier)
# Keep alive for 60 seconds (test will kill us when done)
time.sleep(60)
`);

		try {
			// Wait for the registration to propagate
			await waitForBgReady();

			await page.goto(`${BASE_URL}/rooms/${ROOM_EXT}`);
			await waitForScene(page);

			// Open modifiers panel
			await page
				.getByRole("button", { name: "Modifier tools" })
				.click();
			await page
				.getByRole("heading", { name: "modifiers" })
				.waitFor({ state: "visible", timeout: 10000 });

			// The custom extension should appear in the Method dropdown
			await page
				.getByRole("combobox", { name: "modifiers Method" })
				.click();
			await page.waitForTimeout(500);

			// Look for the custom extension in the dropdown options
			// The name may be displayed as "E2ETestModifier" or "E2E Test Modifier"
			const option = page
				.getByRole("option")
				.filter({ hasText: /E2E/i });
			await expect(option.first()).toBeVisible({ timeout: 15000 });

			await page.screenshot({
				path: "e2e/screenshots/registration-extension.png",
			});
		} finally {
			bg.kill();
		}
	});

	test("register_fs enables file browser with browse and load", async ({
		page,
	}) => {
		CLI(`rooms create --room-id ${ROOM_FS}`);
		PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${ROOM_FS}')
del vis[:]
vis.append(ase.Atoms('H2', positions=[(0,0,0),(0.74,0,0)]))
`);

		// Spawn background process that registers a filesystem with a 2-frame trajectory
		const bg = spawnPY(`
import time
import fsspec
from zndraw import ZnDraw

vis = ZnDraw(url='${BASE_URL}', room='${ROOM_FS}')
fs = fsspec.filesystem('memory')
# Create a directory with a 2-frame water trajectory
fs.mkdir('/testdir')
with fs.open('/testdir/water.xyz', 'w') as f:
    f.write('3\\nframe 0\\nO 0 0 0\\nH 0.96 0 0\\nH -0.24 0.93 0\\n')
    f.write('3\\nframe 1\\nO 0.1 0 0\\nH 1.06 0 0\\nH -0.14 0.93 0\\n')
vis.register_fs(fs, name='e2e-memfs')
# Keep alive to serve provider requests and tasks
time.sleep(120)
`);

		try {
			// Wait for registration to propagate, then retry with reload
			await waitForBgReady(10000);

			// Retry up to 3 times — the file browser fetches fs list on page load
			for (let attempt = 0; attempt < 3; attempt++) {
				await page.goto(`${BASE_URL}/rooms/${ROOM_FS}/files`);
				await page.waitForTimeout(2000);

				const noFs = page.getByText(
					"No filesystems are currently registered",
				);
				const isNoFs = await noFs.isVisible().catch(() => false);
				if (!isNoFs) break;

				// Still showing "no filesystems" — wait and retry
				if (attempt < 2) await waitForBgReady(3000);
			}

			// The "No filesystems" message should NOT be visible
			await expect(
				page.getByText("No filesystems are currently registered"),
			).not.toBeVisible({ timeout: 5000 });

			// The filesystem name should appear (single provider = inline text)
			await expect(page.getByText("e2e-memfs")).toBeVisible({
				timeout: 10000,
			});

			await page.screenshot({
				path: "e2e/screenshots/registration-filesystem.png",
			});

			// --- Browse: navigate into /testdir ---
			// The root listing should show the "testdir" directory
			await expect(
				page.getByRole("button", { name: /testdir/i }),
			).toBeVisible({ timeout: 15000 });
			await page.getByRole("button", { name: /testdir/i }).click();

			// --- Browse: verify file listing ---
			// The file "water.xyz" should appear in the directory listing
			await expect(
				page.getByRole("button", { name: /water\.xyz/i }),
			).toBeVisible({ timeout: 15000 });

			await page.screenshot({
				path: "e2e/screenshots/registration-filesystem-browse.png",
			});

			// --- Load: click the file to open LoadFileDialog ---
			await page.getByRole("button", { name: /water\.xyz/i }).click();
			await expect(
				page.getByText("Load File from Remote Filesystem"),
			).toBeVisible({ timeout: 5000 });
			await expect(page.getByText("File: water.xyz")).toBeVisible();

			// --- Load: submit with default "Current room" target ---
			await page.getByRole("button", { name: "Load" }).click();

			// After submitting, the UI navigates back to the room page
			await page.waitForURL(`**/rooms/${ROOM_FS}`, { timeout: 15000 });
			await waitForScene(page);

			// The worker executes LoadFile in the background — wait for frames
			// Room had 1 frame (H2), loading adds 2 frames (water) → total 3
			await expect(page.getByText("/ 3")).toBeVisible({
				timeout: 30000,
			});

			await page.screenshot({
				path: "e2e/screenshots/registration-filesystem-loaded.png",
			});
		} finally {
			bg.kill();
		}
	});

	test("mount serves frames on demand", async ({ page }) => {
		CLI(`rooms create --room-id ${ROOM_MOUNT}`);

		// Spawn background process that mounts a frame source
		const bg = spawnPY(`
import time
import ase
from zndraw import ZnDraw

vis = ZnDraw(url='${BASE_URL}', room='${ROOM_MOUNT}')
# Room must be empty for mount
del vis[:]
# Create a list of 50 frames
frames = [ase.Atoms('H2', positions=[(0,0,0),(0.74+i*0.01,0,0)]) for i in range(50)]
vis.mount(frames)
# Keep alive to serve frame requests
time.sleep(60)
`);

		try {
			// Wait for mount
			await waitForBgReady();

			await page.goto(`${BASE_URL}/rooms/${ROOM_MOUNT}`);
			await waitForScene(page);

			// Frame count should show 50
			await expect(page.getByText("/ 50")).toBeVisible({
				timeout: 15000,
			});

			// Navigate to frame 25
			CLI(`step set ${ROOM_MOUNT} 25`);
			await expect(page.getByText("26 / 50")).toBeVisible({
				timeout: 10000,
			});

			await page.screenshot({
				path: "e2e/screenshots/registration-mount.png",
			});
		} finally {
			bg.kill();
		}
	});
});
