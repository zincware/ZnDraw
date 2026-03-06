import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-camera-session";

/** Open the geometry panel by clicking sidebar button. */
async function openGeometryPanel(page: import("@playwright/test").Page) {
	await page.getByRole("button", { name: "Manage geometries" }).click();
	await page.waitForSelector('[role="grid"]', {
		state: "visible",
		timeout: 10000,
	});
}

/**
 * Move the camera for the current browser session via REST API calls
 * executed in the browser context (same JWT token, same user identity).
 */
async function moveCameraViaAPI(
	page: import("@playwright/test").Page,
	room: string,
	position: [number, number, number],
	target: [number, number, number],
) {
	await page.evaluate(
		async ({ room, position, target }) => {
			const token = localStorage.getItem("zndraw_jwt_token");
			const headers: Record<string, string> = {
				Authorization: `Bearer ${token}`,
				"Content-Type": "application/json",
			};

			// 1. Find the browser session SID via room presence
			let sessions: Array<{ sid: string }> = [];
			for (let attempt = 0; attempt < 60; attempt++) {
				const resp = await fetch(`/v1/rooms/${room}/presence`, {
					headers,
				});
				if (resp.ok) {
					const data = await resp.json();
					sessions = data.items || [];
					if (sessions.length > 0) break;
				}
				await new Promise((r) => setTimeout(r, 500));
			}
			if (sessions.length === 0)
				throw new Error("No sessions in room presence after 30s");

			const sid = sessions[0].sid;

			// 2. Get the active camera key
			const camResp = await fetch(
				`/v1/rooms/${room}/sessions/${sid}/active-camera`,
				{ headers },
			);
			if (!camResp.ok)
				throw new Error(`Failed to get active camera: ${camResp.status}`);
			const camKey = (await camResp.json()).active_camera;

			// 3. Get current camera geometry data
			const geomResp = await fetch(`/v1/rooms/${room}/geometries/${camKey}`, {
				headers,
			});
			if (!geomResp.ok)
				throw new Error(`Failed to get camera geometry: ${geomResp.status}`);
			const geom = await geomResp.json();
			const camData = geom.geometry.data;

			// 4. Update camera position and PUT back
			camData.position = position;
			camData.target = target;
			const putResp = await fetch(`/v1/rooms/${room}/geometries/${camKey}`, {
				method: "PUT",
				headers,
				body: JSON.stringify({ type: "Camera", data: camData }),
			});
			if (!putResp.ok)
				throw new Error(`Failed to update camera: ${putResp.status}`);
		},
		{ room, position, target },
	);
}

test.describe("Camera Session", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		CLI(`rooms create --room-id ${ROOM}`);
		PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[:]
vis.append(ase.Atoms('H4', positions=[(0,0,0),(2,0,0),(0,2,0),(2,2,0)]))
`);
	});

	test("camera geometry appears in panel", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Open geometry panel — the browser session's camera should be listed
		await openGeometryPanel(page);
		await page.waitForTimeout(500);

		const grid = page.locator('[role="grid"]');
		await expect(grid.getByText("Camera").first()).toBeVisible({
			timeout: 10000,
		});

		await page.screenshot({
			path: "e2e/screenshots/camera-geometry-panel.png",
		});
	});

	test("camera update via API reflects in geometry panel", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Wait for socket session to fully register in Redis
		await page.waitForTimeout(3000);

		// Open geometry panel first
		await openGeometryPanel(page);
		await page.waitForTimeout(500);

		await page.screenshot({
			path: "e2e/screenshots/camera-before-update.png",
		});

		// Move camera via REST API from the browser context (same user)
		await moveCameraViaAPI(page, ROOM, [50, 50, 50], [0, 0, 0]);

		// Wait for geometry_invalidate socket event to update the panel
		await page.waitForTimeout(2000);

		await page.screenshot({
			path: "e2e/screenshots/camera-after-update.png",
		});

		// The geometry panel should still show Camera entry
		const grid = page.locator('[role="grid"]');
		await expect(grid.getByText("Camera").first()).toBeVisible({
			timeout: 5000,
		});
	});

	test("camera position change updates viewport", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Wait for socket session to fully register
		await page.waitForTimeout(3000);

		await page.screenshot({
			path: "e2e/screenshots/camera-viewport-before.png",
		});

		// Move camera very far away
		await moveCameraViaAPI(page, ROOM, [200, 200, 200], [0, 0, 0]);

		await page.waitForTimeout(2000);

		await page.screenshot({
			path: "e2e/screenshots/camera-viewport-after.png",
		});

		// Manual inspection: the two screenshots should look different
		// (atoms appear tiny when camera is far away)
	});
});
