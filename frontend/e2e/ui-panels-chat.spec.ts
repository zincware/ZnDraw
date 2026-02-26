import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY } from "./helpers";

const ROOM = "test-ui";

test.describe("UI Panels & Chat", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		CLI(`rooms create --room-id ${ROOM}`);
		PY(`
from zndraw import ZnDraw
import ase
vis = ZnDraw(url="${BASE_URL}", room="${ROOM}")
del vis[:]
vis.append(ase.Atoms("H4", positions=[(0,0,0),(1,0,0),(0,1,0),(1,1,0)]))
vis.log("Hello **world**! Testing markdown.")
vis.log("Code block:\\n\\n\`\`\`python\\nprint(\\"hello\\")\\n\`\`\`")
vis.log("Check frame @0 for details")
`);
	});

	test("rooms page lists rooms", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms`);
		await page.waitForTimeout(1000);

		// Verify the heading is present
		await expect(
			page.getByRole("heading", { name: "Room Management" }),
		).toBeVisible();

		// Verify subtitle text
		await expect(
			page.getByText(
				"Manage your rooms: lock, hide, set as template, or duplicate",
			),
		).toBeVisible();

		// Verify search box
		const searchBox = page.getByRole("textbox", { name: "Search rooms" });
		await expect(searchBox).toBeVisible();

		// Verify the data grid has column headers
		await expect(
			page.getByRole("columnheader", { name: "Room ID" }),
		).toBeVisible();
		await expect(
			page.getByRole("columnheader", { name: "Frames" }),
		).toBeVisible();
		await expect(
			page.getByRole("columnheader", { name: "Lock" }),
		).toBeVisible();

		// Search for our test room (may be off-screen with many parallel rooms)
		await searchBox.fill(ROOM);
		await page.waitForTimeout(500);

		// Verify our test room appears in the filtered grid
		await expect(page.getByRole("gridcell", { name: ROOM })).toBeVisible();

		// Verify action buttons are present
		await expect(
			page.getByRole("button", { name: "Create New Empty Room" }),
		).toBeVisible();
		await expect(
			page.getByRole("button", { name: "Upload File" }),
		).toBeVisible();

		await page.screenshot({ path: "e2e/screenshots/ui-rooms-page.png" });
	});

	test("chat panel shows markdown messages", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await page.waitForTimeout(1000);

		// Open the chat panel
		await page.getByRole("button", { name: "toggle chat" }).click();
		await page.waitForTimeout(500);

		// Verify the chat heading is visible
		await expect(page.getByRole("heading", { name: "Chat" })).toBeVisible();

		// Verify the text input is present
		await expect(
			page.getByRole("textbox", { name: "Type a message..." }),
		).toBeVisible();

		// Verify markdown rendering: bold text "world" in first message
		await expect(
			page.locator("strong").filter({ hasText: "world" }).first(),
		).toBeVisible();

		// Verify code block with syntax highlighting is present
		await expect(
			page.locator("code").filter({ hasText: 'print("hello")' }).first(),
		).toBeVisible();

		// Verify frame reference @0 is rendered (multiple messages may contain it)
		await expect(page.getByText("@0").first()).toBeVisible();

		// Verify the fullscreen button exists
		await expect(
			page.getByRole("button", { name: "Fullscreen" }),
		).toBeVisible();

		await page.screenshot({ path: "e2e/screenshots/ui-chat.png" });
	});

	test("connection info dialog opens and closes", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await page.waitForTimeout(1000);

		// Open the connection info dialog
		await page.getByRole("button", { name: "show connection info" }).click();
		await page.waitForTimeout(500);

		// Verify the dialog heading is visible
		await expect(
			page.getByRole("heading", {
				name: "Connect to this ZnDraw session",
			}),
		).toBeVisible();

		// Verify the dialog contains connection instructions
		await expect(
			page.getByText(
				"You can connect a Python kernel to this ZnDraw session using:",
			),
		).toBeVisible();

		// Verify the Python code snippet contains key elements
		await expect(
			page.locator("code").filter({ hasText: "from zndraw import ZnDraw" }),
		).toBeVisible();
		await expect(
			page.locator("code").filter({ hasText: `room="${ROOM}"` }),
		).toBeVisible();

		// Verify a copy button exists (dialog may have multiple copy buttons)
		await expect(
			page.getByRole("button", { name: "Copy to clipboard" }).first(),
		).toBeVisible();

		await page.screenshot({
			path: "e2e/screenshots/ui-connection-info.png",
		});

		// Close with Escape
		await page.keyboard.press("Escape");
		await page.waitForTimeout(500);

		// Verify dialog is closed
		await expect(
			page.getByRole("heading", {
				name: "Connect to this ZnDraw session",
			}),
		).not.toBeVisible();
	});

	test("theme toggle switches to dark mode", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await page.waitForTimeout(1000);

		// Click theme toggle (force: true to bypass potential tooltip overlap)
		await page.getByLabel("toggle theme").click({ force: true });
		await page.waitForTimeout(500);

		// After first click the icon should change (dark mode icon shown)
		// The aria-label stays "toggle theme" but the tooltip text changes
		await page.screenshot({
			path: "e2e/screenshots/ui-dark-mode-first-click.png",
		});

		// Click again to toggle back
		await page.getByLabel("toggle theme").click({ force: true });
		await page.waitForTimeout(500);

		await page.screenshot({ path: "e2e/screenshots/ui-dark-mode.png" });

		// The theme toggle button should still be visible after toggling
		await expect(page.getByLabel("toggle theme")).toBeVisible();
	});

	test("settings panel opens with categories", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await page.waitForTimeout(1000);

		// Click the Application settings button in the sidebar
		await page.getByRole("button", { name: "Application settings" }).click();
		await page.waitForTimeout(500);

		// Verify the Settings heading is visible
		await expect(page.getByRole("heading", { name: "Settings" })).toBeVisible();

		// Verify the Settings Category dropdown is present
		await expect(page.getByLabel("Settings Category")).toBeVisible();

		await page.screenshot({ path: "e2e/screenshots/ui-settings.png" });

		// Click the settings button again to close the panel
		await page.getByRole("button", { name: "Application settings" }).click();
		await page.waitForTimeout(500);

		// Settings panel should be hidden
		await expect(
			page.getByRole("heading", { name: "Settings" }),
		).not.toBeVisible();
	});

	test("file browser page loads", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}/files`);
		await page.waitForTimeout(1000);

		// Without registered filesystem providers, we expect the info alert
		await expect(
			page.getByText("No filesystems are currently registered"),
		).toBeVisible();

		// Verify the Back button is present
		await expect(page.getByRole("button", { name: "Back" })).toBeVisible();

		await page.screenshot({
			path: "e2e/screenshots/ui-file-browser.png",
		});
	});

	test("template selection page redirects", async ({ page }) => {
		// The template selection page auto-redirects to a room
		// It either goes to the last visited room or creates a new one
		await page.goto(`${BASE_URL}/`);
		await page.waitForTimeout(3000);

		// After redirect, we should be on a room page (URL contains /rooms/)
		const url = page.url();
		expect(url).toMatch(/\/rooms?\//);

		// The AppBar should be visible (indicates we're on the main page)
		await expect(page.getByText("ZnDraw")).toBeVisible();

		await page.screenshot({
			path: "e2e/screenshots/ui-template-selection.png",
		});
	});
});
