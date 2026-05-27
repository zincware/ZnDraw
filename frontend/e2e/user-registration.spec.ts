import { test, expect } from "@playwright/test";
import { BASE_URL } from "./helpers";

const DISPLAY_NAME_RE = /^[a-z][a-z0-9-]{2,63}$/;

test.describe("User registration with display_name suggestion", () => {
	test("dialog pre-fills a server-suggested display name", async ({ page }) => {
		const suggestionRequest = page.waitForRequest((req) =>
			req.url().includes("/v1/users/available-display-name"),
		);

		await page.goto(BASE_URL);
		await page.getByRole("button", { name: /register/i }).click();

		await suggestionRequest;

		const displayNameField = page.getByLabel("Display name");
		await expect(displayNameField).toBeVisible();

		const value = await displayNameField.inputValue();
		expect(value).toMatch(DISPLAY_NAME_RE);
	});

	test("regenerate button fetches a new suggestion", async ({ page }) => {
		await page.goto(BASE_URL);
		await page.getByRole("button", { name: /register/i }).click();

		const displayNameField = page.getByLabel("Display name");
		await expect(displayNameField).toBeVisible();
		const firstValue = await displayNameField.inputValue();

		const regenerateRequest = page.waitForRequest((req) =>
			req.url().includes("/v1/users/available-display-name"),
		);
		await page
			.getByRole("button", { name: /regenerate suggestion/i })
			.click();
		await regenerateRequest;

		await expect
			.poll(async () => displayNameField.inputValue(), { timeout: 5_000 })
			.not.toEqual(firstValue);
		expect(await displayNameField.inputValue()).toMatch(DISPLAY_NAME_RE);
	});

	test("successful registration lands on a display-name room URL", async ({
		page,
	}) => {
		await page.goto(BASE_URL);
		await page.getByRole("button", { name: /register/i }).click();

		const email = `e2e-${Date.now()}-${Math.random()
			.toString(36)
			.slice(2, 8)}@example.com`;
		const password = "very-strong-passw0rd";

		const displayNameField = page.getByLabel("Display name");
		await expect(displayNameField).toBeVisible();
		const suggestedDisplayName = await displayNameField.inputValue();
		expect(suggestedDisplayName).toMatch(DISPLAY_NAME_RE);

		await page.getByLabel("Email").fill(email);
		await page.getByLabel("Password", { exact: true }).fill(password);
		await page.getByLabel("Confirm Password").fill(password);

		await page.getByRole("button", { name: /^register$/i }).click();

		await page.waitForURL(
			new RegExp(`/rooms/${suggestedDisplayName}/[a-zA-Z0-9_-]+`),
			{ timeout: 10_000 },
		);
	});
});
