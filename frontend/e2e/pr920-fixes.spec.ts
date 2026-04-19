import { test, expect } from "@playwright/test";

const ROOM = "playwright-room";

test.describe("PR #920 fixes", () => {
	test("chat unread badge accumulates while hidden and clears on switch-back", async ({
		page,
		browser,
	}) => {
		await page.goto(`/rooms/${ROOM}`);
		await page.getByTestId("activity-icon-chat").click();
		await expect(
			page
				.getByTestId("sidebar-zone-right")
				.or(page.getByTestId("sidebar-zone-left")),
		).toBeVisible();
		// Switch away (hide chat).
		await page.getByTestId("activity-icon-rooms").click();
		// Send a message from a peer session.
		const ctx2 = await browser.newContext();
		const peer = await ctx2.newPage();
		await peer.goto(`/rooms/${ROOM}`);
		await peer.getByTestId("activity-icon-chat").click();
		const peerInput = peer
			.getByTestId("chat-input")
			.or(peer.locator('[data-testid="chat-input"] input'));
		await peerInput.fill("hello from peer");
		await peer.getByTestId("chat-send").click();
		await ctx2.close();
		// Unread badge becomes visible on the chat icon while we're elsewhere.
		const badge = page
			.getByTestId("activity-icon-chat")
			.locator(".MuiBadge-badge");
		await expect(badge).toBeVisible({ timeout: 10000 });
		// Switch back to chat — the unread reset effect runs and clears the badge.
		await page.getByTestId("activity-icon-chat").click();
		await expect(badge).toBeHidden();
	});
});
