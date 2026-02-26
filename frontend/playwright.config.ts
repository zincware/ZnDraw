import { defineConfig } from "@playwright/test";

export default defineConfig({
	testDir: "./e2e",
	timeout: 90_000,
	expect: { timeout: 15_000 },
	fullyParallel: true,
	workers: 4,
	retries: 0,
	reporter: [["html", { open: "never" }]],
	use: {
		baseURL: process.env.ZNDRAW_URL || "http://localhost:8000",
		screenshot: "only-on-failure",
		trace: "retain-on-failure",
		viewport: { width: 1280, height: 720 },
	},
	projects: [
		{
			name: "chromium",
			use: { browserName: "chromium", headless: false },
		},
	],
});
