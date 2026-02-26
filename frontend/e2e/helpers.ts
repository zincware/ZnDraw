import { execSync, spawn, type ChildProcess } from "child_process";
import { writeFileSync, unlinkSync } from "fs";
import { tmpdir } from "os";
import { join } from "path";
import type { Page, Browser, BrowserContext } from "@playwright/test";

export const BASE_URL = process.env.ZNDRAW_URL || "http://localhost:8000";

/** Run a zndraw-cli command and return stdout. */
export function CLI(cmd: string): string {
	return execSync(`uv run zndraw-cli --url ${BASE_URL} ${cmd}`, {
		encoding: "utf-8",
	});
}

/** Write a Python script to a temp file, run it, return stdout. */
export function PY(code: string): string {
	const tmp = join(
		tmpdir(),
		`zndraw-e2e-${Date.now()}-${Math.random().toString(36).slice(2, 8)}.py`,
	);
	writeFileSync(tmp, code);
	try {
		return execSync(`uv run python ${tmp}`, { encoding: "utf-8" });
	} finally {
		unlinkSync(tmp);
	}
}

/** Wait for the 3D scene canvas to be visible. */
export async function waitForScene(page: Page): Promise<void> {
	await page.waitForTimeout(1000);
	await page.waitForSelector("canvas", { state: "visible", timeout: 15000 });
}

/** Spawn a Python script as a background process. Returns the child process. */
export function spawnPY(code: string): ChildProcess {
	const tmp = join(
		tmpdir(),
		`zndraw-e2e-bg-${Date.now()}-${Math.random().toString(36).slice(2, 8)}.py`,
	);
	writeFileSync(tmp, code);
	const child = spawn("uv", ["run", "python", tmp], {
		stdio: ["pipe", "pipe", "pipe"],
		detached: false,
	});
	child.stderr?.on("data", (data: Buffer) => {
		console.error(`[bg-py stderr] ${data.toString().trim()}`);
	});
	child.on("exit", () => {
		try {
			unlinkSync(tmp);
		} catch {}
	});
	return child;
}

/** Wait for a background process to be ready (registration propagation). */
export async function waitForBgReady(ms: number = 8000): Promise<void> {
	await new Promise((r) => setTimeout(r, ms));
}

/**
 * Create a fresh browser context with isolated auth.
 * Each context auto-creates a guest user on first page load
 * (the frontend calls POST /v1/auth/guest when no token is in localStorage).
 */
export async function createAuthContext(
	browser: Browser,
): Promise<BrowserContext> {
	return browser.newContext({
		viewport: { width: 1280, height: 720 },
	});
}
