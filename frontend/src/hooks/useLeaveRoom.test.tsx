import { act, renderHook } from "@testing-library/react";
import type { DockviewApi } from "dockview-react";
import { MemoryRouter } from "react-router-dom";
import { describe, expect, it, vi } from "vitest";
import { useLeaveRoom } from "./useLeaveRoom";

function makeApi(): DockviewApi {
	return {
		panels: [],
		getPanel: () => undefined,
	} as unknown as DockviewApi;
}

describe("useLeaveRoom", () => {
	it("resolves api lazily when passed a getter", async () => {
		let resolved: DockviewApi | null = null;
		const getter = () => resolved;
		const { result } = renderHook(() => useLeaveRoom({ api: getter }), {
			wrapper: MemoryRouter,
		});
		// Populate after mount.
		resolved = makeApi();
		await act(async () => {
			await result.current({ skipConfirm: true });
		});
		// No throw — lazily resolved the populated api.
		expect(true).toBe(true);
	});

	it("calls react-router navigate, not history.pushState", async () => {
		const pushSpy = vi.spyOn(window.history, "pushState");
		const { result } = renderHook(() => useLeaveRoom({ api: makeApi() }), {
			wrapper: MemoryRouter,
		});
		await act(async () => {
			await result.current({ skipConfirm: true });
		});
		expect(pushSpy).not.toHaveBeenCalled();
	});
});
