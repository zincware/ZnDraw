import type { DockviewApi } from "dockview-react";
import { beforeEach, describe, expect, it } from "vitest";
import { useDockviewApi } from "../dockviewApiStore";

describe("dockviewApiStore", () => {
	beforeEach(() => {
		useDockviewApi.getState().setApi(null);
	});

	it("defaults to null", () => {
		expect(useDockviewApi.getState().api).toBeNull();
	});

	it("setApi updates state", () => {
		const fakeApi = { panels: [] } as unknown as DockviewApi;
		useDockviewApi.getState().setApi(fakeApi);
		expect(useDockviewApi.getState().api).toBe(fakeApi);
	});

	it("setApi(null) clears", () => {
		const fakeApi = {} as unknown as DockviewApi;
		useDockviewApi.getState().setApi(fakeApi);
		useDockviewApi.getState().setApi(null);
		expect(useDockviewApi.getState().api).toBeNull();
	});
});
