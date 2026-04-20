import { beforeEach, describe, expect, it } from "vitest";
import { useAppStore } from "../../../store";

describe("activityBarSlice", () => {
	beforeEach(() => {
		useAppStore.getState().resetLayout();
	});

	it("initial state has all tool panels on left bar", () => {
		const s = useAppStore.getState();
		expect(s.leftBarIcons).toEqual([
			"selections",
			"modifiers",
			"analysis",
			"geometries",
			"plots-browser",
			"rooms",
			"filesystem",
			"chat",
		]);
		expect(s.rightBarIcons).toEqual([]);
		expect(s.bottomBarIcons).toEqual([]);
	});

	it("initial state has all active panels as null", () => {
		const s = useAppStore.getState();
		expect(s.activeLeft).toBeNull();
		expect(s.activeRight).toBeNull();
		expect(s.activeBottom).toBeNull();
	});

	it("moveIconToBar moves an icon from one bar to another", () => {
		const s = useAppStore.getState();
		s.moveIconToBar("filesystem", "right");
		const next = useAppStore.getState();
		expect(next.leftBarIcons).not.toContain("filesystem");
		expect(next.rightBarIcons).toContain("filesystem");
	});

	it("moveIconToBar with index inserts at specified position", () => {
		const s = useAppStore.getState();
		s.moveIconToBar("filesystem", "right", 0);
		const next = useAppStore.getState();
		expect(next.rightBarIcons[0]).toBe("filesystem");
	});

	it("moveIconToBar clears active state in source bar if icon was active", () => {
		const s = useAppStore.getState();
		s.toggleActive("left", "filesystem");
		expect(useAppStore.getState().activeLeft).toBe("filesystem");
		s.moveIconToBar("filesystem", "right");
		const next = useAppStore.getState();
		expect(next.activeLeft).toBeNull();
		expect(next.rightBarIcons).toContain("filesystem");
	});

	it("dropIconOnPanel removes icon from bar and sets target bar active", () => {
		const s = useAppStore.getState();
		const initialLeft = [...s.leftBarIcons];
		const victim = initialLeft[0]; // "selections"
		s.dropIconOnPanel(victim, "right");
		const next = useAppStore.getState();
		expect(next.leftBarIcons).not.toContain(victim);
		expect(next.rightBarIcons).toContain(victim);
		expect(next.activeRight).toBe(victim);
	});

	it("dropIconOnPanel clears active in source bar if icon was active", () => {
		const s = useAppStore.getState();
		s.toggleActive("left", "selections");
		expect(useAppStore.getState().activeLeft).toBe("selections");
		s.dropIconOnPanel("selections", "right");
		const next = useAppStore.getState();
		expect(next.activeLeft).toBeNull();
		expect(next.activeRight).toBe("selections");
	});

	it("resetLayout restores initial state", () => {
		const s = useAppStore.getState();
		s.moveIconToBar("filesystem", "right");
		s.toggleActive("left", "selections");
		s.resetLayout();
		const next = useAppStore.getState();
		expect(next.leftBarIcons).toEqual([
			"selections",
			"modifiers",
			"analysis",
			"geometries",
			"plots-browser",
			"rooms",
			"filesystem",
			"chat",
		]);
		expect(next.rightBarIcons).toEqual([]);
		expect(next.bottomBarIcons).toEqual([]);
		expect(next.activeLeft).toBeNull();
		expect(next.activeRight).toBeNull();
		expect(next.activeBottom).toBeNull();
	});

	it("toggleActive opens a panel if closed", () => {
		const s = useAppStore.getState();
		s.toggleActive("left", "selections");
		const next = useAppStore.getState();
		expect(next.activeLeft).toBe("selections");
	});

	it("toggleActive closes a panel if already open", () => {
		const s = useAppStore.getState();
		s.toggleActive("left", "selections");
		expect(useAppStore.getState().activeLeft).toBe("selections");
		s.toggleActive("left", "selections");
		const next = useAppStore.getState();
		expect(next.activeLeft).toBeNull();
	});

	it("toggleActive switches active panel on same bar", () => {
		const s = useAppStore.getState();
		s.toggleActive("left", "selections");
		expect(useAppStore.getState().activeLeft).toBe("selections");
		s.toggleActive("left", "modifiers");
		const next = useAppStore.getState();
		expect(next.activeLeft).toBe("modifiers");
	});

	it("toggleActive maintains independent active state across bars", () => {
		const s = useAppStore.getState();
		s.moveIconToBar("filesystem", "right");
		s.toggleActive("left", "selections");
		s.toggleActive("right", "filesystem");
		const next = useAppStore.getState();
		expect(next.activeLeft).toBe("selections");
		expect(next.activeRight).toBe("filesystem");
	});
});
