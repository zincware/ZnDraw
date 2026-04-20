import { act, renderHook } from "@testing-library/react";
import type { DragEvent } from "react";
import { afterEach, describe, expect, it } from "vitest";
import { useAppStore } from "../../store";
import { useDragHover } from "../useDragHover";

// Minimal fake drag event that carries the panel MIME type so the hook
// accepts it (onDragEnter/onDragOver filter on DRAG_MIME).
const DRAG_MIME = "application/x-zndraw-panel-id";
const fakeDragEvent = {
	preventDefault() {},
	dataTransfer: {
		types: [DRAG_MIME],
		dropEffect: "",
	},
} as unknown as DragEvent;

describe("useDragHover", () => {
	afterEach(() => {
		// Reset store state so tests are independent.
		useAppStore.setState({ dragHoverBar: null });
	});

	it("sets dragHoverBar on enter", () => {
		const { result } = renderHook(() => useDragHover("left"));
		expect(result.current.isHovered).toBe(false);
		act(() => {
			result.current.dragHandlers.onDragEnter(fakeDragEvent);
		});
		expect(useAppStore.getState().dragHoverBar).toBe("left");
	});

	it("clears dragHoverBar on leave after microtask", async () => {
		const { result } = renderHook(() => useDragHover("right"));
		act(() => {
			result.current.dragHandlers.onDragEnter(fakeDragEvent);
		});
		expect(useAppStore.getState().dragHoverBar).toBe("right");
		act(() => {
			result.current.dragHandlers.onDragLeave();
		});
		await new Promise((r) => setTimeout(r, 5));
		expect(useAppStore.getState().dragHoverBar).toBeNull();
	});

	it("depth counter handles nested enter/leave", async () => {
		const { result } = renderHook(() => useDragHover("bottom"));
		act(() => {
			result.current.dragHandlers.onDragEnter(fakeDragEvent);
			result.current.dragHandlers.onDragEnter(fakeDragEvent);
			result.current.dragHandlers.onDragLeave();
		});
		expect(useAppStore.getState().dragHoverBar).toBe("bottom");
		act(() => {
			result.current.dragHandlers.onDragLeave();
		});
		await new Promise((r) => setTimeout(r, 5));
		expect(useAppStore.getState().dragHoverBar).toBeNull();
	});
});
