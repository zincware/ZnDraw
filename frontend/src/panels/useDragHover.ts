import type { DragEvent } from "react";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import type { BarPosition } from "./registry";

const DRAG_MIME = "application/x-zndraw-panel-id";

/**
 * Drag-hover state machine for an ActivityBar / SidebarZone / BottomZone.
 *
 * Tracks `dragDepth` to reliably detect leaving the subtree — nested
 * children generate sibling enter/leave pairs that would otherwise
 * confuse a naive implementation. Uses setTimeout(0) on the zero
 * crossing so an immediate re-enter on an adjacent child re-arms the
 * counter before we clear ``dragHoverBar``.
 *
 * Only responds to drags carrying the panel-id MIME type so that
 * file-system and other browser drags are ignored.
 */
export function useDragHover(position: BarPosition) {
	const dragDepth = useRef(0);
	const dragHoverBar = useAppStore((s) => s.dragHoverBar);
	const setDragHoverBar = useAppStore((s) => s.setDragHoverBar);
	const isHovered = dragHoverBar === position;

	const onDragEnter = useCallback(
		(e: DragEvent) => {
			if (!e.dataTransfer.types.includes(DRAG_MIME)) return;
			e.preventDefault();
			dragDepth.current++;
			if (dragDepth.current === 1) setDragHoverBar(position);
		},
		[position, setDragHoverBar],
	);

	const onDragLeave = useCallback(() => {
		dragDepth.current = Math.max(0, dragDepth.current - 1);
		if (
			dragDepth.current === 0 &&
			useAppStore.getState().dragHoverBar === position
		) {
			setTimeout(() => {
				if (
					dragDepth.current === 0 &&
					useAppStore.getState().dragHoverBar === position
				) {
					setDragHoverBar(null);
				}
			}, 0);
		}
	}, [position, setDragHoverBar]);

	const onDragOver = useCallback((e: DragEvent) => {
		if (e.dataTransfer.types.includes(DRAG_MIME)) {
			e.preventDefault();
			e.dataTransfer.dropEffect = "move";
		}
	}, []);

	return { isHovered, dragHandlers: { onDragEnter, onDragLeave, onDragOver } };
}
