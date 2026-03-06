import { Box } from "@mui/material";
import { useCallback, useEffect, useRef, useState } from "react";
import { useAppStore } from "../store";

interface FrameSelectionInputProps {
	frameCount: number;
	containerWidth: number;
	onPreviewChange: (range: [number, number] | null) => void;
}

/**
 * Transparent overlay that captures Shift+Click/Drag for frame range selection.
 *
 * - Shift+Click: select range from current frame to clicked frame
 * - Shift+Drag: paint a contiguous range
 * - Ctrl/Cmd+Shift: additive (union with existing selection)
 *
 * Uses window-level listeners during drag so the cursor can leave the overlay.
 */
const FrameSelectionInput: React.FC<FrameSelectionInputProps> = ({
	frameCount,
	containerWidth,
	onPreviewChange,
}) => {
	const [shiftHeld, setShiftHeld] = useState(false);
	const dragStartFrame = useRef<number | null>(null);
	const overlayRef = useRef<HTMLDivElement>(null);

	const positionToFrame = useCallback(
		(clientX: number): number => {
			if (!overlayRef.current || frameCount <= 1) return 0;
			const rect = overlayRef.current.getBoundingClientRect();
			const x = Math.max(0, Math.min(clientX - rect.left, rect.width));
			return Math.round((x / rect.width) * (frameCount - 1));
		},
		[frameCount],
	);

	// Track Shift key state globally
	useEffect(() => {
		const onKeyDown = (e: KeyboardEvent) => {
			if (e.key === "Shift") setShiftHeld(true);
		};
		const onKeyUp = (e: KeyboardEvent) => {
			if (e.key === "Shift") {
				setShiftHeld(false);
				// Cancel any in-progress drag
				if (dragStartFrame.current !== null) {
					dragStartFrame.current = null;
					onPreviewChange(null);
				}
			}
		};
		window.addEventListener("keydown", onKeyDown);
		window.addEventListener("keyup", onKeyUp);
		return () => {
			window.removeEventListener("keydown", onKeyDown);
			window.removeEventListener("keyup", onKeyUp);
		};
	}, [onPreviewChange]);

	const commitSelection = useCallback(
		(endFrame: number, startFrame: number, ctrlKey: boolean) => {
			const currentFrame = useAppStore.getState().currentFrame;
			const frame_selection = useAppStore.getState().frame_selection;

			const wasDrag = startFrame !== endFrame;
			const rangeStart = wasDrag ? startFrame : currentFrame;

			const start = Math.min(rangeStart, endFrame);
			const end = Math.max(rangeStart, endFrame);
			let newRange = Array.from(
				{ length: end - start + 1 },
				(_, i) => start + i,
			);

			// Ctrl/Cmd+Shift = additive (union)
			if (ctrlKey && frame_selection) {
				const merged = new Set([...frame_selection, ...newRange]);
				newRange = [...merged].sort((a, b) => a - b);
			}

			useAppStore.getState().updateFrameSelection(newRange);
		},
		[],
	);

	const handleMouseDown = useCallback(
		(e: React.MouseEvent) => {
			if (!shiftHeld) return;
			e.preventDefault();
			e.stopPropagation();

			const startFrame = positionToFrame(e.clientX);
			dragStartFrame.current = startFrame;
			onPreviewChange([startFrame, startFrame]);

			const onWindowMouseMove = (ev: MouseEvent) => {
				if (dragStartFrame.current === null) return;
				const frame = positionToFrame(ev.clientX);
				onPreviewChange([
					Math.min(dragStartFrame.current, frame),
					Math.max(dragStartFrame.current, frame),
				]);
			};

			const onWindowMouseUp = (ev: MouseEvent) => {
				window.removeEventListener("mousemove", onWindowMouseMove);
				window.removeEventListener("mouseup", onWindowMouseUp);

				if (dragStartFrame.current === null) return;
				const endFrame = positionToFrame(ev.clientX);
				commitSelection(
					endFrame,
					dragStartFrame.current,
					ev.ctrlKey || ev.metaKey,
				);
				onPreviewChange(null);
				dragStartFrame.current = null;
			};

			window.addEventListener("mousemove", onWindowMouseMove);
			window.addEventListener("mouseup", onWindowMouseUp);
		},
		[shiftHeld, positionToFrame, onPreviewChange, commitSelection],
	);

	if (containerWidth === 0 || frameCount <= 1) return null;

	return (
		<Box
			ref={overlayRef}
			onMouseDown={handleMouseDown}
			sx={{
				position: "absolute",
				top: 0,
				left: 0,
				right: 0,
				bottom: 0,
				zIndex: 1,
				cursor: shiftHeld ? "crosshair" : "default",
				pointerEvents: shiftHeld ? "auto" : "none",
			}}
		/>
	);
};

export default FrameSelectionInput;
