import { Box, keyframes } from "@mui/material";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import {
	BOTTOM_MAX_PX,
	BOTTOM_MIN_PX,
	PANELS,
	type PanelId,
} from "./registry";

const DRAG_MIME = "application/x-zndraw-panel-id";

const shimmer = keyframes`
	0%, 100% { background-color: rgba(25, 118, 210, 0.12); }
	50% { background-color: rgba(25, 118, 210, 0.28); }
`;

export function BottomZone() {
	const active = useAppStore((s) => s.activeBottom);
	const height = useAppStore((s) => s.bottomHeight);
	const setBarSize = useAppStore((s) => s.setBarSize);
	const isDragActive = useAppStore((s) => s.isPanelDragActive);
	const hoverBar = useAppStore((s) => s.dragHoverBar);
	const setHoverBar = useAppStore((s) => s.setDragHoverBar);
	const dropIconOnPanel = useAppStore((s) => s.dropIconOnPanel);
	const setPanelDragActive = useAppStore((s) => s.setPanelDragActive);
	const zoneRef = useRef<HTMLDivElement | null>(null);
	const dragDepth = useRef(0);

	const isHovered = hoverBar === "bottom";

	const onPointerDown = useCallback(
		(e: React.PointerEvent<HTMLDivElement>) => {
			e.preventDefault();
			(e.target as HTMLElement).setPointerCapture(e.pointerId);
			const startY = e.clientY;
			const startHeight =
				zoneRef.current?.getBoundingClientRect().height ?? height;

			const onMove = (ev: PointerEvent) => {
				const delta = ev.clientY - startY;
				setBarSize("bottom", startHeight - delta);
			};
			const onUp = (ev: PointerEvent) => {
				(e.target as HTMLElement).releasePointerCapture(ev.pointerId);
				window.removeEventListener("pointermove", onMove);
				window.removeEventListener("pointerup", onUp);
				window.removeEventListener("pointercancel", onUp);
			};
			window.addEventListener("pointermove", onMove);
			window.addEventListener("pointerup", onUp);
			window.addEventListener("pointercancel", onUp);
		},
		[setBarSize, height],
	);

	const onDragOver = useCallback((e: React.DragEvent) => {
		if (e.dataTransfer.types.includes(DRAG_MIME)) {
			e.preventDefault();
			e.dataTransfer.dropEffect = "move";
		}
	}, []);

	const onDragEnter = useCallback(
		(e: React.DragEvent) => {
			if (!e.dataTransfer.types.includes(DRAG_MIME)) return;
			dragDepth.current++;
			if (dragDepth.current === 1) setHoverBar("bottom");
		},
		[setHoverBar],
	);

	const onDragLeave = useCallback(() => {
		dragDepth.current = Math.max(0, dragDepth.current - 1);
		if (
			dragDepth.current === 0 &&
			useAppStore.getState().dragHoverBar === "bottom"
		) {
			setTimeout(() => {
				if (
					dragDepth.current === 0 &&
					useAppStore.getState().dragHoverBar === "bottom"
				) {
					setHoverBar(null);
				}
			}, 0);
		}
	}, [setHoverBar]);

	const onDrop = useCallback(
		(e: React.DragEvent) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			dragDepth.current = 0;
			setPanelDragActive(false);
			if (!id) return;
			e.preventDefault();
			dropIconOnPanel(id as PanelId, "bottom");
		},
		[dropIconOnPanel, setPanelDragActive],
	);

	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	const safeHeight = Math.min(BOTTOM_MAX_PX, Math.max(BOTTOM_MIN_PX, height));

	return (
		<Box
			ref={zoneRef}
			data-testid="bottom-zone"
			data-drop-hover={isHovered}
			onDragEnter={onDragEnter}
			onDragOver={onDragOver}
			onDragLeave={onDragLeave}
			onDrop={onDrop}
			sx={{
				position: "relative",
				height: safeHeight,
				flexShrink: 0,
				borderTop: 1,
				borderColor: "divider",
				display: "flex",
				flexDirection: "column",
				overflow: "hidden",
				bgcolor: "background.default",
			}}
		>
			<Component />
			{isDragActive && (
				<Box
					data-testid="bottom-drop-overlay"
					sx={{
						position: "absolute",
						inset: 0,
						pointerEvents: "none",
						zIndex: 2,
						border: 2,
						borderStyle: "solid",
						borderColor: "primary.main",
						transition: "background-color 120ms ease",
						...(isHovered
							? { bgcolor: "rgba(25, 118, 210, 0.32)" }
							: { animation: `${shimmer} 1.2s ease-in-out infinite` }),
					}}
				/>
			)}
			<Box
				data-testid="bottom-resize"
				onPointerDown={onPointerDown}
				sx={{
					position: "absolute",
					left: 0,
					right: 0,
					top: -4,
					height: 8,
					cursor: "row-resize",
					zIndex: 3,
					"&:hover": { bgcolor: "action.hover" },
					touchAction: "none",
				}}
			/>
		</Box>
	);
}
