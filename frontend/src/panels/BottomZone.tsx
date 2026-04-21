import { Box } from "@mui/material";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { shimmer } from "./dragStyles";
import { BOTTOM_MAX_PX, BOTTOM_MIN_PX, PANELS, type PanelId } from "./registry";
import { useDragHover } from "./useDragHover";

const DRAG_MIME = "application/x-zndraw-panel-id";

export function BottomZone() {
	const active = useAppStore((s) => s.activeBottom);
	const height = useAppStore((s) => s.bottomHeight);
	const setBarSize = useAppStore((s) => s.setBarSize);
	const isDragActive = useAppStore((s) => s.isPanelDragActive);
	const dropIconOnPanel = useAppStore((s) => s.dropIconOnPanel);
	const setPanelDragActive = useAppStore((s) => s.setPanelDragActive);
	const zoneRef = useRef<HTMLDivElement | null>(null);

	const { isHovered, dragHandlers } = useDragHover("bottom");

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

	const onDrop = useCallback(
		(e: React.DragEvent) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
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
			{...dragHandlers}
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
