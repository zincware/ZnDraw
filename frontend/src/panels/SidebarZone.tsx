import { Box } from "@mui/material";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { shimmer } from "./dragStyles";
import {
	type BarPosition,
	PANELS,
	type PanelId,
	SIDEBAR_MAX_PX,
	SIDEBAR_MIN_PX,
} from "./registry";
import { useDragHover } from "./useDragHover";

const DRAG_MIME = "application/x-zndraw-panel-id";

interface SidebarZoneProps {
	position: Exclude<BarPosition, "bottom">;
}

export function SidebarZone({ position }: SidebarZoneProps) {
	const active = useAppStore((s) =>
		position === "left" ? s.activeLeft : s.activeRight,
	);
	const width = useAppStore((s) =>
		position === "left" ? s.leftWidth : s.rightWidth,
	);
	const setBarSize = useAppStore((s) => s.setBarSize);
	const isDragActive = useAppStore((s) => s.isPanelDragActive);
	const dropIconOnPanel = useAppStore((s) => s.dropIconOnPanel);
	const setPanelDragActive = useAppStore((s) => s.setPanelDragActive);
	const zoneRef = useRef<HTMLDivElement | null>(null);

	const { isHovered, dragHandlers } = useDragHover(position);

	const onPointerDown = useCallback(
		(e: React.PointerEvent<HTMLDivElement>) => {
			e.preventDefault();
			(e.target as HTMLElement).setPointerCapture(e.pointerId);
			const startX = e.clientX;
			const startWidth =
				zoneRef.current?.getBoundingClientRect().width ?? width;

			const onMove = (ev: PointerEvent) => {
				const delta = ev.clientX - startX;
				const next =
					position === "left" ? startWidth + delta : startWidth - delta;
				setBarSize(position, next);
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
		[position, setBarSize, width],
	);

	const onDrop = useCallback(
		(e: React.DragEvent) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			setPanelDragActive(false);
			if (!id) return;
			e.preventDefault();
			// Drop on the panel → move AND open the dropped icon's panel.
			dropIconOnPanel(id as PanelId, position);
		},
		[dropIconOnPanel, position, setPanelDragActive],
	);

	// Panel is only visible when something is active. Drag doesn't force-open
	// an empty panel.
	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	const safeWidth = Math.min(SIDEBAR_MAX_PX, Math.max(SIDEBAR_MIN_PX, width));

	return (
		<Box
			ref={zoneRef}
			data-testid={`sidebar-zone-${position}`}
			data-drop-hover={isHovered}
			{...dragHandlers}
			onDrop={onDrop}
			sx={{
				position: "relative",
				width: safeWidth,
				flexShrink: 0,
				borderColor: "divider",
				borderRight: position === "left" ? 1 : 0,
				borderLeft: position === "right" ? 1 : 0,
				display: "flex",
				flexDirection: "column",
				overflow: "hidden",
				bgcolor: "background.default",
			}}
		>
			<Component />
			{isDragActive && (
				<Box
					data-testid={`sidebar-drop-overlay-${position}`}
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
				data-testid={`sidebar-resize-${position}`}
				onPointerDown={onPointerDown}
				sx={{
					position: "absolute",
					top: 0,
					bottom: 0,
					width: 8,
					[position === "left" ? "right" : "left"]: -4,
					cursor: "col-resize",
					zIndex: 3,
					"&:hover": { bgcolor: "action.hover" },
					touchAction: "none",
				}}
			/>
		</Box>
	);
}
