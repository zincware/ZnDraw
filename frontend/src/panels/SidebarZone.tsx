import { Box } from "@mui/material";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import {
	type BarPosition,
	PANELS,
	SIDEBAR_MAX_PX,
	SIDEBAR_MIN_PX,
} from "./registry";

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
	const zoneRef = useRef<HTMLDivElement | null>(null);

	const onPointerDown = useCallback(
		(e: React.PointerEvent<HTMLDivElement>) => {
			e.preventDefault();
			(e.target as HTMLElement).setPointerCapture(e.pointerId);
			const startX = e.clientX;
			const startWidth = zoneRef.current?.getBoundingClientRect().width ?? width;

			const onMove = (ev: PointerEvent) => {
				const delta = ev.clientX - startX;
				const next = position === "left" ? startWidth + delta : startWidth - delta;
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

	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	// Clamp on render too, in case SSR-ish state had stale value.
	const safeWidth = Math.min(SIDEBAR_MAX_PX, Math.max(SIDEBAR_MIN_PX, width));

	return (
		<Box
			ref={zoneRef}
			data-testid={`sidebar-zone-${position}`}
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
			<Box
				data-testid={`sidebar-resize-${position}`}
				onPointerDown={onPointerDown}
				sx={{
					position: "absolute",
					top: 0,
					bottom: 0,
					width: 8,
					// inner edge: right for left sidebar, left for right sidebar.
					[position === "left" ? "right" : "left"]: -4,
					cursor: "col-resize",
					zIndex: 1,
					"&:hover": { bgcolor: "action.hover" },
					touchAction: "none",
				}}
			/>
		</Box>
	);
}
