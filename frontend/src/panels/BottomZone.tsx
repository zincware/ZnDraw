import { Box } from "@mui/material";
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { BOTTOM_MAX_PX, BOTTOM_MIN_PX, PANELS } from "./registry";

export function BottomZone() {
	const active = useAppStore((s) => s.activeBottom);
	const height = useAppStore((s) => s.bottomHeight);
	const setBarSize = useAppStore((s) => s.setBarSize);
	const zoneRef = useRef<HTMLDivElement | null>(null);

	const onPointerDown = useCallback(
		(e: React.PointerEvent<HTMLDivElement>) => {
			e.preventDefault();
			(e.target as HTMLElement).setPointerCapture(e.pointerId);
			const startY = e.clientY;
			const startHeight =
				zoneRef.current?.getBoundingClientRect().height ?? height;

			const onMove = (ev: PointerEvent) => {
				const delta = ev.clientY - startY;
				setBarSize("bottom", startHeight - delta); // dragging up grows the zone
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

	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	const safeHeight = Math.min(BOTTOM_MAX_PX, Math.max(BOTTOM_MIN_PX, height));

	return (
		<Box
			ref={zoneRef}
			data-testid="bottom-zone"
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
					zIndex: 1,
					"&:hover": { bgcolor: "action.hover" },
					touchAction: "none",
				}}
			/>
		</Box>
	);
}
