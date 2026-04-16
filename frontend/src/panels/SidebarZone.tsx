import { Box } from "@mui/material";
import { useAppStore } from "../store";
import type { BarPosition } from "./registry";
import { PANELS } from "./registry";

interface SidebarZoneProps {
	position: Exclude<BarPosition, "bottom">;
}

const WIDTH = 320;

export function SidebarZone({ position }: SidebarZoneProps) {
	const active = useAppStore((s) =>
		position === "left" ? s.activeLeft : s.activeRight,
	);
	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	return (
		<Box
			data-testid={`sidebar-zone-${position}`}
			sx={{
				width: WIDTH,
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
		</Box>
	);
}
