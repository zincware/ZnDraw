import { Box } from "@mui/material";
import { useAppStore } from "../store";
import { PANELS } from "./registry";

const HEIGHT = 260;

export function BottomZone() {
	const active = useAppStore((s) => s.activeBottom);
	if (!active) return null;
	const def = PANELS[active];
	if (def.kind !== "tool") return null;
	const Component = def.component;

	return (
		<Box
			data-testid="bottom-zone"
			sx={{
				height: HEIGHT,
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
		</Box>
	);
}
