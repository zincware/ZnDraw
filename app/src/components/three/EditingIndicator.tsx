import { Box, Chip, Typography } from "@mui/material";
import EditIcon from "@mui/icons-material/Edit";
import OpenWithIcon from "@mui/icons-material/OpenWith";
import RotateLeftIcon from "@mui/icons-material/RotateLeft";
import AspectRatioIcon from "@mui/icons-material/AspectRatio";
import { useAppStore } from "../../store";

const transformModeIcons = {
	translate: <OpenWithIcon fontSize="small" />,
	rotate: <RotateLeftIcon fontSize="small" />,
	scale: <AspectRatioIcon fontSize="small" />,
};

const transformModeLabels = {
	translate: "Translate",
	rotate: "Rotate",
	scale: "Scale",
};

const axisColors = {
	x: "#f44336", // Red
	y: "#4caf50", // Green
	z: "#2196f3", // Blue
};

/**
 * EditingIndicator shows when the user is in geometry editing mode
 * Similar to DrawingIndicator but for transform controls editing
 */
export default function EditingIndicator() {
	const mode = useAppStore((state) => state.mode);
	const transformMode = useAppStore((state) => state.transformMode);
	const editingSelectedAxis = useAppStore((state) => state.editingSelectedAxis);

	// Only show if editing mode is active
	if (mode !== "editing") {
		return null;
	}

	return (
		<Box
			sx={{
				position: "fixed",
				bottom: 80,
				right: 16,
				zIndex: (theme) => theme.zIndex.tooltip,
				display: "flex",
				flexDirection: "column",
				alignItems: "flex-end",
				gap: 0.5,
			}}
		>
			{/* Axis indicator - shows when holding X/Y/Z */}
			{editingSelectedAxis && (
				<Chip
					label={
						<Typography variant="body2" sx={{ fontWeight: "bold" }}>
							{editingSelectedAxis.toUpperCase()} axis
						</Typography>
					}
					size="small"
					sx={{
						backgroundColor: axisColors[editingSelectedAxis],
						color: "white",
						fontWeight: "bold",
						animation: "pulse 0.5s ease-in-out",
						"@keyframes pulse": {
							"0%": { transform: "scale(1)" },
							"50%": { transform: "scale(1.05)" },
							"100%": { transform: "scale(1)" },
						},
					}}
				/>
			)}
			<Chip
				icon={<EditIcon />}
				label={
					<Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
						{transformModeIcons[transformMode]}
						<Typography variant="body2">
							{transformModeLabels[transformMode]}
						</Typography>
					</Box>
				}
				color="secondary"
				size="medium"
				sx={{
					"& .MuiChip-label": {
						px: 1,
					},
				}}
			/>
			<Typography
				variant="caption"
				sx={{ color: "text.secondary", fontSize: "0.7rem" }}
			>
				{editingSelectedAxis
					? "↑/↓ to adjust (Shift for larger steps)"
					: "T: cycle modes | X/Y/Z: select axis"}
			</Typography>
		</Box>
	);
}
