import { Box, Paper, Typography, useTheme } from "@mui/material";
import { usePropertyInspectorSettings } from "../../hooks/usePropertyInspectorSettings";
import { useAppStore } from "../../store";
import { formatPropertyValue } from "../../utils/propertyFormatting";

/**
 * StaticInfoBox - Displays general scene information in the top-right corner
 * of the viewer.
 *
 * Shows:
 * - Number of particles in the current frame
 * - Current frame number
 * - Total frame count
 * - Selection count (if any particles are selected)
 */
export default function StaticInfoBox() {
	const theme = useTheme();
	const geometries = useAppStore((state) => state.geometries);
	const particleCount = useAppStore((state) => state.particleCount);
	const selections = useAppStore((state) => state.selections);
	const playing = useAppStore((state) => state.playing);
	const fps = useAppStore((state) => state.fps);
	const frameLoadTime = useAppStore((state) => state.frameLoadTime);

	// Get active state from PropertyInspector geometry
	const isActive = geometries["property-inspector"]?.data?.active ?? false;

	// Performance optimization: only fetch global properties when boxes are visible
	const { enabledProperties, propertyValues, isEnabled } =
		usePropertyInspectorSettings({
			category: "global",
			enabled: isActive,
		});

	if (!isActive) return null;

	const selectionCount = selections?.particles?.length || 0;

	return (
		<Box
			sx={{
				position: "absolute",
				top: 12,
				right: 12,
				width: 240,
				zIndex: 1300,
				pointerEvents: "auto",
			}}
		>
			<Paper
				elevation={3}
				sx={{
					width: "100%",
					display: "flex",
					flexDirection: "column",
					backgroundColor:
						theme.palette.mode === "dark"
							? "rgba(0, 0, 0, 0.8)"
							: "rgba(255, 255, 255, 0.9)",
					backdropFilter: "blur(10px)",
					border: `1px solid ${
						theme.palette.mode === "dark"
							? "rgba(255, 255, 255, 0.1)"
							: "rgba(0, 0, 0, 0.1)"
					}`,
				}}
			>
				<Box
					sx={{
						padding: 2,
					}}
				>
					<Typography
						variant="h6"
						sx={{
							color: theme.palette.text.primary,
							fontWeight: "bold",
							marginBottom: 1,
							fontSize: "0.9rem",
						}}
					>
						Scene Info
					</Typography>
					<Box sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}>
						<Typography
							variant="body2"
							sx={{ color: theme.palette.text.secondary }}
						>
							Particles:{" "}
							<strong style={{ color: theme.palette.text.primary }}>
								{particleCount}
							</strong>
						</Typography>
						{selectionCount > 0 && (
							<Typography
								variant="body2"
								sx={{ color: theme.palette.warning.main }}
							>
								Selected: <strong>{selectionCount}</strong>
							</Typography>
						)}
						{playing && fps !== null && (
							<Typography
								variant="body2"
								sx={{ color: theme.palette.success.main }}
							>
								FPS: <strong>{fps.toFixed(1)}</strong>
							</Typography>
						)}
						{!playing && frameLoadTime !== null && (
							<Typography
								variant="body2"
								sx={{ color: theme.palette.info.main }}
							>
								Load: <strong>{frameLoadTime}ms</strong>
							</Typography>
						)}
						{isEnabled &&
							propertyValues.map((query, index) => {
								const key = enabledProperties[index];

								// Skip if loading or error
								if (query.isLoading || query.isError || !query.data?.value) {
									return null;
								}

								const value = query.data.value;

								// Format value using utility function (no particleId for global props)
								const displayValue = formatPropertyValue(value);

								return (
									<Typography
										key={key}
										variant="caption"
										sx={{
											color: theme.palette.text.secondary,
											fontFamily: "monospace",
											fontSize: "0.75rem",
										}}
									>
										{key}: {displayValue}
									</Typography>
								);
							})}
					</Box>
				</Box>
			</Paper>
		</Box>
	);
}
