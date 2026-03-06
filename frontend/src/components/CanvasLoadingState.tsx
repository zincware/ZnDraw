import { Box, CircularProgress } from "@mui/material";
import { useTheme } from "@mui/material/styles";

/**
 * Loading state displayed while the Canvas is initializing.
 * Shows a centered spinner.
 */
export function CanvasLoadingState() {
	const theme = useTheme();

	return (
		<Box
			sx={{
				width: "100%",
				height: "calc(100vh - 64px)",
				display: "flex",
				alignItems: "center",
				justifyContent: "center",
				bgcolor: theme.palette.background.default,
			}}
		>
			<CircularProgress />
		</Box>
	);
}
