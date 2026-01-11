import { Box, Typography, Button } from "@mui/material";
import { useTheme } from "@mui/material/styles";
import type { InitializationError } from "../store";

interface CanvasErrorStateProps {
	error: InitializationError;
}

/**
 * Error state displayed when Canvas initialization fails.
 * Shows the error message and a refresh button.
 */
export function CanvasErrorState({ error }: CanvasErrorStateProps) {
	const theme = useTheme();

	return (
		<Box
			sx={{
				width: "100%",
				height: "calc(100vh - 64px)",
				display: "flex",
				flexDirection: "column",
				alignItems: "center",
				justifyContent: "center",
				bgcolor: theme.palette.background.default,
				gap: 2,
			}}
		>
			<Typography variant="h6" color="error">
				{error.message}
			</Typography>
			{error.details && (
				<Typography variant="body2" color="text.secondary">
					{error.details}
				</Typography>
			)}
			<Button
				variant="contained"
				onClick={() => window.location.reload()}
				sx={{ mt: 2 }}
			>
				Refresh Page
			</Button>
		</Box>
	);
}
