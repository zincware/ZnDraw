import { lazy, Suspense } from "react";
import { useWindowManagerStore } from "../stores/windowManagerStore";
import { useShallow } from "zustand/react/shallow";
import { Box, CircularProgress } from "@mui/material";

// Lazy load FigureWindow - Plotly (4.8MB) only loads when a figure is opened
const FigureWindow = lazy(() => import("./FigureWindow"));

function WindowManager() {
	const windowIds = useWindowManagerStore(
		useShallow((state) => Object.keys(state.openWindows)),
	);

	if (windowIds.length === 0) return null;

	return (
		<Suspense
			fallback={
				<Box
					sx={{
						position: "fixed",
						top: "50%",
						left: "50%",
						transform: "translate(-50%, -50%)",
						zIndex: 9999,
					}}
				>
					<CircularProgress />
				</Box>
			}
		>
			{windowIds.map((windowId) => (
				<FigureWindow key={windowId} windowId={windowId} />
			))}
		</Suspense>
	);
}

export default WindowManager;
