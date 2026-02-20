import { Addchart as AddChartIcon } from "@mui/icons-material";
import { IconButton } from "@mui/material";
import Tooltip from "@mui/material/Tooltip";
import { useMemo } from "react";
import { useShallow } from "zustand/react/shallow";
import { useFigureList } from "../hooks/useFigures";
import { useWindowManagerStore } from "../stores/windowManagerStore";

function AddPlotButton() {
	// 1. Get all figure keys available on the server
	const { data: figureListResponse } = useFigureList();
	const allServerKeys = figureListResponse?.items || [];

	// 2. Get all figure keys currently being displayed in open windows
	// Use a selector that only returns the figureKeys array to avoid unnecessary re-renders
	const displayedKeys = useWindowManagerStore(
		useShallow((state) =>
			Object.values(state.openWindows).map((win) => win.figureKey),
		),
	);

	// 3. Find the first server key that is not currently displayed
	const nextKeyToShow = useMemo(() => {
		const displayedSet = new Set(displayedKeys);
		return allServerKeys.find((key) => !displayedSet.has(key));
	}, [allServerKeys, displayedKeys]);

	const { openWindow } = useWindowManagerStore();

	const handleOpenNextPlot = () => {
		if (nextKeyToShow) {
			openWindow(nextKeyToShow);
		}
	};

	return (
		<Tooltip title={"Add a new plot window"}>
			<IconButton
				color="inherit"
				aria-label="add plot"
				onClick={handleOpenNextPlot}
				disabled={!nextKeyToShow} // Button is disabled if no available plots are left
			>
				<AddChartIcon />
			</IconButton>
		</Tooltip>
	);
}

export default AddPlotButton;
