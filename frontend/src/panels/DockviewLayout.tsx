import { Box, Typography } from "@mui/material";
import type { DockviewApi, DockviewReadyEvent } from "dockview-react";
import { DockviewReact } from "dockview-react";
import "dockview-react/dist/styles/dockview.css";
import { useCallback, useRef, useState } from "react";
import { PANELS } from "./registry";
import { ViewerView } from "./ViewerView";

const components = {
	viewer: ViewerView,
	// plotView registered in Task 7
};

const DRAG_MIME_PLOT = "application/x-zndraw-plot-key";

let sharedApi: DockviewApi | null = null;

export function getDockviewApi(): DockviewApi | null {
	return sharedApi;
}

function addViewerPanel(api: DockviewApi) {
	api.addPanel({
		id: "viewer",
		component: "viewer",
		title: PANELS.viewer.kind === "view" ? PANELS.viewer.title : "Viewer",
	});
}

export function DockviewLayout() {
	const apiRef = useRef<DockviewApi | null>(null);
	const [isEmpty, setIsEmpty] = useState(false);

	const onReady = useCallback((event: DockviewReadyEvent) => {
		apiRef.current = event.api;
		sharedApi = event.api;

		addViewerPanel(event.api);
		setIsEmpty(event.api.panels.length === 0);

		event.api.onUnhandledDragOverEvent((e) => {
			const dt = (e.nativeEvent as DragEvent | undefined)?.dataTransfer;
			if (dt?.types.includes(DRAG_MIME_PLOT)) e.accept();
		});

		event.api.onDidRemovePanel(() => {
			setIsEmpty(event.api.panels.length === 0);
		});
		event.api.onDidAddPanel(() => {
			setIsEmpty(event.api.panels.length === 0);
		});
	}, []);

	return (
		<Box sx={{ flexGrow: 1, position: "relative", minWidth: 0, minHeight: 0 }}>
			<DockviewReact
				className="dockview-theme-light"
				onReady={onReady}
				components={components}
				floatingGroupBounds="boundedWithinViewport"
			/>
			{isEmpty && (
				<Box
					data-testid="dockview-welcome"
					sx={{
						position: "absolute",
						inset: 0,
						display: "flex",
						alignItems: "center",
						justifyContent: "center",
						pointerEvents: "none",
					}}
				>
					<Typography variant="body1" color="text.secondary">
						No room selected. Open the Rooms panel to pick one.
					</Typography>
				</Box>
			)}
		</Box>
	);
}
