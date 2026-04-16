import { Box, Typography } from "@mui/material";
import type {
	DockviewApi,
	DockviewDidDropEvent,
	DockviewReadyEvent,
} from "dockview-react";
import { DockviewReact } from "dockview-react";
import "dockview-react/dist/styles/dockview.css";
import { useCallback, useRef, useState } from "react";
import { PlotView } from "./PlotView";
import { PANELS } from "./registry";
import { ViewerView } from "./ViewerView";

const components = {
	viewer: ViewerView,
	plotView: PlotView,
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

function onDidDrop(event: DockviewDidDropEvent) {
	const dt = (event.nativeEvent as DragEvent | undefined)?.dataTransfer;
	const key = dt?.getData(DRAG_MIME_PLOT);
	if (!key) return;
	const id = `plot-${key}`;
	if (event.api.getPanel(id)) {
		event.api.getPanel(id)?.api.setActive();
		return;
	}
	event.api.addPanel({
		id,
		component: "plotView",
		title: key,
		params: { figureKey: key },
		position: {
			referenceGroup: event.group ?? undefined,
			direction: "within",
		},
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
				onDidDrop={onDidDrop}
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
