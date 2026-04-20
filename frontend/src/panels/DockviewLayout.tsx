import { Box, Typography } from "@mui/material";
import { useColorScheme } from "@mui/material/styles";
import type {
	DockviewApi,
	DockviewDidDropEvent,
	DockviewReadyEvent,
} from "dockview-react";
import { DockviewReact, themeDark, themeLight } from "dockview-react";
import "dockview-react/dist/styles/dockview.css";
import "./dockview-mui.css";
import { useCallback, useEffect, useRef, useState } from "react";
import { useDockviewApi } from "../stores/dockviewApiStore";
import { GroupActions } from "./groupActions";
import { PlotView } from "./PlotView";
import { openPlotTab } from "./plotViewFactory";
import { PANELS } from "./registry";
import { ViewerView } from "./ViewerView";

const components = {
	viewer: ViewerView,
	plotView: PlotView,
};

const DRAG_MIME_PLOT = "application/x-zndraw-plot-key";

function addViewerPanel(api: DockviewApi) {
	api.addPanel({
		id: "viewer",
		component: "viewer",
		title: PANELS.viewer.kind === "view" ? PANELS.viewer.title : "Viewer",
	});
}

export function resetDockview(api: DockviewApi): void {
	for (const p of api.panels) p.api.close();
	addViewerPanel(api);
}

export function ensureViewerPanel(api: DockviewApi): void {
	if (!api.getPanel("viewer")) addViewerPanel(api);
}

function onDidDrop(event: DockviewDidDropEvent) {
	const dt = (event.nativeEvent as DragEvent | undefined)?.dataTransfer;
	const key = dt?.getData(DRAG_MIME_PLOT);
	if (!key) return;
	openPlotTab(event.api, key, {
		referenceGroupId: event.group?.id,
		direction: "within",
	});
}

export function DockviewLayout() {
	const apiRef = useRef<DockviewApi | null>(null);
	const disposablesRef = useRef<Array<{ dispose(): void }>>([]);
	const [isEmpty, setIsEmpty] = useState(false);
	const { mode, systemMode } = useColorScheme();
	const resolvedMode = mode === "system" ? systemMode : mode;
	const dockTheme = resolvedMode === "dark" ? themeDark : themeLight;

	const onReady = useCallback((event: DockviewReadyEvent) => {
		apiRef.current = event.api;
		useDockviewApi.getState().setApi(event.api);

		addViewerPanel(event.api);
		setIsEmpty(event.api.panels.length === 0);

		disposablesRef.current.push(
			event.api.onUnhandledDragOverEvent((e) => {
				const dt = (e.nativeEvent as DragEvent | undefined)?.dataTransfer;
				if (dt?.types.includes(DRAG_MIME_PLOT)) e.accept();
			}),
			event.api.onDidRemovePanel(() => {
				setIsEmpty(event.api.panels.length === 0);
			}),
			event.api.onDidAddPanel(() => {
				setIsEmpty(event.api.panels.length === 0);
			}),
		);
	}, []);

	useEffect(() => {
		return () => {
			for (const d of disposablesRef.current) d.dispose();
			disposablesRef.current = [];
			useDockviewApi.getState().setApi(null);
		};
	}, []);

	return (
		<Box sx={{ flexGrow: 1, position: "relative", minWidth: 0, minHeight: 0 }}>
			<Box sx={{ position: "absolute", inset: 0 }}>
				<DockviewReact
					theme={dockTheme}
					onReady={onReady}
					components={components}
					onDidDrop={onDidDrop}
					floatingGroupBounds="boundedWithinViewport"
					rightHeaderActionsComponent={GroupActions}
				/>
			</Box>
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
