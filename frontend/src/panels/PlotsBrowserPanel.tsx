import CloseIcon from "@mui/icons-material/Close";
import {
	Box,
	IconButton,
	List,
	ListItem,
	ListItemButton,
	ListItemText,
	Typography,
} from "@mui/material";
import { useEffect, useMemo, useState } from "react";
import { useFigureList } from "../hooks/useFigures";
import { useDockviewApi } from "../stores/dockviewApiStore";
import { closePlotTab, openPlotTab } from "./plotViewFactory";

const DRAG_MIME_PLOT = "application/x-zndraw-plot-key";

function useOpenPlotKeys(): Set<string> {
	const api = useDockviewApi((s) => s.api);
	const [version, setVersion] = useState(0);

	useEffect(() => {
		if (!api) return;
		const disposables = [
			api.onDidAddPanel(() => setVersion((v) => v + 1)),
			api.onDidRemovePanel(() => setVersion((v) => v + 1)),
		];
		return () => {
			for (const d of disposables) d.dispose();
		};
	}, [api]);

	// biome-ignore lint/correctness/useExhaustiveDependencies: version is an intentional tick trigger — it forces recompute when panels are added/removed
	return useMemo(() => {
		if (!api) return new Set<string>();
		return new Set(
			api.panels
				.filter((p) => p.id.startsWith("plot-"))
				.map((p) => p.id.slice("plot-".length)),
		);
	}, [api, version]);
}

export function PlotsBrowserPanel() {
	const api = useDockviewApi((s) => s.api);
	const { data, isLoading } = useFigureList();
	const openKeys = useOpenPlotKeys();
	const allKeys = useMemo(() => data?.items ?? [], [data]);

	const onRowClick = (key: string) => {
		if (!api) return;
		openPlotTab(api, key);
	};

	const onRowClose = (e: React.MouseEvent, key: string) => {
		e.stopPropagation();
		if (!api) return;
		closePlotTab(api, key);
	};

	const onDragStart = (e: React.DragEvent, key: string) => {
		e.dataTransfer.setData(DRAG_MIME_PLOT, key);
		e.dataTransfer.effectAllowed = "copy";
	};

	if (isLoading) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography>Loading plots…</Typography>
			</Box>
		);
	}

	if (allKeys.length === 0) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography color="text.secondary">No plots available</Typography>
			</Box>
		);
	}

	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Typography variant="overline" sx={{ px: 2, pt: 1 }}>
				Plots
			</Typography>
			<List dense sx={{ flexGrow: 1, overflow: "auto" }}>
				{allKeys.map((k: string) => {
					const isOpen = openKeys.has(k);
					return (
						<ListItem
							key={k}
							data-testid={`plot-row-${k}`}
							data-open={isOpen ? "true" : "false"}
							disablePadding
							draggable
							onDragStart={(e) => onDragStart(e, k)}
							sx={{
								"& .plot-row-close": { visibility: "hidden" },
								"&:hover .plot-row-close": {
									visibility: isOpen ? "visible" : "hidden",
								},
							}}
							secondaryAction={
								isOpen ? (
									<IconButton
										className="plot-row-close"
										edge="end"
										size="small"
										aria-label={`Close ${k}`}
										onClick={(e) => onRowClose(e, k)}
									>
										<CloseIcon fontSize="small" />
									</IconButton>
								) : null
							}
						>
							<ListItemButton onClick={() => onRowClick(k)}>
								<Box
									sx={{
										width: 8,
										height: 8,
										borderRadius: "50%",
										bgcolor: isOpen ? "primary.main" : "transparent",
										border: 1,
										borderColor: "primary.main",
										mr: 1.5,
									}}
								/>
								<ListItemText primary={k} />
							</ListItemButton>
						</ListItem>
					);
				})}
			</List>
		</Box>
	);
}
