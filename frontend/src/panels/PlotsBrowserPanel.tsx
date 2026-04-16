import CloseIcon from "@mui/icons-material/Close";
import {
	Box,
	Divider,
	IconButton,
	List,
	ListItem,
	ListItemButton,
	ListItemText,
	Typography,
} from "@mui/material";
import { useEffect, useMemo, useState } from "react";
import { useFigureList } from "../hooks/useFigures";
import { getDockviewApi } from "./DockviewLayout";
import { closePlotTab, openPlotTab } from "./plotViewFactory";

const DRAG_MIME_PLOT = "application/x-zndraw-plot-key";

function useOpenPlotKeys(): Set<string> {
	const [version, setVersion] = useState(0);

	useEffect(() => {
		const api = getDockviewApi();
		if (!api) return;
		const add = api.onDidAddPanel(() => setVersion((v) => v + 1));
		const remove = api.onDidRemovePanel(() => setVersion((v) => v + 1));
		return () => {
			add.dispose();
			remove.dispose();
		};
	}, []);

	return useMemo(() => {
		const api = getDockviewApi();
		if (!api) return new Set<string>();
		return new Set(
			api.panels
				.filter((p) => p.id.startsWith("plot-"))
				.map((p) => p.id.slice("plot-".length)),
		);
		// version forces recomputation when dockview panels change
		// biome-ignore lint/correctness/useExhaustiveDependencies: version is a tick trigger
	}, [version]);
}

export function PlotsBrowserPanel() {
	const { data, isLoading } = useFigureList();
	const openKeys = useOpenPlotKeys();
	const allKeys = useMemo(() => data?.items ?? [], [data]);

	const onRowClick = (key: string) => {
		const api = getDockviewApi();
		if (!api) return;
		openPlotTab(api, key);
	};

	const onRowClose = (key: string) => {
		const api = getDockviewApi();
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

	const openKeysArr = allKeys.filter((k: string) => openKeys.has(k));

	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Typography variant="overline" sx={{ px: 2, pt: 1 }}>
				Available Plots
			</Typography>
			<List dense sx={{ flexGrow: 1, overflow: "auto" }}>
				{allKeys.map((k: string) => (
					<ListItem
						key={k}
						disablePadding
						draggable
						onDragStart={(e) => onDragStart(e, k)}
					>
						<ListItemButton onClick={() => onRowClick(k)}>
							<Box
								sx={{
									width: 8,
									height: 8,
									borderRadius: "50%",
									bgcolor: openKeys.has(k) ? "primary.main" : "transparent",
									border: 1,
									borderColor: "primary.main",
									mr: 1.5,
								}}
							/>
							<ListItemText primary={k} />
						</ListItemButton>
					</ListItem>
				))}
			</List>
			{openKeysArr.length > 0 && (
				<>
					<Divider />
					<Typography variant="overline" sx={{ px: 2, pt: 1 }}>
						Currently Open ({openKeysArr.length})
					</Typography>
					<List dense>
						{openKeysArr.map((k: string) => (
							<ListItem
								key={k}
								secondaryAction={
									<IconButton
										edge="end"
										size="small"
										onClick={() => onRowClose(k)}
									>
										<CloseIcon fontSize="small" />
									</IconButton>
								}
							>
								<ListItemText primary={k} />
							</ListItem>
						))}
					</List>
				</>
			)}
		</Box>
	);
}
