import { Box, Typography } from "@mui/material";
import type { IDockviewPanelProps } from "dockview-react";

const Stub = ({ label }: { label: string }) => (
	<Box sx={{ p: 2 }}>
		<Typography variant="h6">{label}</Typography>
		<Typography variant="body2" color="text.secondary">
			Placeholder — implementation in a later task.
		</Typography>
	</Box>
);

export const StubSelections = () => <Stub label="Selections" />;
export const StubModifiers = () => <Stub label="Modifiers" />;
export const StubAnalysis = () => <Stub label="Analysis" />;
export const StubGeometries = () => <Stub label="Geometries" />;
export const StubPlotsBrowser = () => <Stub label="Plots Browser" />;
export const StubRooms = () => <Stub label="Rooms" />;
export const StubFilesystem = () => <Stub label="Files" />;
export const StubChat = () => <Stub label="Chat" />;

export const StubViewerView = (_props: IDockviewPanelProps) => (
	<Stub label="3D Viewer (stub)" />
);
