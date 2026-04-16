import Analytics from "@mui/icons-material/Analytics";
import Build from "@mui/icons-material/Build";
import Category from "@mui/icons-material/Category";
import Chat from "@mui/icons-material/Chat";
import FilterCenterFocus from "@mui/icons-material/FilterCenterFocus";
import Folder from "@mui/icons-material/Folder";
import MeetingRoom from "@mui/icons-material/MeetingRoom";
import ShowChart from "@mui/icons-material/ShowChart";
import type { SvgIconComponent } from "@mui/icons-material";
import type { IDockviewPanelProps } from "dockview-react";
import type { ComponentType } from "react";
import { FilesystemPanel } from "./FilesystemPanel";
import { PlotsBrowserPanel } from "./PlotsBrowserPanel";
import { RoomsPanel } from "./RoomsPanel";
import {
	StubAnalysis,
	StubChat,
	StubGeometries,
	StubModifiers,
	StubSelections,
	StubViewerView,
} from "./stubs";

export type PanelId =
	| "selections"
	| "modifiers"
	| "analysis"
	| "geometries"
	| "plots-browser"
	| "rooms"
	| "filesystem"
	| "chat"
	| "viewer";

export type BarPosition = "left" | "right" | "bottom";

type PanelDefault = {
	bar: BarPosition | "editor";
	active?: boolean;
	order?: number;
};

type ToolPanelDef = {
	kind: "tool";
	icon: SvgIconComponent;
	label: string;
	component: ComponentType;
	default: PanelDefault;
};

type ViewPanelDef = {
	kind: "view";
	title: string;
	component: ComponentType<IDockviewPanelProps>;
	closable: boolean;
	cascadeOnClose?: (string | RegExp)[];
	onClose?: () => void;
};

export type PanelDef = ToolPanelDef | ViewPanelDef;

export const PANELS: Record<PanelId, PanelDef> = {
	selections: {
		kind: "tool",
		icon: FilterCenterFocus,
		label: "Selections",
		component: StubSelections,
		default: { bar: "left", order: 0 },
	},
	modifiers: {
		kind: "tool",
		icon: Build,
		label: "Modifiers",
		component: StubModifiers,
		default: { bar: "left", order: 1 },
	},
	analysis: {
		kind: "tool",
		icon: Analytics,
		label: "Analysis",
		component: StubAnalysis,
		default: { bar: "left", order: 2 },
	},
	geometries: {
		kind: "tool",
		icon: Category,
		label: "Geometries",
		component: StubGeometries,
		default: { bar: "left", order: 3 },
	},
	"plots-browser": {
		kind: "tool",
		icon: ShowChart,
		label: "Plots",
		component: PlotsBrowserPanel,
		default: { bar: "left", order: 4 },
	},
	rooms: {
		kind: "tool",
		icon: MeetingRoom,
		label: "Rooms",
		component: RoomsPanel,
		default: { bar: "left", order: 5 },
	},
	filesystem: {
		kind: "tool",
		icon: Folder,
		label: "Files",
		component: FilesystemPanel,
		default: { bar: "left", order: 6 },
	},
	chat: {
		kind: "tool",
		icon: Chat,
		label: "Chat",
		component: StubChat,
		default: { bar: "right", order: 0 },
	},
	viewer: {
		kind: "view",
		title: "3D Viewer",
		component: StubViewerView,
		closable: true,
		cascadeOnClose: [/^plot-.+$/],
		// onClose is assigned in a later task (requires useLeaveRoom hook)
	},
};

export const TOOL_PANEL_IDS = (Object.entries(PANELS) as [PanelId, PanelDef][])
	.filter(([, def]) => def.kind === "tool")
	.map(([id]) => id);

export function getDefaultsForBar(bar: BarPosition): PanelId[] {
	return (Object.entries(PANELS) as [PanelId, PanelDef][])
		.filter(([, def]) => def.kind === "tool" && def.default.bar === bar)
		.sort(
			([, a], [, b]) =>
				(a.kind === "tool" ? (a.default.order ?? 0) : 0) -
				(b.kind === "tool" ? (b.default.order ?? 0) : 0),
		)
		.map(([id]) => id);
}
