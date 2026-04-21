import type { SvgIconComponent } from "@mui/icons-material";
import Analytics from "@mui/icons-material/Analytics";
import Build from "@mui/icons-material/Build";
import Category from "@mui/icons-material/Category";
import Chat from "@mui/icons-material/Chat";
import FilterCenterFocus from "@mui/icons-material/FilterCenterFocus";
import Folder from "@mui/icons-material/Folder";
import MeetingRoom from "@mui/icons-material/MeetingRoom";
import ShowChart from "@mui/icons-material/ShowChart";
import type { IDockviewPanelProps } from "dockview-react";
import type { ComponentType } from "react";
import { GeometryPanel } from "../components/geometry/GeometryPanel";
import { SecondaryPanel } from "../components/SecondaryPanel";
import { SelectionsPanel } from "../components/SelectionsPanel";
import ChatPanel from "./ChatPanel";
import { FilesystemPanel } from "./FilesystemPanel";
import { PlotsBrowserPanel } from "./PlotsBrowserPanel";
import { RoomsPanel } from "./RoomsPanel";
import { ViewerView } from "./ViewerView";

const ModifiersPanelWrapper = () => <SecondaryPanel panelTitle="modifiers" />;
const AnalysisPanelWrapper = () => <SecondaryPanel panelTitle="analysis" />;

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
		component: SelectionsPanel,
		default: { bar: "left", order: 0 },
	},
	modifiers: {
		kind: "tool",
		icon: Build,
		label: "Modifiers",
		component: ModifiersPanelWrapper,
		default: { bar: "left", order: 1 },
	},
	analysis: {
		kind: "tool",
		icon: Analytics,
		label: "Analysis",
		component: AnalysisPanelWrapper,
		default: { bar: "left", order: 2 },
	},
	geometries: {
		kind: "tool",
		icon: Category,
		label: "Geometries",
		component: GeometryPanel,
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
		component: ChatPanel,
		default: { bar: "left", order: 7 },
	},
	viewer: {
		kind: "view",
		title: "3D Viewer",
		component: ViewerView,
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

// Sidebar / bottom-zone sizes (px). Session-only — stored in activityBarSlice.
export const SIDEBAR_DEFAULT_PX = 600;
export const SIDEBAR_MIN_PX = 200;
export const SIDEBAR_MAX_PX = 900;

export const BOTTOM_DEFAULT_PX = 260;
export const BOTTOM_MIN_PX = 120;
export const BOTTOM_MAX_PX = 560;
