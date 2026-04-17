import type { StateCreator } from "zustand";
import type { AppState } from "../../store";
import {
	type BarPosition,
	BOTTOM_DEFAULT_PX,
	BOTTOM_MAX_PX,
	BOTTOM_MIN_PX,
	getDefaultsForBar,
	type PanelId,
	PANELS,
	SIDEBAR_DEFAULT_PX,
	SIDEBAR_MAX_PX,
	SIDEBAR_MIN_PX,
} from "../../panels/registry";

export interface ActivityBarSlice {
	leftBarIcons: PanelId[];
	rightBarIcons: PanelId[];
	bottomBarIcons: PanelId[];
	activeLeft: PanelId | null;
	activeRight: PanelId | null;
	activeBottom: PanelId | null;

	leftWidth: number;
	rightWidth: number;
	bottomHeight: number;

	isPanelDragActive: boolean;
	dragHoverBar: BarPosition | null;

	moveIconToBar: (id: PanelId, bar: BarPosition, index?: number) => void;
	dropIconOnPanel: (id: PanelId, bar: BarPosition) => void;
	toggleActive: (bar: BarPosition, id: PanelId) => void;
	resetLayout: () => void;
	setBarSize: (bar: BarPosition, px: number) => void;
	setPanelDragActive: (active: boolean) => void;
	setDragHoverBar: (bar: BarPosition | null) => void;
}

function initialState() {
	return {
		leftBarIcons: getDefaultsForBar("left"),
		rightBarIcons: getDefaultsForBar("right"),
		bottomBarIcons: getDefaultsForBar("bottom"),
		activeLeft: null as PanelId | null,
		activeRight: null as PanelId | null,
		activeBottom: null as PanelId | null,
		leftWidth: SIDEBAR_DEFAULT_PX,
		rightWidth: SIDEBAR_DEFAULT_PX,
		bottomHeight: BOTTOM_DEFAULT_PX,
		isPanelDragActive: false,
		dragHoverBar: null as BarPosition | null,
	};
}

const BAR_KEY: Record<
	BarPosition,
	"leftBarIcons" | "rightBarIcons" | "bottomBarIcons"
> = {
	left: "leftBarIcons",
	right: "rightBarIcons",
	bottom: "bottomBarIcons",
};

const ACTIVE_KEY: Record<
	BarPosition,
	"activeLeft" | "activeRight" | "activeBottom"
> = {
	left: "activeLeft",
	right: "activeRight",
	bottom: "activeBottom",
};

export const createActivityBarSlice: StateCreator<
	AppState,
	[],
	[],
	ActivityBarSlice
> = (set) => ({
	...initialState(),

	moveIconToBar: (id, bar, index) =>
		set((state) => {
			const sourceBar: BarPosition | null = state.leftBarIcons.includes(id)
				? "left"
				: state.rightBarIcons.includes(id)
					? "right"
					: state.bottomBarIcons.includes(id)
						? "bottom"
						: null;
			if (!sourceBar) return {};
			if (PANELS[id].kind !== "tool") return {};

			const patch: Partial<ActivityBarSlice> = {};

			// Remove from source bar
			const sourceKey = BAR_KEY[sourceBar];
			const sourceList = state[sourceKey].filter((x) => x !== id);
			patch[sourceKey] = sourceList as PanelId[] & never;

			// Clear active in source bar if this icon was active
			if (state[ACTIVE_KEY[sourceBar]] === id) {
				patch[ACTIVE_KEY[sourceBar]] = null as never;
			}

			// Insert into target bar (same or different)
			const targetKey = BAR_KEY[bar];
			const currentTarget = sourceBar === bar ? sourceList : state[targetKey];
			const nextTarget = [...currentTarget];
			const insertAt =
				index === undefined || index < 0 || index > nextTarget.length
					? nextTarget.length
					: index;
			nextTarget.splice(insertAt, 0, id);
			patch[targetKey] = nextTarget as PanelId[] & never;

			return patch as ActivityBarSlice;
		}),

	dropIconOnPanel: (id, bar) =>
		set((state) => {
			// Move the icon, then force the target bar's active panel to the
			// dropped icon (so its content opens on drop).
			const sourceBar: BarPosition | null = state.leftBarIcons.includes(id)
				? "left"
				: state.rightBarIcons.includes(id)
					? "right"
					: state.bottomBarIcons.includes(id)
						? "bottom"
						: null;
			if (!sourceBar) return {};
			if (PANELS[id].kind !== "tool") return {};

			const patch: Partial<ActivityBarSlice> = {};

			const sourceKey = BAR_KEY[sourceBar];
			const sourceList = state[sourceKey].filter((x) => x !== id);
			patch[sourceKey] = sourceList as PanelId[] & never;

			if (state[ACTIVE_KEY[sourceBar]] === id) {
				patch[ACTIVE_KEY[sourceBar]] = null as never;
			}

			const targetKey = BAR_KEY[bar];
			const currentTarget = sourceBar === bar ? sourceList : state[targetKey];
			if (!currentTarget.includes(id)) {
				patch[targetKey] = [...currentTarget, id] as PanelId[] & never;
			}

			patch[ACTIVE_KEY[bar]] = id as never;

			return patch as ActivityBarSlice;
		}),

	toggleActive: (bar, id) =>
		set((state) => {
			const key = ACTIVE_KEY[bar];
			const current = state[key];
			return { [key]: current === id ? null : id } as Partial<ActivityBarSlice>;
		}),

	resetLayout: () => set(initialState()),

	setBarSize: (bar, px) =>
		set(() => {
			if (bar === "bottom") {
				const clamped = Math.min(BOTTOM_MAX_PX, Math.max(BOTTOM_MIN_PX, px));
				return { bottomHeight: clamped };
			}
			const clamped = Math.min(SIDEBAR_MAX_PX, Math.max(SIDEBAR_MIN_PX, px));
			return bar === "left" ? { leftWidth: clamped } : { rightWidth: clamped };
		}),

	setPanelDragActive: (active) =>
		set(() =>
			active
				? { isPanelDragActive: true }
				: { isPanelDragActive: false, dragHoverBar: null },
		),

	setDragHoverBar: (bar) => set(() => ({ dragHoverBar: bar })),
});
