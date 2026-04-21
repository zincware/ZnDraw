import type { DockviewApi } from "dockview-react";

export function plotPanelId(figureKey: string): string {
	return `plot-${figureKey}`;
}

const PLOT_ID_PREFIX = "plot-";
const VIEWER_PANEL_ID = "viewer";

/**
 * Resolves which group a new plot tab should open in.
 *
 * Rules (first match wins):
 * 1. If any existing `plot-*` panel is open, pick the group containing
 *    the plot whose id sorts first lexicographically. Stable across
 *    sessions and independent of insertion order.
 * 2. Else, place a new group to the right of the viewer's group.
 * 3. Else (no viewer), let dockview decide.
 */
function resolvePlotPosition(
	api: DockviewApi,
): { referenceGroup: string; direction: "within" | "right" } | undefined {
	const plotPanels = api.panels
		.filter((p) => p.id.startsWith(PLOT_ID_PREFIX))
		.sort((a, b) => a.id.localeCompare(b.id));

	if (plotPanels.length > 0) {
		const group = plotPanels[0].group;
		if (group) {
			return { referenceGroup: group.id, direction: "within" };
		}
	}

	const viewer = api.getPanel(VIEWER_PANEL_ID);
	if (viewer?.group) {
		return { referenceGroup: viewer.group.id, direction: "right" };
	}

	return undefined;
}

export function openPlotTab(
	api: DockviewApi,
	figureKey: string,
	opts?: {
		referenceGroupId?: string;
		direction?: "right" | "left" | "above" | "below" | "within";
	},
): void {
	const id = plotPanelId(figureKey);
	const existing = api.getPanel(id);
	if (existing) {
		existing.api.setActive();
		return;
	}

	// Explicit caller override (used by the onDidDrop path in DockviewLayout
	// where the user dropped onto a specific target).
	if (opts?.referenceGroupId) {
		api.addPanel({
			id,
			component: "plotView",
			title: figureKey,
			params: { figureKey },
			position: {
				referenceGroup: opts.referenceGroupId,
				direction: opts.direction ?? "within",
			},
		});
		return;
	}

	const position = resolvePlotPosition(api);
	api.addPanel({
		id,
		component: "plotView",
		title: figureKey,
		params: { figureKey },
		position,
	});
}

export function closePlotTab(api: DockviewApi, figureKey: string): void {
	const panel = api.getPanel(plotPanelId(figureKey));
	panel?.api.close();
}
