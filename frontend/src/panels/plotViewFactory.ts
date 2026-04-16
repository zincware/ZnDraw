import type { DockviewApi } from "dockview-react";

export function plotPanelId(figureKey: string): string {
	return `plot-${figureKey}`;
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
	const referenceGroupId = opts?.referenceGroupId;
	api.addPanel({
		id,
		component: "plotView",
		title: figureKey,
		params: { figureKey },
		position: referenceGroupId
			? {
					referenceGroup: referenceGroupId,
					direction: opts?.direction ?? "within",
				}
			: undefined,
	});
}

export function closePlotTab(api: DockviewApi, figureKey: string): void {
	const panel = api.getPanel(plotPanelId(figureKey));
	panel?.api.close();
}
