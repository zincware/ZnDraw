import { useDockviewApi } from "../../stores/dockviewApiStore";
import { closePlotTab, openPlotTab } from "../../panels/plotViewFactory";
import type { HandlerContext } from "./types";

/** Socket event payload for the `figure_invalidate` event. */
export interface FigureInvalidateEvent {
	key: string;
	operation?: "set" | "delete";
}

/**
 * Creates figure handler functions.
 *
 * Drives the dockview API via `useDockviewApi.getState().api` to open/close plot tabs
 * when figures are created, updated, or deleted server-side.
 */
export function createFigureHandlers(ctx: HandlerContext) {
	function onFiguresInvalidate(data: FigureInvalidateEvent) {
		if (!data.key) return;

		const operation = data.operation || "set";

		if (operation === "delete") {
			// Step 1: Remove the figure data from the cache
			ctx.queryClient.removeQueries({
				queryKey: ["figures", ctx.roomId, "detail", data.key],
			});

			// Step 2: Close the dockview panel displaying this figure, if open
			const api = useDockviewApi.getState().api;
			if (api) closePlotTab(api, data.key);

			// Step 3: Invalidate the figures list to update the UI
			ctx.queryClient.invalidateQueries({
				queryKey: ["figures", ctx.roomId, "list"],
			});
		} else if (operation === "set") {
			// Step 1: Invalidate the data so it's fresh when needed
			ctx.queryClient.invalidateQueries({
				queryKey: ["figures", ctx.roomId, "detail", data.key],
			});

			// Step 2: Invalidate the list in case this is a new figure
			ctx.queryClient.invalidateQueries({
				queryKey: ["figures", ctx.roomId, "list"],
			});

			// Step 3: Auto-open as a tab — openPlotTab already no-ops if open.
			const api = useDockviewApi.getState().api;
			if (api) openPlotTab(api, data.key);
		}
	}

	return { onFiguresInvalidate };
}
