import { useWindowManagerStore } from "../../stores/windowManagerStore";
import type { HandlerContext } from "./types";

/** Socket event payload for the `figure_invalidate` event. */
export interface FigureInvalidateEvent {
	key: string;
	operation?: "set" | "delete";
}

const MAX_AUTO_OPEN_WINDOWS = 5;

/**
 * Creates figure handler functions.
 *
 * Uses `useWindowManagerStore.getState()` for imperative access to window
 * state (not a hook call). Window open/close operations interact with
 * windowManagerStore, not appStore.
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

			// Step 2: Close any open windows displaying this figure
			const windowsToClose = Object.entries(
				useWindowManagerStore.getState().openWindows,
			)
				.filter(([_, window]) => window.figureKey === data.key)
				.map(([windowId]) => windowId);

			windowsToClose.forEach((windowId) => {
				useWindowManagerStore.getState().closeWindow(windowId);
			});

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

			// Step 3: Auto-open window if not already displaying this figure
			const openWindowsState = useWindowManagerStore.getState().openWindows;
			const isWindowDisplayingFigure = Object.values(openWindowsState).some(
				(window) => window.figureKey === data.key,
			);

			const openWindowCount = Object.keys(openWindowsState).length;

			if (!isWindowDisplayingFigure) {
				if (openWindowCount < MAX_AUTO_OPEN_WINDOWS) {
					ctx.openWindow(data.key);
				}
			}
		}
	}

	return { onFiguresInvalidate };
}
