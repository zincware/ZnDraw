import type { DockviewApi } from "dockview-react";
import { useCallback } from "react";
import { useAppStore } from "../store";

interface LeaveRoomOptions {
	api: DockviewApi | null;
	showConfirm?: (message: string) => Promise<boolean>;
}

export function useLeaveRoom({ api, showConfirm }: LeaveRoomOptions) {
	const showSnackbar = useAppStore((s) => s.showSnackbar);

	return useCallback(
		async (opts?: { skipConfirm?: boolean }) => {
			if (!api) return;
			const plotPanels = api.panels.filter((p) => p.id.startsWith("plot-"));
			const hasPlots = plotPanels.length > 0;

			if (hasPlots && !opts?.skipConfirm) {
				const msg = `Leave room? ${plotPanels.length} plot(s) will close.`;
				const ok = showConfirm ? await showConfirm(msg) : window.confirm(msg);
				if (!ok) return;
			}

			for (const panel of plotPanels) {
				panel.api.close();
			}
			const viewer = api.getPanel("viewer");
			if (viewer) viewer.api.close();

			window.history.pushState({}, "", "/");
			showSnackbar("Left room", "info");
		},
		[api, showConfirm, showSnackbar],
	);
}
