import type { DockviewApi } from "dockview-react";
import { useCallback } from "react";
import { useNavigate } from "react-router-dom";
import { useAppStore } from "../store";

type ApiResolver = DockviewApi | null | (() => DockviewApi | null);

interface LeaveRoomOptions {
	api: ApiResolver;
	showConfirm?: (message: string) => Promise<boolean>;
}

export function useLeaveRoom({ api, showConfirm }: LeaveRoomOptions) {
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const navigate = useNavigate();

	return useCallback(
		async (opts?: { skipConfirm?: boolean }) => {
			const resolved = typeof api === "function" ? api() : api;
			if (!resolved) return;
			const plotPanels = resolved.panels.filter((p) =>
				p.id.startsWith("plot-"),
			);
			const hasPlots = plotPanels.length > 0;

			if (hasPlots && !opts?.skipConfirm) {
				const msg = `Leave room? ${plotPanels.length} plot(s) will close.`;
				const ok = showConfirm ? await showConfirm(msg) : window.confirm(msg);
				if (!ok) return;
			}

			for (const panel of plotPanels) {
				panel.api.close();
			}
			const viewer = resolved.getPanel("viewer");
			if (viewer) viewer.api.close();

			navigate("/");
			showSnackbar("Left room", "info");
		},
		[api, showConfirm, showSnackbar, navigate],
	);
}
