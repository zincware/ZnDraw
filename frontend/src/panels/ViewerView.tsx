import { Box } from "@mui/material";
import type { IDockviewPanelProps } from "dockview-react";
import { useEffect, useRef, useState } from "react";
import MyScene from "../components/Canvas";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import { useDockviewApi } from "../stores/dockviewApiStore";

export function ViewerView(props: IDockviewPanelProps) {
	const { api } = props;
	const [remountKey, setRemountKey] = useState(0);
	const previousLocation = useRef(api.location.type);
	const dockApi = useDockviewApi((s) => s.api);
	const leaveRoom = useLeaveRoom({ api: dockApi });
	const leaveRoomRef = useRef(leaveRoom);
	useEffect(() => {
		leaveRoomRef.current = leaveRoom;
	}, [leaveRoom]);

	useEffect(() => {
		const locationDisposable = api.onDidLocationChange((e) => {
			if (e.location.type !== previousLocation.current) {
				previousLocation.current = e.location.type;
				setRemountKey((k) => k + 1);
			}
		});

		return () => {
			locationDisposable.dispose();
		};
	}, [api]);

	// When the viewer panel is removed (user closed the tab), run the
	// leave-room cascade. DockviewPanelApi does not expose onDidDispose,
	// so we listen on the top-level DockviewApi's onDidRemovePanel and
	// filter by panel id. The cascade closes any remaining plot-* panels
	// and resets the URL to "/".
	useEffect(() => {
		if (!dockApi) return;
		const removeDisposable = dockApi.onDidRemovePanel((panel) => {
			if (panel.id === api.id) {
				void leaveRoomRef.current({ skipConfirm: true });
			}
		});
		return () => {
			removeDisposable.dispose();
		};
	}, [dockApi, api.id]);

	return (
		<Box
			key={remountKey}
			data-testid="viewer-view"
			sx={{
				position: "relative",
				width: "100%",
				height: "100%",
				overflow: "hidden",
			}}
		>
			<MyScene />
		</Box>
	);
}
