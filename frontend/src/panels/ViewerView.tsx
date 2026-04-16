import { Box } from "@mui/material";
import type { IDockviewPanelProps } from "dockview-react";
import { useEffect, useRef, useState } from "react";
import MyScene from "../components/Canvas";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import { getDockviewApi } from "./DockviewLayout";

export function ViewerView(props: IDockviewPanelProps) {
	const { api } = props;
	const [remountKey, setRemountKey] = useState(0);
	const previousLocation = useRef(api.location.type);
	const leaveRoom = useLeaveRoom({ api: getDockviewApi() });
	const leaveRoomRef = useRef(leaveRoom);
	leaveRoomRef.current = leaveRoom;

	useEffect(() => {
		const locationDisposable = api.onDidLocationChange((e) => {
			if (e.location.type !== previousLocation.current) {
				previousLocation.current = e.location.type;
				setRemountKey((k) => k + 1);
			}
		});

		// When the viewer panel is removed (user closed the tab), run the
		// leave-room cascade. DockviewPanelApi does not expose onDidDispose,
		// so we listen on the top-level DockviewApi's onDidRemovePanel and
		// filter by panel id. The cascade closes any remaining plot-* panels
		// and resets the URL to "/".
		const dockApi = getDockviewApi();
		const removeDisposable = dockApi?.onDidRemovePanel((panel) => {
			if (panel.id === api.id) {
				void leaveRoomRef.current({ skipConfirm: true });
			}
		});

		return () => {
			locationDisposable.dispose();
			removeDisposable?.dispose();
		};
	}, [api]);

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
