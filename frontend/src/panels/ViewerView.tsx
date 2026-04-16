import { Box } from "@mui/material";
import type { IDockviewPanelProps } from "dockview-react";
import { useEffect, useRef, useState } from "react";
import MyScene from "../components/Canvas";

export function ViewerView(props: IDockviewPanelProps) {
	const { api } = props;
	const [remountKey, setRemountKey] = useState(0);
	const previousLocation = useRef(api.location.type);

	useEffect(() => {
		const disposable = api.onDidLocationChange((e) => {
			if (e.location.type !== previousLocation.current) {
				previousLocation.current = e.location.type;
				setRemountKey((k) => k + 1);
			}
		});
		return () => disposable.dispose();
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
