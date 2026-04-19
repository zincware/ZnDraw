import FullscreenIcon from "@mui/icons-material/Fullscreen";
import FullscreenExitIcon from "@mui/icons-material/FullscreenExit";
import OpenInNewIcon from "@mui/icons-material/OpenInNew";
import { Box, IconButton, Tooltip } from "@mui/material";
import type { IDockviewHeaderActionsProps } from "dockview-react";
import { useEffect, useState } from "react";

export function GroupActions(props: IDockviewHeaderActionsProps) {
	const { api, containerApi, group } = props;
	const [isMax, setIsMax] = useState<boolean>(() => api.isMaximized());

	useEffect(() => {
		const update = () => setIsMax(api.isMaximized());
		const disposable = containerApi.onDidMaximizedGroupChange(update);
		update();
		return () => disposable.dispose();
	}, [api, containerApi]);

	const onPopout = () => {
		containerApi.addPopoutGroup(group);
	};

	const onToggleMaximize = () => {
		if (api.isMaximized()) {
			api.exitMaximized();
		} else {
			api.maximize();
		}
	};

	return (
		<Box sx={{ display: "flex", alignItems: "center", px: 0.5 }}>
			<Tooltip title="Pop out">
				<IconButton
					size="small"
					onClick={onPopout}
					data-testid="group-popout"
					aria-label="Pop out"
				>
					<OpenInNewIcon fontSize="small" />
				</IconButton>
			</Tooltip>
			<Tooltip title={isMax ? "Exit full-screen" : "Full-screen"}>
				<IconButton
					size="small"
					onClick={onToggleMaximize}
					data-testid="group-maximize"
					aria-label={isMax ? "Exit full-screen" : "Full-screen"}
				>
					{isMax ? (
						<FullscreenExitIcon fontSize="small" />
					) : (
						<FullscreenIcon fontSize="small" />
					)}
				</IconButton>
			</Tooltip>
		</Box>
	);
}
