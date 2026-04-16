import Badge from "@mui/material/Badge";
import { Box, IconButton, Tooltip } from "@mui/material";
import { useCallback } from "react";
import { useShallow } from "zustand/react/shallow";
import { useAppStore } from "../store";
import { type BarPosition, type PanelId, PANELS } from "./registry";

const DRAG_MIME = "application/x-zndraw-panel-id";

const BAR_SX = {
	left: { flexDirection: "column", width: 48, borderRight: 1 },
	right: { flexDirection: "column", width: 48, borderLeft: 1 },
	bottom: { flexDirection: "row", height: 40, borderTop: 1, width: "100%" },
} as const;

const ACTIVE_INDICATOR: Record<
	BarPosition,
	(color: string) => Record<string, string>
> = {
	left: (c) => ({ borderLeft: `2px solid ${c}` }),
	right: (c) => ({ borderRight: `2px solid ${c}` }),
	bottom: (c) => ({ borderTop: `2px solid ${c}` }),
};

interface ActivityBarProps {
	position: BarPosition;
}

export function ActivityBar({ position }: ActivityBarProps) {
	const icons = useAppStore(
		useShallow((s) =>
			position === "left"
				? s.leftBarIcons
				: position === "right"
					? s.rightBarIcons
					: s.bottomBarIcons,
		),
	);
	const active = useAppStore((s) =>
		position === "left"
			? s.activeLeft
			: position === "right"
				? s.activeRight
				: s.activeBottom,
	);
	const toggleActive = useAppStore((s) => s.toggleActive);
	const moveIconToBar = useAppStore((s) => s.moveIconToBar);
	const chatUnread = useAppStore((s) => s.chatUnreadCount);

	const onDragStart = useCallback(
		(e: React.DragEvent, id: PanelId) => {
			e.dataTransfer.setData(DRAG_MIME, id);
			e.dataTransfer.effectAllowed = "move";
		},
		[],
	);

	const onDragOver = useCallback((e: React.DragEvent) => {
		if (e.dataTransfer.types.includes(DRAG_MIME)) {
			e.preventDefault();
			e.dataTransfer.dropEffect = "move";
		}
	}, []);

	const onDropOnBar = useCallback(
		(e: React.DragEvent) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			if (!id) return;
			e.preventDefault();
			moveIconToBar(id as PanelId, position);
		},
		[moveIconToBar, position],
	);

	const onDropOnIcon = useCallback(
		(e: React.DragEvent, overIdx: number) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			if (!id) return;
			e.preventDefault();
			e.stopPropagation();
			moveIconToBar(id as PanelId, position, overIdx);
		},
		[moveIconToBar, position],
	);

	return (
		<Box
			data-testid={`activity-bar-${position}`}
			onDragOver={onDragOver}
			onDrop={onDropOnBar}
			sx={{
				display: "flex",
				borderColor: "divider",
				bgcolor: "background.paper",
				alignItems: "center",
				...BAR_SX[position],
			}}
		>
			{icons.map((id, idx) => {
				const def = PANELS[id];
				if (def.kind !== "tool") return null;
				const Icon = def.icon;
				const isActive = active === id;
				const iconElement = <Icon />;
				const iconWithBadge =
					id === "chat" && !isActive && chatUnread > 0 ? (
						<Badge badgeContent={chatUnread} color="error" max={99}>
							{iconElement}
						</Badge>
					) : (
						iconElement
					);
				return (
					<Tooltip
						key={id}
						title={def.label}
						placement={
							position === "left"
								? "right"
								: position === "right"
									? "left"
									: "top"
						}
					>
						<IconButton
							data-testid={`activity-icon-${id}`}
							draggable
							onDragStart={(e) => onDragStart(e, id)}
							onDragOver={onDragOver}
							onDrop={(e) => onDropOnIcon(e, idx)}
							onClick={() => toggleActive(position, id)}
							color={isActive ? "primary" : "default"}
							sx={{
								borderRadius: 0,
								m: position === "bottom" ? "0 2px" : "2px 0",
								...(isActive
									? ACTIVE_INDICATOR[position](
											"var(--mui-palette-primary-main)",
										)
									: {}),
							}}
						>
							{iconWithBadge}
						</IconButton>
					</Tooltip>
				);
			})}
		</Box>
	);
}
