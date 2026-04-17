import Badge from "@mui/material/Badge";
import { Box, IconButton, Tooltip, keyframes } from "@mui/material";
import { useCallback, useMemo, useRef } from "react";
import { useShallow } from "zustand/react/shallow";
import { useHasFilesystemProviders } from "../hooks/useHasFilesystemProviders";
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

const shimmer = keyframes`
	0%, 100% { background-color: rgba(25, 118, 210, 0.12); }
	50% { background-color: rgba(25, 118, 210, 0.28); }
`;

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
	const isDragActive = useAppStore((s) => s.isPanelDragActive);
	const hoverBar = useAppStore((s) => s.dragHoverBar);
	const setHoverBar = useAppStore((s) => s.setDragHoverBar);
	const setPanelDragActive = useAppStore((s) => s.setPanelDragActive);

	const isHovered = hoverBar === position;
	const dragDepth = useRef(0);

	const hasFilesystemProviders = useHasFilesystemProviders();
	const visibleIcons = useMemo(
		() =>
			hasFilesystemProviders
				? icons
				: icons.filter((id) => id !== "filesystem"),
		[icons, hasFilesystemProviders],
	);

	const onDragStart = useCallback((e: React.DragEvent, id: PanelId) => {
		e.dataTransfer.setData(DRAG_MIME, id);
		e.dataTransfer.effectAllowed = "move";
	}, []);

	const onDragOver = useCallback((e: React.DragEvent) => {
		if (e.dataTransfer.types.includes(DRAG_MIME)) {
			e.preventDefault();
			e.dataTransfer.dropEffect = "move";
		}
	}, []);

	const onDragEnter = useCallback(
		(e: React.DragEvent) => {
			if (!e.dataTransfer.types.includes(DRAG_MIME)) return;
			dragDepth.current++;
			if (dragDepth.current === 1) setHoverBar(position);
		},
		[position, setHoverBar],
	);

	const onDragLeave = useCallback(() => {
		dragDepth.current = Math.max(0, dragDepth.current - 1);
		if (
			dragDepth.current === 0 &&
			useAppStore.getState().dragHoverBar === position
		) {
			setTimeout(() => {
				if (
					dragDepth.current === 0 &&
					useAppStore.getState().dragHoverBar === position
				) {
					setHoverBar(null);
				}
			}, 0);
		}
	}, [position, setHoverBar]);

	const onDropOnBar = useCallback(
		(e: React.DragEvent) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			dragDepth.current = 0;
			setPanelDragActive(false);
			if (!id) return;
			e.preventDefault();
			// Drop on the bar → move only. Don't touch the panel's active state.
			moveIconToBar(id as PanelId, position);
		},
		[moveIconToBar, position, setPanelDragActive],
	);

	const onDropOnIcon = useCallback(
		(e: React.DragEvent, overIdx: number) => {
			const id = e.dataTransfer.getData(DRAG_MIME) as PanelId | "";
			dragDepth.current = 0;
			setPanelDragActive(false);
			if (!id) return;
			e.preventDefault();
			e.stopPropagation();
			moveIconToBar(id as PanelId, position, overIdx);
		},
		[moveIconToBar, position, setPanelDragActive],
	);

	// Hide the bar entirely when it has no icons and no drag is in progress.
	if (visibleIcons.length === 0 && !isDragActive) return null;

	const dragBgSx = isDragActive
		? isHovered
			? { bgcolor: "rgba(25, 118, 210, 0.32)" }
			: { animation: `${shimmer} 1.2s ease-in-out infinite` }
		: {};

	return (
		<Box
			data-testid={`activity-bar-${position}`}
			data-drop-hover={isHovered}
			onDragEnter={onDragEnter}
			onDragOver={onDragOver}
			onDragLeave={onDragLeave}
			onDrop={onDropOnBar}
			sx={{
				display: "flex",
				borderColor: "divider",
				bgcolor: "background.paper",
				alignItems: "center",
				transition: "background-color 120ms ease",
				...BAR_SX[position],
				...dragBgSx,
			}}
		>
			{visibleIcons.map((id, idx) => {
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
