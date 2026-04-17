import Badge from "@mui/material/Badge";
import { Box, IconButton, Tooltip, keyframes } from "@mui/material";
import { useCallback, useEffect, useState } from "react";
import { useShallow } from "zustand/react/shallow";
import { useAppStore } from "../store";
import { type BarPosition, type PanelId, PANELS } from "./registry";

const DRAG_MIME = "application/x-zndraw-panel-id";

const BAR_SX = {
	left: { flexDirection: "column", width: 48, borderRight: 1 },
	right: { flexDirection: "column", width: 48, borderLeft: 1 },
	bottom: { flexDirection: "row", height: 40, borderTop: 1, width: "100%" },
} as const;

const SLIVER_IDLE_SX = {
	left: { width: 4 },
	right: { width: 4 },
	bottom: { height: 4, width: "100%" },
} as const;

// Hot zone (dragging) widens the empty bar to 56 px so it's an easy drop target.
const SLIVER_HOT_SX = {
	left: { width: 56 },
	right: { width: 56 },
	bottom: { height: 56, width: "100%" },
} as const;

const ACTIVE_INDICATOR: Record<
	BarPosition,
	(color: string) => Record<string, string>
> = {
	left: (c) => ({ borderLeft: `2px solid ${c}` }),
	right: (c) => ({ borderRight: `2px solid ${c}` }),
	bottom: (c) => ({ borderTop: `2px solid ${c}` }),
};

const pulse = keyframes`
	0%, 100% { opacity: 0.5; }
	50% { opacity: 1; }
`;

type SliverState = "full" | "sliver" | "hot" | "over-zone";

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

	const [isDragActive, setIsDragActive] = useState(false);
	const [isOverZone, setIsOverZone] = useState(false);

	useEffect(() => {
		const onGlobalDragStart = (e: DragEvent) => {
			if (e.dataTransfer?.types.includes(DRAG_MIME)) {
				setIsDragActive(true);
			}
		};
		const onGlobalDragEnd = () => {
			setIsDragActive(false);
			setIsOverZone(false);
		};
		window.addEventListener("dragstart", onGlobalDragStart, { passive: true });
		window.addEventListener("dragend", onGlobalDragEnd, { passive: true });
		return () => {
			window.removeEventListener("dragstart", onGlobalDragStart);
			window.removeEventListener("dragend", onGlobalDragEnd);
		};
	}, []);

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

	const state: SliverState =
		icons.length === 0
			? isDragActive
				? isOverZone
					? "over-zone"
					: "hot"
				: "sliver"
			: "full";

	if (state !== "full") {
		const sizeSx = state === "sliver" ? SLIVER_IDLE_SX[position] : SLIVER_HOT_SX[position];
		const isHot = state === "hot";
		const isOver = state === "over-zone";
		const hintByPosition: Record<BarPosition, string> = {
			left: "Drop to dock left",
			right: "Drop to dock right",
			bottom: "Drop to dock bottom",
		};
		return (
			<Box
				data-testid={`activity-bar-${position}`}
				data-sliver-state={state}
				onDragOver={(e) => {
					onDragOver(e);
					if (isDragActive) setIsOverZone(true);
				}}
				onDragLeave={() => setIsOverZone(false)}
				onDrop={(e) => {
					setIsOverZone(false);
					onDropOnBar(e);
				}}
				sx={{
					display: "flex",
					flexShrink: 0,
					alignItems: "center",
					justifyContent: "center",
					position: "relative",
					...sizeSx,
					bgcolor: isOver
						? "rgba(25, 118, 210, 0.35)"
						: isHot
							? "rgba(25, 118, 210, 0.15)"
							: "transparent",
					borderColor: "primary.main",
					borderStyle: isOver ? "solid" : isHot ? "dashed" : "none",
					borderWidth:
						isHot || isOver
							? position === "left"
								? "0 2px 0 0"
								: position === "right"
									? "0 0 0 2px"
									: "2px 0 0 0"
							: 0,
					animation: isHot ? `${pulse} 1s ease-in-out infinite` : "none",
				}}
			>
				{isOver && (
					<Box
						component="span"
						sx={{
							fontSize: 11,
							fontWeight: 500,
							color: "primary.main",
							whiteSpace: "nowrap",
							transform:
								position === "left"
									? "rotate(-90deg)"
									: position === "right"
										? "rotate(90deg)"
										: "none",
							pointerEvents: "none",
						}}
					>
						{hintByPosition[position]}
					</Box>
				)}
			</Box>
		);
	}

	return (
		<Box
			data-testid={`activity-bar-${position}`}
			data-sliver-state="full"
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
