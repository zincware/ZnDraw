import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import SearchIcon from "@mui/icons-material/Search";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import {
	Box,
	IconButton,
	InputAdornment,
	List,
	ListItemButton,
	ListItemText,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import {
	createRoom,
	type Room,
	setDefaultRoom,
	updateRoom,
	uploadTrajectory,
} from "../myapi/client";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";
import { useDockviewApi } from "../stores/dockviewApiStore";
import { extractDetail } from "../utils/errors";
import { RoomRowMenu } from "./roomRowMenu";
import { RoomsHeaderActions } from "./roomsHeaderActions";

/**
 * Compact rooms list with search, header actions, per-row template/lock
 * toggles, and a per-row overflow menu. Reactive via the existing roomsStore.
 */
export function RoomsPanel() {
	const rooms = useRoomsStore((s) => s.roomsArray);
	const loading = useRoomsStore((s) => s.loading);
	const fetchRooms = useRoomsStore((s) => s.fetchRooms);
	const currentRoomId = useAppStore((s) => s.roomId);
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const dockApi = useDockviewApi((s) => s.api);
	const leaveRoom = useLeaveRoom({ api: dockApi });
	const navigate = useNavigate();
	const [query, setQuery] = useState("");

	const handleFiles = useCallback(
		async (files: File[]) => {
			const newRoomId = crypto.randomUUID();
			try {
				await createRoom({ room_id: newRoomId });
				// Cascade-close the current room before navigating so plot tabs
				// and viewer state from the prior room don't leak into the new one.
				await leaveRoom({ skipConfirm: true });
				for (const file of files) {
					await uploadTrajectory(newRoomId, file);
				}
				showSnackbar(`Room created with ${files.length} file(s)`, "success");
				navigate(`/rooms/${newRoomId}`);
			} catch (error) {
				showSnackbar(extractDetail(error, "Upload failed"), "error");
			}
		},
		[navigate, showSnackbar, leaveRoom],
	);

	const {
		isDragging,
		handleDragOver,
		handleDragEnter,
		handleDragLeave,
		handleDrop,
	} = useDragAndDrop(handleFiles);

	useEffect(() => {
		fetchRooms();
	}, [fetchRooms]);

	const filtered = useMemo(() => {
		const q = query.trim().toLowerCase();
		if (!q) return rooms;
		return rooms.filter(
			(r) =>
				r.id.toLowerCase().includes(q) ||
				r.description?.toLowerCase().includes(q),
		);
	}, [rooms, query]);

	const switchToRoom = async (id: string) => {
		if (id === currentRoomId) return;
		await leaveRoom({ skipConfirm: true });
		navigate(`/rooms/${id}`);
	};

	return (
		<Box
			data-testid="rooms-panel"
			onDragOver={handleDragOver}
			onDragEnter={handleDragEnter}
			onDragLeave={handleDragLeave}
			onDrop={handleDrop}
			sx={{
				display: "flex",
				flexDirection: "column",
				height: "100%",
				position: "relative",
			}}
		>
			<Box
				sx={{
					display: "flex",
					alignItems: "center",
					justifyContent: "space-between",
					px: 1.5,
					pt: 1,
					pb: 0.5,
				}}
			>
				<Typography
					variant="overline"
					sx={{ color: "text.secondary", letterSpacing: "0.08em" }}
				>
					Rooms
				</Typography>
				<RoomsHeaderActions />
			</Box>
			<Box sx={{ px: 1.5, pb: 1 }}>
				<TextField
					size="small"
					fullWidth
					value={query}
					onChange={(e) => setQuery(e.target.value)}
					placeholder="Search rooms…"
					data-testid="rooms-search"
					slotProps={{
						input: {
							startAdornment: (
								<InputAdornment position="start">
									<SearchIcon fontSize="small" />
								</InputAdornment>
							),
						},
					}}
				/>
			</Box>
			{loading && (
				<Box sx={{ px: 2 }}>
					<Typography variant="body2" color="text.secondary">
						Loading rooms…
					</Typography>
				</Box>
			)}
			{!loading && filtered.length === 0 && (
				<Box sx={{ px: 2, pt: 1 }}>
					<Typography variant="body2" color="text.secondary">
						{query
							? "No rooms match the search."
							: "No rooms available — create one above."}
					</Typography>
				</Box>
			)}
			<List dense sx={{ flexGrow: 1, minHeight: 0, overflow: "auto", pt: 0 }}>
				{filtered.map((r) => (
					<RoomsListRow
						key={r.id}
						room={r}
						selected={r.id === currentRoomId}
						onSelect={() => switchToRoom(r.id)}
					/>
				))}
			</List>
			{isDragging && (
				<Box
					data-testid="rooms-drop-overlay"
					sx={{
						position: "absolute",
						inset: 0,
						bgcolor: "rgba(25, 118, 210, 0.2)",
						border: 2,
						borderStyle: "dashed",
						borderColor: "primary.main",
						display: "flex",
						alignItems: "center",
						justifyContent: "center",
						pointerEvents: "none",
						zIndex: 2,
					}}
				>
					<Typography variant="body2" color="primary.main" fontWeight={500}>
						Drop to create a new room
					</Typography>
				</Box>
			)}
		</Box>
	);
}

interface RoomsListRowProps {
	room: Room;
	selected: boolean;
	onSelect: () => void;
}

function RoomsListRow({ room, selected, onSelect }: RoomsListRowProps) {
	const showSnackbar = useAppStore((s) => s.showSnackbar);

	const primary = room.description?.trim() || room.id;
	const secondary = `${room.frame_count} frame${
		room.frame_count === 1 ? "" : "s"
	}`;

	const onToggleTemplate = async (e: React.MouseEvent) => {
		e.stopPropagation();
		try {
			await setDefaultRoom(room.is_default ? null : room.id);
			useRoomsStore
				.getState()
				.updateRoom(room.id, { is_default: !room.is_default });
			showSnackbar(
				room.is_default ? "Template cleared" : "Set as template",
				"success",
			);
		} catch {
			showSnackbar("Failed to update template", "error");
		}
	};

	const onToggleLock = async (e: React.MouseEvent) => {
		e.stopPropagation();
		try {
			await updateRoom(room.id, { locked: !room.locked });
			showSnackbar(room.locked ? "Room unlocked" : "Room locked", "success");
		} catch {
			showSnackbar("Failed to update lock", "error");
		}
	};

	return (
		<ListItemButton
			data-testid={`rooms-row-${room.id}`}
			selected={selected}
			onClick={onSelect}
			sx={{
				gap: 0.5,
				alignItems: "center",
				pl: 0.5,
				pr: 0.5,
				borderLeft: 2,
				borderColor: selected ? "primary.main" : "transparent",
			}}
		>
			<Tooltip title={room.is_default ? "Template room" : "Set as template"}>
				<IconButton
					size="small"
					data-testid={`rooms-row-template-${room.id}`}
					onClick={onToggleTemplate}
				>
					{room.is_default ? (
						<StarIcon fontSize="small" color="primary" />
					) : (
						<StarBorderIcon fontSize="small" />
					)}
				</IconButton>
			</Tooltip>
			<ListItemText
				primary={primary}
				secondary={secondary}
				slotProps={{
					primary: {
						variant: "body2",
						fontWeight: 500,
						noWrap: true,
					},
					secondary: {
						variant: "caption",
						color: "text.secondary",
					},
				}}
			/>
			<Tooltip
				title={
					room.locked ? "Locked (click to unlock)" : "Unlocked (click to lock)"
				}
			>
				<IconButton
					size="small"
					data-testid={`rooms-row-lock-${room.id}`}
					onClick={onToggleLock}
				>
					{room.locked ? (
						<LockIcon fontSize="small" color="error" />
					) : (
						<LockOpenIcon fontSize="small" color="success" />
					)}
				</IconButton>
			</Tooltip>
			<RoomRowMenu room={room} />
		</ListItemButton>
	);
}
