import SearchIcon from "@mui/icons-material/Search";
import {
	Box,
	InputAdornment,
	List,
	ListItemButton,
	ListItemText,
	TextField,
	Typography,
} from "@mui/material";
import { useEffect, useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";
import { getDockviewApi } from "./DockviewLayout";
import { RoomsHeaderActions } from "./roomsHeaderActions";
import { RoomRowMenu } from "./roomRowMenu";

/**
 * VS-Code-style compact rooms list with search, header actions, and
 * a per-row actions menu. Reactive via the existing roomsStore.
 */
export function RoomsPanel() {
	const rooms = useRoomsStore((s) => s.roomsArray);
	const loading = useRoomsStore((s) => s.loading);
	const fetchRooms = useRoomsStore((s) => s.fetchRooms);
	const currentRoomId = useAppStore((s) => s.roomId);
	const leaveRoom = useLeaveRoom({ api: getDockviewApi() });
	const navigate = useNavigate();
	const [query, setQuery] = useState("");

	useEffect(() => {
		fetchRooms();
	}, [fetchRooms]);

	const filtered = useMemo(() => {
		const q = query.trim().toLowerCase();
		if (!q) return rooms;
		return rooms.filter(
			(r) =>
				r.id.toLowerCase().includes(q) ||
				(r.description && r.description.toLowerCase().includes(q)),
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
			sx={{ display: "flex", flexDirection: "column", height: "100%" }}
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
			<List dense sx={{ flexGrow: 1, overflow: "auto", pt: 0 }}>
				{filtered.map((r) => (
					<ListItemButton
						key={r.id}
						data-testid={`rooms-row-${r.id}`}
						selected={r.id === currentRoomId}
						onClick={() => switchToRoom(r.id)}
						sx={{
							alignItems: "center",
							pl: 1.5,
							pr: 0.5,
							borderLeft: 2,
							borderColor:
								r.id === currentRoomId ? "primary.main" : "transparent",
						}}
					>
						<ListItemText
							primary={r.description ?? r.id}
							secondary={`${r.id.slice(0, 8)} · ${r.frame_count} frame${
								r.frame_count === 1 ? "" : "s"
							}`}
							primaryTypographyProps={{
								variant: "body2",
								fontWeight: 500,
								noWrap: true,
							}}
							secondaryTypographyProps={{
								variant: "caption",
								color: "text.secondary",
							}}
						/>
						<RoomRowMenu room={r} />
					</ListItemButton>
				))}
			</List>
		</Box>
	);
}
