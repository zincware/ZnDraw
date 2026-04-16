import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import {
	Box,
	IconButton,
	List,
	ListItem,
	ListItemButton,
	ListItemText,
	Typography,
} from "@mui/material";
import { useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { useLeaveRoom } from "../hooks/useLeaveRoom";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";
import { getDockviewApi } from "./DockviewLayout";

/**
 * RoomsPanel lists all rooms available on the server and lets the user switch
 * between them.
 *
 * Minimum viable feature set for v1:
 * - List rooms (sorted by ``roomsStore``).
 * - Highlight the current room.
 * - Click a row to switch — performs the leave-room cascade (close plots,
 *   close viewer) then navigates to the selected room.
 * - Copy link button per row.
 *
 * Notes
 * -----
 * Per-room star / leave / delete actions are intentionally not included here;
 * they remain in the AppBar ``RoomManagementMenu`` until a later task.
 */
export function RoomsPanel() {
	const rooms = useRoomsStore((s) => s.roomsArray);
	const loading = useRoomsStore((s) => s.loading);
	const fetchRooms = useRoomsStore((s) => s.fetchRooms);
	const currentRoomId = useAppStore((s) => s.roomId);
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const leaveRoom = useLeaveRoom({ api: getDockviewApi() });
	const navigate = useNavigate();

	useEffect(() => {
		fetchRooms();
	}, [fetchRooms]);

	const switchToRoom = async (roomId: string) => {
		if (roomId === currentRoomId) return;
		// User explicitly chose to switch — skip the confirm prompt.
		await leaveRoom({ skipConfirm: true });
		navigate(`/rooms/${roomId}`);
	};

	const copyLink = (roomId: string) => {
		navigator.clipboard.writeText(`${window.location.origin}/rooms/${roomId}`);
		showSnackbar("Link copied", "success");
	};

	if (loading) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography>Loading rooms…</Typography>
			</Box>
		);
	}
	if (rooms.length === 0) {
		return (
			<Box sx={{ p: 2 }}>
				<Typography color="text.secondary">
					No rooms available — create one via the Filesystem panel.
				</Typography>
			</Box>
		);
	}
	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Typography variant="overline" sx={{ px: 2, pt: 1 }}>
				Rooms
			</Typography>
			<List dense sx={{ flexGrow: 1, overflow: "auto" }}>
				{rooms.map((r) => (
					<ListItem
						key={r.id}
						secondaryAction={
							<IconButton size="small" onClick={() => copyLink(r.id)}>
								<ContentCopyIcon fontSize="small" />
							</IconButton>
						}
					>
						<ListItemButton
							selected={r.id === currentRoomId}
							onClick={() => switchToRoom(r.id)}
						>
							<ListItemText
								primary={r.description ?? r.id}
								secondary={r.id}
							/>
						</ListItemButton>
					</ListItem>
				))}
			</List>
		</Box>
	);
}
