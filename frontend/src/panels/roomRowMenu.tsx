import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import DeleteIcon from "@mui/icons-material/Delete";
import DownloadIcon from "@mui/icons-material/Download";
import DuplicateIcon from "@mui/icons-material/FileCopy";
import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import MoreVertIcon from "@mui/icons-material/MoreVert";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import {
	IconButton,
	ListItemIcon,
	ListItemText,
	Menu,
	MenuItem,
	Tooltip,
} from "@mui/material";
import { useState } from "react";
import DuplicateRoomDialog from "../components/DuplicateRoomDialog";
import {
	downloadFrames,
	type Room,
	setDefaultRoom,
	updateRoom,
} from "../myapi/client";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";

interface Props {
	room: Room;
}

export function RoomRowMenu({ room }: Props) {
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const rooms = useRoomsStore((s) => s.roomsArray);
	const [anchor, setAnchor] = useState<HTMLElement | null>(null);
	const [duplicateOpen, setDuplicateOpen] = useState(false);
	const open = Boolean(anchor);

	const onSetTemplate = async () => {
		setAnchor(null);
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

	const onToggleLock = async () => {
		setAnchor(null);
		try {
			await updateRoom(room.id, { locked: !room.locked });
			showSnackbar(room.locked ? "Room unlocked" : "Room locked", "success");
		} catch {
			showSnackbar("Failed to update lock", "error");
		}
	};

	const onDuplicate = () => {
		setAnchor(null);
		setDuplicateOpen(true);
	};

	const onCopyLink = async () => {
		setAnchor(null);
		await navigator.clipboard.writeText(
			`${window.location.origin}/rooms/${room.id}`,
		);
		showSnackbar("Link copied", "success");
	};

	const onDownload = () => {
		setAnchor(null);
		downloadFrames({ roomId: room.id });
		showSnackbar("Downloading all frames", "success");
	};

	return (
		<>
			<IconButton
				size="small"
				data-testid={`room-row-menu-${room.id}`}
				onClick={(e) => setAnchor(e.currentTarget)}
				onMouseDown={(e) => e.stopPropagation()}
				aria-label="Room actions"
			>
				<MoreVertIcon fontSize="small" />
			</IconButton>
			<Menu anchorEl={anchor} open={open} onClose={() => setAnchor(null)}>
				<MenuItem onClick={onSetTemplate}>
					<ListItemIcon>
						{room.is_default ? <StarIcon /> : <StarBorderIcon />}
					</ListItemIcon>
					<ListItemText>
						{room.is_default ? "Remove template" : "Set as template"}
					</ListItemText>
				</MenuItem>
				<MenuItem onClick={onDuplicate}>
					<ListItemIcon>
						<DuplicateIcon />
					</ListItemIcon>
					<ListItemText>Duplicate room</ListItemText>
				</MenuItem>
				<MenuItem onClick={onToggleLock}>
					<ListItemIcon>
						{room.locked ? <LockOpenIcon /> : <LockIcon />}
					</ListItemIcon>
					<ListItemText>{room.locked ? "Unlock" : "Lock"}</ListItemText>
				</MenuItem>
				<MenuItem onClick={onCopyLink}>
					<ListItemIcon>
						<ContentCopyIcon />
					</ListItemIcon>
					<ListItemText>Copy link</ListItemText>
				</MenuItem>
				<MenuItem onClick={onDownload}>
					<ListItemIcon>
						<DownloadIcon />
					</ListItemIcon>
					<ListItemText>Download frames</ListItemText>
				</MenuItem>
				<Tooltip title="Deleting rooms is not yet supported by the backend.">
					<span>
						<MenuItem disabled>
							<ListItemIcon>
								<DeleteIcon />
							</ListItemIcon>
							<ListItemText>Delete</ListItemText>
						</MenuItem>
					</span>
				</Tooltip>
			</Menu>
			<DuplicateRoomDialog
				open={duplicateOpen}
				sourceRoomId={room.id}
				sourceDescription={room.description ?? room.id}
				existingRoomIds={rooms.map((r) => r.id)}
				onClose={() => setDuplicateOpen(false)}
			/>
		</>
	);
}
