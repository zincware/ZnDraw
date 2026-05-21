import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import DownloadIcon from "@mui/icons-material/Download";
import FolderOpenIcon from "@mui/icons-material/FolderOpen";
import LockIcon from "@mui/icons-material/Lock";
import MoreVertIcon from "@mui/icons-material/MoreVert";
import PowerSettingsNewIcon from "@mui/icons-material/PowerSettingsNew";
import ShareIcon from "@mui/icons-material/Share";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import Button from "@mui/material/Button";
import Chip from "@mui/material/Chip";
import Dialog from "@mui/material/Dialog";
import DialogActions from "@mui/material/DialogActions";
import DialogContent from "@mui/material/DialogContent";
import DialogTitle from "@mui/material/DialogTitle";
import IconButton from "@mui/material/IconButton";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import Menu from "@mui/material/Menu";
import MenuItem from "@mui/material/MenuItem";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import { useEffect, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { composeRoomAddress } from "../utils/roomAddress";
import {
	type RoomDetail,
	type Visibility,
	downloadFrames,
	getRoom,
	listProviders,
	setDefaultRoom,
	shutdownServer,
	updateRoom,
} from "../myapi/client";
import { useRoomsStore } from "../roomsStore";
import { socket } from "../socket";
import { useAppStore } from "../store";
import DuplicateRoomDialog from "./DuplicateRoomDialog";
import ShareDialog from "./ShareDialog";
import VisibilitySelector from "./VisibilitySelector";

/**
 * RoomManagementMenu provides room management actions in the AppBar:
 * - Lock/Unlock room
 * - Set as template
 * - Duplicate room
 * - Go to room list
 */
export default function RoomManagementMenu() {
	const { ownerId, roomName } = useParams<{ ownerId: string; roomName: string }>();
	const roomId = ownerId && roomName ? composeRoomAddress(ownerId, roomName) : undefined;
	const navigate = useNavigate();
	// Use individual selectors to prevent unnecessary re-renders
	const userName = useAppStore((state) => state.user?.email ?? null);
	const userId = useAppStore((state) => state.user?.id ?? null);
	const isAdmin = useAppStore((state) => state.user?.is_superuser ?? false);
	const currentFrame = useAppStore((state) => state.currentFrame);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const userLock = useAppStore((state) => state.userLock);
	const userLockMessage = useAppStore((state) => state.userLockMessage);

	const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
	const [roomDetail, setRoomDetail] = useState<RoomDetail | null>(null);
	const [filesystemAvailable, setFilesystemAvailable] = useState(false);
	const [duplicateOpen, setDuplicateOpen] = useState(false);
	const [shutdownDialog, setShutdownDialog] = useState(false);
	const [shareOpen, setShareOpen] = useState(false);
	const [visibility, setVisibility] = useState<Visibility>("public");

	// Subscribe to rooms from Zustand store (triggers re-render on changes)
	const rooms = useRoomsStore((state) => state.roomsArray);

	// Subscribe to current room (triggers re-render when this room updates)
	const currentRoomFromStore = useRoomsStore((state) =>
		roomId ? state.getRoom(roomId) : undefined,
	);

	// Fetch rooms on mount
	useEffect(() => {
		useRoomsStore.getState().fetchRooms();
	}, []);

	// Update local roomDetail state when room changes in store
	useEffect(() => {
		if (currentRoomFromStore) {
			setRoomDetail(currentRoomFromStore);
		}
	}, [currentRoomFromStore]);

	// Sync visibility from store
	useEffect(() => {
		if (currentRoomFromStore?.visibility) {
			setVisibility(currentRoomFromStore.visibility);
		}
	}, [currentRoomFromStore?.visibility]);

	// Derive isDefault from room data (either from store or local state)
	const isDefault =
		roomDetail?.is_default ?? currentRoomFromStore?.is_default ?? false;

	// Check if filesystem providers are available
	useEffect(() => {
		const checkFilesystems = async () => {
			if (!roomId) return;
			try {
				const providers = await listProviders(roomId, "filesystem");
				setFilesystemAvailable(providers.length > 0);
			} catch {
				setFilesystemAvailable(false);
			}
		};

		checkFilesystems();

		const handleProvidersUpdate = () => {
			checkFilesystems();
		};

		socket.on("providers_invalidate", handleProvidersUpdate);

		return () => {
			socket.off("providers_invalidate", handleProvidersUpdate);
		};
	}, [roomId]);

	const menuOpen = Boolean(anchorEl);

	const isEditLocked = userLock !== null;

	const handleOpenMenu = async (event: React.MouseEvent<HTMLElement>) => {
		setAnchorEl(event.currentTarget);

		// Fetch room details when opening menu
		if (roomId) {
			try {
				const detail = await getRoom(roomId);
				setRoomDetail(detail);
			} catch (err) {
				console.error("Failed to fetch room details:", err);
			}
		}
	};

	const handleCloseMenu = () => {
		setAnchorEl(null);
	};

	const handleLockIconClick = () => {
		if (isEditLocked) {
			showSnackbar(
				`Locked by ${userLock || "another user"}: ${userLockMessage || "in use"}`,
				"info",
			);
		}
	};

	const handleVisibility = async (v: Visibility) => {
		if (!roomId) return;
		setVisibility(v);
		try {
			await updateRoom(roomId, { visibility: v });
			useRoomsStore.getState().updateRoom(roomId, { visibility: v });
			showSnackbar(`Visibility updated to ${v}`, "success");
		} catch {
			showSnackbar("Failed to update visibility", "error");
			if (currentRoomFromStore?.visibility) {
				setVisibility(currentRoomFromStore.visibility);
			}
		}
	};

	const canManage =
		isAdmin ||
		(userId !== null && ownerId === userId);

	const handleToggleDefault = async () => {
		if (!roomId) return;

		try {
			await setDefaultRoom(isDefault ? null : roomId);

			// Update the room in the store immediately for instant UI feedback
			useRoomsStore.getState().updateRoom(roomId, { is_default: !isDefault });

			showSnackbar(
				isDefault ? "Template cleared" : "Set as template",
				"success",
			);
		} catch {
			showSnackbar("Failed to update template", "error");
			// Revert the optimistic update on error
			useRoomsStore.getState().updateRoom(roomId, { is_default: isDefault });
		}
		handleCloseMenu();
	};

	const handleOpenDuplicateDialog = () => {
		setDuplicateOpen(true);
		handleCloseMenu();
	};

	const handleGoToFilesystem = () => {
		if (!ownerId || !roomName) return;
		navigate(`/rooms/${ownerId}/${roomName}/files`);
		handleCloseMenu();
	};

	const handleOpenShutdownDialog = () => {
		setShutdownDialog(true);
		handleCloseMenu();
	};

	const handleCloseShutdownDialog = () => {
		setShutdownDialog(false);
	};

	const handleShutdownServer = async () => {
		setShutdownDialog(false);
		// Call shutdown - server will shut down before it can respond, so don't show error
		shutdownServer().catch(() => {
			// Expected to fail since server shuts down immediately
		});
	};

	const handleDownloadCurrentFrame = () => {
		if (!roomId) return;

		// Send the actual current frame index from frontend state
		downloadFrames({
			roomId,
			indices: [currentFrame],
		});

		handleCloseMenu();

		showSnackbar("Downloading current frame as ExtendedXYZ", "success");
	};

	const handleDownloadAllFrames = () => {
		if (!roomId) return;

		// No parameters = download all frames
		downloadFrames({
			roomId,
		});

		handleCloseMenu();

		showSnackbar("Downloading all frames as ExtendedXYZ", "success");
	};

	if (!roomId) {
		return null;
	}

	// Build tooltip text for edit-lock icon
	const getLockTooltip = () => {
		if (isEditLocked) {
			const user = userLock || "Someone";
			const action = userLockMessage || "using this room";
			return `${user}: ${action}`;
		}
		return "";
	};

	return (
		<>
			{/* Edit-lock indicator — shown when another user holds the edit lock */}
			{isEditLocked && (
				<Tooltip title={getLockTooltip()} arrow>
					<IconButton
						size="small"
						onClick={handleLockIconClick}
						sx={{ color: "warning.main", mr: 0.5 }}
					>
						<LockIcon />
					</IconButton>
				</Tooltip>
			)}

			{/* Template room indicator - always visible */}
			{isDefault && (
				<Chip
					icon={<StarIcon />}
					label="Template"
					color="primary"
					size="small"
					sx={{ mr: 1 }}
				/>
			)}

			{/* Room menu button */}
			<Tooltip title="Room menu">
				<IconButton
					color="inherit"
					aria-label="room menu"
					onClick={handleOpenMenu}
				>
					<MoreVertIcon />
				</IconButton>
			</Tooltip>

			{/* Menu with room management options */}
			<Menu
				anchorEl={anchorEl}
				open={menuOpen}
				onClose={handleCloseMenu}
				anchorOrigin={{
					vertical: "bottom",
					horizontal: "right",
				}}
				transformOrigin={{
					vertical: "top",
					horizontal: "right",
				}}
			>
				<MenuItem disableRipple sx={{ display: "block", px: 2, py: 1 }}>
					<Typography variant="caption" sx={{ mb: 0.5, display: "block" }}>
						Visibility
					</Typography>
					<VisibilitySelector
						value={visibility}
						onChange={handleVisibility}
						disabled={!canManage}
					/>
				</MenuItem>

				{isAdmin && (
					<MenuItem onClick={handleToggleDefault}>
						<ListItemIcon>
							{isDefault ? <StarBorderIcon /> : <StarIcon />}
						</ListItemIcon>
						<ListItemText>
							{isDefault ? "Remove as Template" : "Set as Template"}
						</ListItemText>
					</MenuItem>
				)}

				{canManage && (
					<MenuItem
						onClick={() => {
							setShareOpen(true);
							handleCloseMenu();
						}}
					>
						<ListItemIcon>
							<ShareIcon fontSize="small" />
						</ListItemIcon>
						<ListItemText>Share…</ListItemText>
					</MenuItem>
				)}

				<MenuItem onClick={handleOpenDuplicateDialog}>
					<ListItemIcon>
						<ContentCopyIcon />
					</ListItemIcon>
					<ListItemText>Duplicate Room</ListItemText>
				</MenuItem>

				<MenuItem onClick={handleDownloadCurrentFrame}>
					<ListItemIcon>
						<DownloadIcon />
					</ListItemIcon>
					<ListItemText>Current Frame (ExtXYZ)</ListItemText>
				</MenuItem>

				<MenuItem onClick={handleDownloadAllFrames}>
					<ListItemIcon>
						<DownloadIcon />
					</ListItemIcon>
					<ListItemText>All Frames (ExtXYZ)</ListItemText>
				</MenuItem>

				{filesystemAvailable && (
					<MenuItem onClick={handleGoToFilesystem}>
						<ListItemIcon>
							<FolderOpenIcon />
						</ListItemIcon>
						<ListItemText>Filesystem</ListItemText>
					</MenuItem>
				)}

				{isAdmin && (
					<MenuItem onClick={handleOpenShutdownDialog}>
						<ListItemIcon>
							<PowerSettingsNewIcon />
						</ListItemIcon>
						<ListItemText>Close ZnDraw</ListItemText>
					</MenuItem>
				)}
			</Menu>

			{/* Shutdown Confirmation Dialog */}
			<Dialog
				open={shutdownDialog}
				onClose={handleCloseShutdownDialog}
				maxWidth="xs"
				fullWidth
			>
				<DialogTitle>Shutdown Server?</DialogTitle>
				<DialogContent>
					<Typography variant="body2" color="text.secondary" sx={{ mt: 1 }}>
						This will close the ZnDraw server for all connected users. Are you
						sure?
					</Typography>
				</DialogContent>
				<DialogActions>
					<Button onClick={handleCloseShutdownDialog}>Cancel</Button>
					<Button
						onClick={handleShutdownServer}
						variant="contained"
						color="error"
					>
						Shutdown
					</Button>
				</DialogActions>
			</Dialog>

			{/* Duplicate Dialog */}
			<DuplicateRoomDialog
				open={duplicateOpen}
				sourceRoomId={roomId || ""}
				sourceDescription={roomDetail?.description || roomId || "room"}
				existingRoomIds={rooms.map((r) => r.room_id)}
				onClose={() => setDuplicateOpen(false)}
			/>

			{/* Share Dialog */}
			{roomId && (
				<ShareDialog
					roomId={roomId}
					open={shareOpen}
					onClose={() => setShareOpen(false)}
				/>
			)}
		</>
	);
}
