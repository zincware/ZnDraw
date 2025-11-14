import { useState, useEffect } from "react";
import { useParams, useNavigate } from "react-router-dom";
import IconButton from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";
import Menu from "@mui/material/Menu";
import MenuItem from "@mui/material/MenuItem";
import ListItemIcon from "@mui/material/ListItemIcon";
import ListItemText from "@mui/material/ListItemText";
import Dialog from "@mui/material/Dialog";
import DialogTitle from "@mui/material/DialogTitle";
import DialogContent from "@mui/material/DialogContent";
import DialogActions from "@mui/material/DialogActions";
import TextField from "@mui/material/TextField";
import Button from "@mui/material/Button";
import Chip from "@mui/material/Chip";
import Typography from "@mui/material/Typography";
import MoreVertIcon from "@mui/icons-material/MoreVert";
import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import ListIcon from "@mui/icons-material/List";
import FolderOpenIcon from "@mui/icons-material/FolderOpen";
import PowerSettingsNewIcon from "@mui/icons-material/PowerSettingsNew";
import CameraAltIcon from "@mui/icons-material/CameraAlt";
import DownloadIcon from "@mui/icons-material/Download";
import CircularProgress from "@mui/material/CircularProgress";
import {
  getRoom,
  updateRoom,
  duplicateRoom,
  setDefaultRoom,
  getDefaultRoom,
  listRooms,
  RoomDetail,
  getFileBrowserConfig,
  shutdownServer,
  downloadFrames,
  listFilesystems,
} from "../myapi/client";
import { takeAndUploadScreenshot } from "../utils/screenshot";
import { useExtensionData } from "../hooks/useSchemas";
import { useAppStore } from "../store";
import { useRoomsStore } from "../roomsStore";
import { socket } from "../socket";
import CloudIcon from "@mui/icons-material/Cloud";

interface DuplicateFormState {
  newRoomId: string;
  description: string;
  error: string | null;
}

/**
 * RoomManagementMenu provides room management actions in the AppBar:
 * - Lock/Unlock room
 * - Set as default
 * - Duplicate room
 * - Go to room list
 */
export default function RoomManagementMenu() {
  const { roomId } = useParams<{ roomId: string }>();
  const navigate = useNavigate();
  // Use individual selectors to prevent unnecessary re-renders
  const userName = useAppStore((state) => state.userName);
  const currentFrame = useAppStore((state) => state.currentFrame);
  const showSnackbar = useAppStore((state) => state.showSnackbar);
  const lockMetadata = useAppStore((state) => state.lockMetadata);

  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const [roomDetail, setRoomDetail] = useState<RoomDetail | null>(null);
  const [isDefault, setIsDefault] = useState(false);
  const [fileBrowserEnabled, setFileBrowserEnabled] = useState(false);
  const [remoteFilesystemsAvailable, setRemoteFilesystemsAvailable] = useState(false);
  const [duplicateDialog, setDuplicateDialog] = useState(false);
  const [duplicateForm, setDuplicateForm] = useState<DuplicateFormState>({
    newRoomId: "",
    description: "",
    error: null,
  });
  const [shutdownDialog, setShutdownDialog] = useState(false);
  const [screenshotLoading, setScreenshotLoading] = useState(false);

  // Fetch camera settings to check preserve_drawing_buffer
  const { data: cameraSettings } = useExtensionData(
    roomId || "",
    "settings",
    "camera",
  );

  // Subscribe to rooms from Zustand store (triggers re-render on changes)
  const rooms = useRoomsStore((state) => state.roomsArray);

  // Subscribe to current room (triggers re-render when this room updates)
  const currentRoomFromStore = useRoomsStore((state) =>
    roomId ? state.getRoom(roomId) : undefined
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

  // Check if file browser is enabled
  useEffect(() => {
    const checkFileBrowser = async () => {
      const config = await getFileBrowserConfig();
      setFileBrowserEnabled(config.enabled);
    };
    checkFileBrowser();
  }, []);

  // Check if remote filesystems are available
  useEffect(() => {
    const checkRemoteFilesystems = async () => {
      if (!roomId) return;
      try {
        const data = await listFilesystems(roomId);
        setRemoteFilesystemsAvailable(data.length > 0);
      } catch (error) {
        // If error, assume no filesystems available
        setRemoteFilesystemsAvailable(false);
      }
    };

    checkRemoteFilesystems();

    // Listen for filesystem updates via Socket.IO
    const handleFilesystemsUpdate = () => {
      checkRemoteFilesystems();
    };

    socket.on("filesystems:update", handleFilesystemsUpdate);

    // Cleanup listener on unmount
    return () => {
      socket.off("filesystems:update", handleFilesystemsUpdate);
    };
  }, [roomId]);

  const menuOpen = Boolean(anchorEl);

  const handleOpenMenu = async (event: React.MouseEvent<HTMLElement>) => {
    setAnchorEl(event.currentTarget);
    
    // Fetch room details when opening menu
    if (roomId) {
      try {
        const [detail, defaultRoom] = await Promise.all([
          getRoom(roomId),
          getDefaultRoom(),
        ]);
        setRoomDetail(detail);
        setIsDefault(defaultRoom.roomId === roomId);
      } catch (err) {
        console.error("Failed to fetch room details:", err);
      }
    }
  };

  const handleCloseMenu = () => {
    setAnchorEl(null);
  };

  const handleToggleLock = async () => {
    if (!roomId || !roomDetail) return;
    
    try {
      await updateRoom(roomId, { locked: !roomDetail.locked });
      setRoomDetail({ ...roomDetail, locked: !roomDetail.locked });
      showSnackbar(roomDetail.locked ? "Room unlocked" : "Room locked", "success");
    } catch (err) {
      showSnackbar("Failed to update lock status", "error");
    }
    handleCloseMenu();
  };

  const handleToggleDefault = async () => {
    if (!roomId) return;
    
    try {
      await setDefaultRoom(isDefault ? null : roomId);
      setIsDefault(!isDefault);
      showSnackbar(isDefault ? "Default room cleared" : "Set as default room", "success");
    } catch (err) {
      showSnackbar("Failed to update default room", "error");
    }
    handleCloseMenu();
  };

  const validateRoomId = (roomIdToCheck: string): string | null => {
    if (!roomIdToCheck) {
      return null; // Empty means auto-generate, which is valid
    }
    
    // Check if room ID already exists
    if (rooms.some((room) => room.id === roomIdToCheck)) {
      return "A room with this ID already exists";
    }
    
    // Basic validation for room ID format (alphanumeric, hyphens, underscores)
    if (!/^[a-zA-Z0-9_-]+$/.test(roomIdToCheck)) {
      return "Room ID can only contain letters, numbers, hyphens, and underscores";
    }
    
    return null;
  };

  const handleOpenDuplicateDialog = () => {
    setDuplicateForm({
      newRoomId: "",
      description: `Copy of ${roomDetail?.description || roomId || "room"}`,
      error: null,
    });
    setDuplicateDialog(true);
    handleCloseMenu();
  };

  const handleCloseDuplicateDialog = () => {
    setDuplicateDialog(false);
    setDuplicateForm({
      newRoomId: "",
      description: "",
      error: null,
    });
  };

  const handleDuplicateRoom = async () => {
    if (!roomId) return;
    
    // Validate room ID if provided
    const validationError = validateRoomId(duplicateForm.newRoomId);
    if (validationError) {
      setDuplicateForm({ ...duplicateForm, error: validationError });
      return;
    }
    
    try {
      const result = await duplicateRoom(roomId, {
        newRoomId: duplicateForm.newRoomId || undefined,
        description: duplicateForm.description,
      });
      
      handleCloseDuplicateDialog();

      // Navigate directly to the new room
      navigate(`/rooms/${result.roomId}`);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : "Failed to duplicate room";
      setDuplicateForm({
        ...duplicateForm,
        error: errorMessage,
      });
    }
  };

  const handleGoToRoomList = () => {
    navigate("/rooms");
    handleCloseMenu();
  };

  const handleGoToFileBrowser = () => {
    navigate("/file-browser");
    handleCloseMenu();
  };

  const handleGoToRemoteFileBrowser = () => {
    if (!roomId) return;
    navigate(`/rooms/${roomId}/remote-files`);
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

  const handleTakeScreenshot = async () => {
    if (!roomId) {
      showSnackbar("No room ID available", "error");
      return;
    }

    // Check if preserve_drawing_buffer is enabled
    if (!cameraSettings?.preserve_drawing_buffer) {
      showSnackbar("Enable 'preserve_drawing_buffer' in Camera settings first", "error");
      return;
    }

    setScreenshotLoading(true);
    handleCloseMenu();

    try {
      const canvas = document.querySelector("canvas");

      if (!canvas) {
        throw new Error("Canvas element not found");
      }

      const result = await takeAndUploadScreenshot(canvas, roomId);

      showSnackbar("Screenshot saved successfully", "success");
    } catch (error) {
      showSnackbar(`Screenshot failed: ${error instanceof Error ? error.message : "Unknown error"}`, "error");
    } finally {
      setScreenshotLoading(false);
    }
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

  return (
    <>
      {/* Metadata lock indicator (yellow) - vis.lock for uploads */}
      {lockMetadata?.locked && !roomDetail?.locked && (
        <Tooltip
          title={
            lockMetadata.userName
              ? `Locked by ${lockMetadata.userName}${lockMetadata.timestamp ? ` (${Math.floor((Date.now() - lockMetadata.timestamp * 1000) / 1000)}s ago)` : ''}`
              : "Room is locked"
          }
          arrow
        >
          <Chip
            icon={<LockIcon />}
            label={lockMetadata.msg || "Uploading"}
            color="warning"
            size="small"
            sx={{ mr: 1 }}
          />
        </Tooltip>
      )}
      
      {/* Room lock indicator (red) - permanent lock */}
      {roomDetail?.locked && (
        <Chip
          icon={<LockIcon />}
          label="Locked"
          color="error"
          size="small"
          sx={{ mr: 1 }}
        />
      )}
      
      {/* Default room indicator - always visible */}
      {isDefault && (
        <Chip
          icon={<StarIcon />}
          label="Default"
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
        <MenuItem onClick={handleToggleLock}>
          <ListItemIcon>
            {roomDetail?.locked ? <LockOpenIcon /> : <LockIcon />}
          </ListItemIcon>
          <ListItemText>
            {roomDetail?.locked ? "Unlock Room" : "Lock Room"}
          </ListItemText>
        </MenuItem>

        <MenuItem onClick={handleToggleDefault}>
          <ListItemIcon>
            {isDefault ? <StarBorderIcon /> : <StarIcon />}
          </ListItemIcon>
          <ListItemText>
            {isDefault ? "Remove as Default" : "Set as Default"}
          </ListItemText>
        </MenuItem>

        <MenuItem onClick={handleOpenDuplicateDialog}>
          <ListItemIcon>
            <ContentCopyIcon />
          </ListItemIcon>
          <ListItemText>Duplicate Room</ListItemText>
        </MenuItem>

        <MenuItem onClick={handleTakeScreenshot} disabled={screenshotLoading}>
          <ListItemIcon>
            {screenshotLoading ? <CircularProgress size={20} /> : <CameraAltIcon />}
          </ListItemIcon>
          <ListItemText>Take Screenshot</ListItemText>
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

        <MenuItem onClick={handleGoToRoomList}>
          <ListItemIcon>
            <ListIcon />
          </ListItemIcon>
          <ListItemText>Go to Room List</ListItemText>
        </MenuItem>

        {fileBrowserEnabled && (
          <MenuItem onClick={handleGoToFileBrowser}>
            <ListItemIcon>
              <FolderOpenIcon />
            </ListItemIcon>
            <ListItemText>Open File Browser</ListItemText>
          </MenuItem>
        )}

        {remoteFilesystemsAvailable && (
          <MenuItem onClick={handleGoToRemoteFileBrowser}>
            <ListItemIcon>
              <CloudIcon />
            </ListItemIcon>
            <ListItemText>Remote Filesystems</ListItemText>
          </MenuItem>
        )}

        <MenuItem onClick={handleOpenShutdownDialog}>
          <ListItemIcon>
            <PowerSettingsNewIcon />
          </ListItemIcon>
          <ListItemText>Close ZnDraw</ListItemText>
        </MenuItem>
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
            This will close the ZnDraw server for all connected users. Are you sure?
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
      <Dialog
        open={duplicateDialog}
        onClose={handleCloseDuplicateDialog}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Duplicate Room</DialogTitle>
        <DialogContent>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 3, mt: 1 }}>
            Duplicating: {roomDetail?.description || roomId}
          </Typography>
          
          <TextField
            margin="dense"
            label="New Room ID (optional)"
            type="text"
            fullWidth
            variant="outlined"
            value={duplicateForm.newRoomId}
            onChange={(e) => {
              const newId = e.target.value;
              setDuplicateForm({
                ...duplicateForm,
                newRoomId: newId,
                error: validateRoomId(newId),
              });
            }}
            helperText={
              duplicateForm.error || 
              "Leave empty to auto-generate a unique ID"
            }
            error={!!duplicateForm.error}
            sx={{ mb: 2 }}
          />
          
          <TextField
            autoFocus
            margin="dense"
            label="Description for new room"
            type="text"
            fullWidth
            variant="outlined"
            value={duplicateForm.description}
            onChange={(e) =>
              setDuplicateForm({ ...duplicateForm, description: e.target.value })
            }
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={handleCloseDuplicateDialog}>Cancel</Button>
          <Button 
            onClick={handleDuplicateRoom} 
            variant="contained"
            disabled={!!duplicateForm.error}
          >
            Duplicate
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
}
