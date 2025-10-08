import { useState, useEffect } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { useQuery } from "@tanstack/react-query";
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
import Snackbar from "@mui/material/Snackbar";
import Alert from "@mui/material/Alert";
import Chip from "@mui/material/Chip";
import Typography from "@mui/material/Typography";
import SettingsIcon from "@mui/icons-material/Settings";
import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import ListIcon from "@mui/icons-material/List";
import FolderOpenIcon from "@mui/icons-material/FolderOpen";
import {
  getRoom,
  updateRoom,
  duplicateRoom,
  setDefaultRoom,
  getDefaultRoom,
  listRooms,
  RoomDetail,
  getFileBrowserConfig,
} from "../myapi/client";

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
  
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
  const [roomDetail, setRoomDetail] = useState<RoomDetail | null>(null);
  const [isDefault, setIsDefault] = useState(false);
  const [fileBrowserEnabled, setFileBrowserEnabled] = useState(false);
  const [duplicateDialog, setDuplicateDialog] = useState(false);
  const [duplicateForm, setDuplicateForm] = useState<DuplicateFormState>({
    newRoomId: "",
    description: "",
    error: null,
  });
  const [snackbar, setSnackbar] = useState<{
    open: boolean;
    message: string;
    severity: "success" | "error";
  }>({ open: false, message: "", severity: "success" });

  // Fetch all rooms to check for duplicates and get current room metadata
  const { data: rooms = [] } = useQuery({
    queryKey: ["rooms"],
    queryFn: listRooms,
    refetchInterval: 5000,
  });
  
  // Update roomDetail when rooms data changes
  useEffect(() => {
    if (roomId && rooms.length > 0) {
      const currentRoom = rooms.find((room) => room.id === roomId);
      if (currentRoom) {
        setRoomDetail(currentRoom);
      }
    }
  }, [roomId, rooms]);

  // Check if file browser is enabled
  useEffect(() => {
    const checkFileBrowser = async () => {
      const config = await getFileBrowserConfig();
      setFileBrowserEnabled(config.enabled);
    };
    checkFileBrowser();
  }, []);

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
      setSnackbar({
        open: true,
        message: roomDetail.locked ? "Room unlocked" : "Room locked",
        severity: "success",
      });
    } catch (err) {
      setSnackbar({
        open: true,
        message: "Failed to update lock status",
        severity: "error",
      });
    }
    handleCloseMenu();
  };

  const handleToggleDefault = async () => {
    if (!roomId) return;
    
    try {
      await setDefaultRoom(isDefault ? null : roomId);
      setIsDefault(!isDefault);
      setSnackbar({
        open: true,
        message: isDefault ? "Default room cleared" : "Set as default room",
        severity: "success",
      });
    } catch (err) {
      setSnackbar({
        open: true,
        message: "Failed to update default room",
        severity: "error",
      });
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
      const userId = crypto.randomUUID();
      navigate(`/rooms/${result.roomId}/${userId}`);
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

  if (!roomId) {
    return null;
  }

  return (
    <>
      {/* Metadata lock indicator (yellow) - vis.lock for uploads */}
      {roomDetail?.metadataLocked && !roomDetail?.locked && (
        <Chip
          icon={<LockIcon />}
          label="Uploading"
          color="warning"
          size="small"
          sx={{ mr: 1 }}
        />
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

      {/* Room settings menu button */}
      <Tooltip title="Room settings">
        <IconButton
          color="inherit"
          aria-label="room settings"
          onClick={handleOpenMenu}
        >
          <SettingsIcon />
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
      </Menu>

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

      {/* Snackbar for notifications */}
      <Snackbar
        open={snackbar.open}
        autoHideDuration={4000}
        onClose={() => setSnackbar({ ...snackbar, open: false })}
        anchorOrigin={{ vertical: "bottom", horizontal: "center" }}
      >
        <Alert
          onClose={() => setSnackbar({ ...snackbar, open: false })}
          severity={snackbar.severity}
          sx={{ width: "100%" }}
        >
          {snackbar.message}
        </Alert>
      </Snackbar>
    </>
  );
}
