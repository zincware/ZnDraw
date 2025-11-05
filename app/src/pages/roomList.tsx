import React, { useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import Box from "@mui/material/Box";
import Container from "@mui/material/Container";
import Typography from "@mui/material/Typography";
import Button from "@mui/material/Button";
import IconButton from "@mui/material/IconButton";
import Tooltip from "@mui/material/Tooltip";
import CircularProgress from "@mui/material/CircularProgress";
import Alert from "@mui/material/Alert";
import Dialog from "@mui/material/Dialog";
import DialogTitle from "@mui/material/DialogTitle";
import DialogContent from "@mui/material/DialogContent";
import DialogActions from "@mui/material/DialogActions";
import TextField from "@mui/material/TextField";
import {
  DataGrid,
  GridColDef,
  GridRenderCellParams,
  GridRowParams,
} from "@mui/x-data-grid";
import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import VisibilityIcon from "@mui/icons-material/Visibility";
import VisibilityOffIcon from "@mui/icons-material/VisibilityOff";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import FolderOpenIcon from "@mui/icons-material/FolderOpen";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import {
  listRooms,
  updateRoom,
  duplicateRoom,
  setDefaultRoom,
  getFileBrowserConfig,
  Room,
} from "../myapi/client";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import DropOverlay from "../components/DropOverlay";
import { useAppStore } from "../store";
import { useRoomsStore } from "../roomsStore";
import { useSocketManager } from "../hooks/useSocketManager";

interface DuplicateDialogState {
  open: boolean;
  roomId: string;
  roomDescription: string;
}

interface DuplicateFormState {
  newRoomId: string;
  description: string;
  error: string | null;
}

export default function RoomListPage() {
  const queryClient = useQueryClient();
  // Use individual selectors to prevent unnecessary re-renders
  const showSnackbar = useAppStore((state) => state.showSnackbar);

  // Connect to socket and join overview:public room
  useSocketManager({ isOverview: true });

  // Subscribe to rooms store (triggers re-render on changes)
  const allRooms = useRoomsStore((state) => state.roomsArray);
  const loading = useRoomsStore((state) => state.loading);
  const roomsError = useRoomsStore((state) => state.error);

  const [searchQuery, setSearchQuery] = useState<string>("");
  const [duplicateDialog, setDuplicateDialog] = useState<DuplicateDialogState>({
    open: false,
    roomId: "",
    roomDescription: "",
  });
  const [duplicateForm, setDuplicateForm] = useState<DuplicateFormState>({
    newRoomId: "",
    description: "",
    error: null,
  });
  const [fileBrowserEnabled, setFileBrowserEnabled] = useState(false);
  const navigate = useNavigate();

  // Drag and drop support
  const { isDragging, handleDragOver, handleDragEnter, handleDragLeave, handleDrop } = useDragAndDrop();

  // File upload ref for button click
  const fileInputRef = React.useRef<HTMLInputElement>(null);

  const handleFileUploadClick = () => {
    fileInputRef.current?.click();
  };

  const handleFileInputChange = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const files = event.target.files;
    if (!files || files.length === 0) return;

    // Use the uploadFile API from drag/drop hook
    const file = files[0];
    const uploadEvent = {
      preventDefault: () => {},
      stopPropagation: () => {},
      dataTransfer: { files: [file] }
    } as any;

    await handleDrop(uploadEvent);

    // Reset input
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  // Check if file browser is enabled
  useEffect(() => {
    const checkFileBrowser = async () => {
      const config = await getFileBrowserConfig();
      setFileBrowserEnabled(config.enabled);
    };
    checkFileBrowser();
  }, []);

  // Fetch rooms on mount (socket updates will keep them synced)
  useEffect(() => {
    useRoomsStore.getState().fetchRooms();
  }, []);

  // Filter rooms based on search query
  const rooms = searchQuery
    ? allRooms.filter(
        (room) =>
          room.id.toLowerCase().includes(searchQuery.toLowerCase()) ||
          (room.description &&
            room.description.toLowerCase().includes(searchQuery.toLowerCase()))
      )
    : allRooms;

  const error = roomsError;
  const isError = !!roomsError;

  const handleUpdateDescription = async (roomId: string, description: string) => {
    try {
      await updateRoom(roomId, { description: description || null });
      // Socket event will update Zustand store automatically
      showSnackbar("Description updated", "success");
    } catch (err) {
      showSnackbar("Failed to update description", "error");
    }
  };

  const handleToggleLock = async (roomId: string, currentLocked: boolean) => {
    try {
      await updateRoom(roomId, { locked: !currentLocked });
      // Socket event will update Zustand store automatically
      showSnackbar(currentLocked ? "Room unlocked" : "Room locked", "success");
    } catch (err) {
      showSnackbar("Failed to update lock status", "error");
    }
  };

  const handleToggleHidden = async (roomId: string, currentHidden: boolean) => {
    try {
      await updateRoom(roomId, { hidden: !currentHidden });
      // Socket event will update Zustand store automatically
      showSnackbar(currentHidden ? "Room shown" : "Room hidden", "success");
    } catch (err) {
      showSnackbar("Failed to update visibility", "error");
    }
  };

  const handleToggleDefault = async (roomId: string, isCurrentlyDefault: boolean) => {
    try {
      await setDefaultRoom(isCurrentlyDefault ? null : roomId);
      // Socket event will update Zustand store automatically
      showSnackbar(isCurrentlyDefault ? "Default room cleared" : "Default room set", "success");
    } catch (err) {
      showSnackbar("Failed to update default room", "error");
    }
  };

  const handleOpenDuplicateDialog = (roomId: string, description: string) => {
    setDuplicateDialog({
      open: true,
      roomId,
      roomDescription: description,
    });
    setDuplicateForm({
      newRoomId: "",
      description: `Copy of ${description || roomId}`,
      error: null,
    });
  };

  const handleCloseDuplicateDialog = () => {
    setDuplicateDialog({ open: false, roomId: "", roomDescription: "" });
    setDuplicateForm({
      newRoomId: "",
      description: "",
      error: null,
    });
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

  const handleDuplicateRoom = async () => {
    // Validate room ID if provided
    const validationError = validateRoomId(duplicateForm.newRoomId);
    if (validationError) {
      setDuplicateForm({ ...duplicateForm, error: validationError });
      return;
    }

    try {
      const result = await duplicateRoom(duplicateDialog.roomId, {
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

  const handleOpenRoom = (roomId: string) => {
    navigate(`/rooms/${roomId}`);
  };

  const columns: GridColDef[] = [
    {
      field: "isDefault",
      headerName: "",
      width: 50,
      sortable: true,
      renderCell: (params: GridRenderCellParams) => (
        <Tooltip title={params.value ? "Default room" : "Set as default"}>
          <IconButton
            size="small"
            onClick={() => handleToggleDefault(params.row.id, params.value)}
          >
            {params.value ? (
              <StarIcon color="primary" />
            ) : (
              <StarBorderIcon />
            )}
          </IconButton>
        </Tooltip>
      ),
    },
    {
      field: "id",
      headerName: "Room ID",
      width: 250,
      sortable: true,
      renderCell: (params: GridRenderCellParams) => (
        <Box
          sx={{
            cursor: "pointer",
            color: "primary.main",
            textDecoration: "underline",
            "&:hover": {
              color: "primary.dark",
            },
          }}
          onClick={() => handleOpenRoom(params.value)}
        >
          {params.value}
        </Box>
      ),
    },
    {
      field: "description",
      headerName: "Description",
      width: 300,
      editable: true,
      sortable: true,
    },
    {
      field: "frameCount",
      headerName: "Frames",
      width: 100,
      type: "number",
      sortable: true,
    },
    {
      field: "locked",
      headerName: "Lock",
      width: 80,
      renderCell: (params: GridRenderCellParams) => {
        const isLocked = params.value;
        const isMetadataLocked = params.row.metadataLocked;
        
        // Priority: permanent lock > metadata lock > unlocked
        if (isLocked) {
          // Permanently locked - red lock
          return (
            <Tooltip title="Locked (immutable)">
              <IconButton
                size="small"
                onClick={() => handleToggleLock(params.row.id, isLocked)}
              >
                <LockIcon color="error" />
              </IconButton>
            </Tooltip>
          );
        } else if (isMetadataLocked) {
          // Uploading - yellow lock (non-clickable)
          return (
            <Tooltip title="Uploading...">
              <LockIcon sx={{ color: "warning.main", fontSize: 24 }} />
            </Tooltip>
          );
        } else {
          // Unlocked - green unlock icon
          return (
            <Tooltip title="Unlocked (editable)">
              <IconButton
                size="small"
                onClick={() => handleToggleLock(params.row.id, isLocked)}
              >
                <LockOpenIcon color="success" />
              </IconButton>
            </Tooltip>
          );
        }
      },
    },
    {
      field: "hidden",
      headerName: "Visible",
      width: 80,
      renderCell: (params: GridRenderCellParams) => (
        <Tooltip title={params.value ? "Hidden" : "Visible"}>
          <IconButton
            size="small"
            onClick={() => handleToggleHidden(params.row.id, params.value)}
          >
            {params.value ? (
              <VisibilityOffIcon color="disabled" />
            ) : (
              <VisibilityIcon color="action" />
            )}
          </IconButton>
        </Tooltip>
      ),
    },
    {
      field: "actions",
      headerName: "Actions",
      width: 80,
      sortable: false,
      renderCell: (params: GridRenderCellParams) => (
        <Tooltip title="Duplicate room">
          <IconButton
            size="small"
            onClick={() =>
              handleOpenDuplicateDialog(
                params.row.id,
                params.row.description || "",
              )
            }
          >
            <ContentCopyIcon />
          </IconButton>
        </Tooltip>
      ),
    },
  ];

  if (loading) {
    return (
      <Container maxWidth="lg">
        <Box
          sx={{
            display: "flex",
            justifyContent: "center",
            alignItems: "center",
            minHeight: "100vh",
          }}
        >
          <CircularProgress />
        </Box>
      </Container>
    );
  }

  if (isError) {
    return (
      <Container maxWidth="lg">
        <Box sx={{ mt: 4 }}>
          <Alert severity="error">
            {error instanceof Error ? error.message : "Unknown error"}
          </Alert>
        </Box>
      </Container>
    );
  }

  return (
    <Container
      maxWidth="lg"
      onDragOver={handleDragOver}
      onDragEnter={handleDragEnter}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
    >
      <DropOverlay isDragging={isDragging} />

      {/* Hidden file input for button upload */}
      <input
        type="file"
        ref={fileInputRef}
        style={{ display: 'none' }}
        onChange={handleFileInputChange}
        accept=".xyz,.extxyz,.pdb,.cif,.h5,.h5md,.hdf5,.gro,.mol,.sdf,.db,.json,.traj,.nc,.car,.xsf,.cube,.vasp,.poscar,.contcar,.xdatcar,.outcar,.xml,.pwi,.pwo,.out,.castep,.cell,.geom,.md,.gjf,.com,.log,.arc,.dmol"
      />

      <Box sx={{ mt: 4, mb: 4 }}>
        <Typography variant="h3" component="h1" gutterBottom>
          Room Management
        </Typography>
        <Typography variant="subtitle1" color="text.secondary" gutterBottom>
          Manage your rooms: lock, hide, set default, or duplicate
        </Typography>

        {/* Search Bar */}
        <Box sx={{ mt: 3, mb: 2 }}>
          <TextField
            fullWidth
            size="small"
            label="Search rooms"
            placeholder="Search by metadata (file path, etc.) - supports regex"
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
            variant="outlined"
          />
        </Box>

        <Box sx={{ height: 600, width: "100%", mt: 3 }}>
          <DataGrid
            rows={rooms}
            columns={columns}
            initialState={{
              pagination: {
                paginationModel: { page: 0, pageSize: 10 },
              },
              sorting: {
                sortModel: [{ field: "isDefault", sort: "desc" }],
              },
            }}
            pageSizeOptions={[5, 10, 25]}
            processRowUpdate={async (newRow, oldRow) => {
              if (newRow.description !== oldRow.description) {
                await handleUpdateDescription(newRow.id, newRow.description);
              }
              return newRow;
            }}
            onProcessRowUpdateError={(error) => {
              showSnackbar("Failed to update row", "error");
            }}
          />
        </Box>

        <Box sx={{ mt: 2, display: "flex", gap: 2 }}>
          <Button
            variant="outlined"
            onClick={() => {
              const roomId = crypto.randomUUID();
              navigate(`/rooms/${roomId}?template=empty`);
            }}
          >
            Create New Empty Room
          </Button>
          <Button
            variant="contained"
            startIcon={<UploadFileIcon />}
            onClick={handleFileUploadClick}
          >
            Upload File
          </Button>
          {fileBrowserEnabled && (
            <Button
              variant="outlined"
              startIcon={<FolderOpenIcon />}
              onClick={() => navigate("/file-browser")}
            >
              Open File Browser
            </Button>
          )}
        </Box>
      </Box>

      {/* Duplicate Dialog */}
      <Dialog 
        open={duplicateDialog.open} 
        onClose={handleCloseDuplicateDialog}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Duplicate Room</DialogTitle>
        <DialogContent>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 3 }}>
            Duplicating: {duplicateDialog.roomDescription || duplicateDialog.roomId}
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

    </Container>
  );
}
