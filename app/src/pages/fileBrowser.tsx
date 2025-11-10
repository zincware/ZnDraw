import React, { useState } from "react";
import { useNavigate } from "react-router-dom";
import { useQuery, useMutation } from "@tanstack/react-query";
import {
  Box,
  Container,
  Typography,
  Paper,
  Breadcrumbs,
  Link,
  CircularProgress,
  Alert,
  Button,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  IconButton,
  Tooltip,
  List,
  ListItem,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Divider,
  AppBar,
  Toolbar,
} from "@mui/material";
import FolderIcon from "@mui/icons-material/Folder";
import InsertDriveFileIcon from "@mui/icons-material/InsertDriveFile";
import CheckCircleIcon from "@mui/icons-material/CheckCircle";
import HomeIcon from "@mui/icons-material/Home";
import ArrowBackIcon from "@mui/icons-material/ArrowBack";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import {
  listDirectory,
  loadFile,
  createRoomFromFile,
  DirectoryListResponse,
  FileItem,
  LoadFileRequest,
  LoadFileAlreadyLoadedResponse,
} from "../myapi/client";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import DropOverlay from "../components/DropOverlay";
import { useAppStore } from "../store";

/**
 * FileBrowser page allows browsing local filesystem and loading files into ZnDraw.
 */
export default function FileBrowserPage() {
  const navigate = useNavigate();
  // Use individual selectors to prevent unnecessary re-renders
  const showSnackbar = useAppStore((state) => state.showSnackbar);
  const [currentPath, setCurrentPath] = useState<string>("");
  const [searchQuery, setSearchQuery] = useState<string>("");
  const [loadDialog, setLoadDialog] = useState<{
    open: boolean;
    file: FileItem | null;
  }>({ open: false, file: null });
  const [fileAlreadyLoadedDialog, setFileAlreadyLoadedDialog] = useState<{
    open: boolean;
    data: LoadFileAlreadyLoadedResponse | null;
    filePath: string;
  }>({ open: false, data: null, filePath: "" });
  const [roomName, setRoomName] = useState<string>("");
  const [sliceParams, setSliceParams] = useState<{
    start: string;
    stop: string;
    step: string;
  }>({ start: "", stop: "", step: "" });

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

  // Query for directory listing
  const {
    data: directoryData,
    isLoading,
    error,
    refetch,
  } = useQuery<DirectoryListResponse>({
    queryKey: ["directory", currentPath, searchQuery],
    queryFn: () => listDirectory(currentPath || undefined, searchQuery || undefined),
    retry: false,
  });

  // Mutation for loading files
  const loadFileMutation = useMutation({
    mutationFn: (request: LoadFileRequest) => loadFile(request),
    onSuccess: (data) => {
      if (data.status === "file_already_loaded") {
        // File already loaded - show dialog with options
        const filePath = currentPath
          ? `${currentPath}/${loadDialog.file?.name}`
          : loadDialog.file?.name || "";
        setFileAlreadyLoadedDialog({
          open: true,
          data: data,
          filePath: filePath,
        });
        setLoadDialog({ open: false, file: null });
      } else {
        // File loading queued
        showSnackbar(`File loading queued in room: ${data.room}`, "success");
        setLoadDialog({ open: false, file: null });

        // Navigate directly to the room
        navigate(`/rooms/${data.room}`);
      }
    },
    onError: (error: any) => {
      showSnackbar(error?.response?.data?.error || "Failed to load file", "error");
    },
  });

  // Mutation for creating room from existing file
  const createRoomMutation = useMutation({
    mutationFn: createRoomFromFile,
    onSuccess: (data) => {
      showSnackbar(`New room '${data.roomId}' created from existing file (no re-upload!)`, "success");
      setFileAlreadyLoadedDialog({ open: false, data: null, filePath: "" });

      // Navigate to the new room
      navigate(`/rooms/${data.roomId}`);
    },
    onError: (error: any) => {
      showSnackbar(error?.response?.data?.error || "Failed to create room from file", "error");
    },
  });

  const handleItemClick = (item: FileItem) => {
    if (item.type === "directory") {
      // Navigate to subdirectory
      const newPath = currentPath ? `${currentPath}/${item.name}` : item.name;
      setCurrentPath(newPath);
    } else {
      // Open load dialog for all files
      // Backend will attempt to read unknown formats via ASE
      setLoadDialog({ open: true, file: item });
      setRoomName(""); // Reset room name
      setSliceParams({ start: "", stop: "", step: "" }); // Reset slice params
    }
  };

  const handleBreadcrumbClick = (index: number) => {
    const pathParts = currentPath.split("/").filter(Boolean);
    const newPath = pathParts.slice(0, index + 1).join("/");
    setCurrentPath(newPath);
  };

  const handleGoToRoot = () => {
    setCurrentPath("");
  };

  const handleGoToParent = () => {
    if (directoryData?.parent !== null) {
      setCurrentPath(directoryData?.parent || "");
    }
  };

  const validateSliceParams = (): string | null => {
    const start = sliceParams.start ? parseInt(sliceParams.start) : null;
    const stop = sliceParams.stop ? parseInt(sliceParams.stop) : null;
    const step = sliceParams.step ? parseInt(sliceParams.step) : null;

    if (sliceParams.start && isNaN(start!)) {
      return "Start must be a valid integer";
    }
    if (sliceParams.stop && isNaN(stop!)) {
      return "Stop must be a valid integer";
    }
    if (sliceParams.step && isNaN(step!)) {
      return "Step must be a valid integer";
    }
    if (step !== null && step <= 0) {
      return "Step must be a positive integer (e.g., 1, 2, 5)";
    }
    if (start !== null && start < 0) {
      return "Start cannot be negative";
    }
    if (stop !== null && stop < 0) {
      return "Stop cannot be negative";
    }
    if (start !== null && stop !== null && start >= stop) {
      return "Start must be less than stop";
    }
    return null;
  };

  const handleLoadFile = () => {
    if (!loadDialog.file) return;

    // Validate slice parameters
    const validationError = validateSliceParams();
    if (validationError) {
      showSnackbar(validationError, "error");
      return;
    }

    const filePath = currentPath
      ? `${currentPath}/${loadDialog.file.name}`
      : loadDialog.file.name;

    const request: LoadFileRequest = {
      path: filePath,
      room: roomName || undefined,
      // Only include slice params if they're non-empty
      ...(sliceParams.start && { start: parseInt(sliceParams.start) }),
      ...(sliceParams.stop && { stop: parseInt(sliceParams.stop) }),
      ...(sliceParams.step && { step: parseInt(sliceParams.step) }),
    };

    loadFileMutation.mutate(request);
  };

  const handleOpenExistingRoom = () => {
    if (!fileAlreadyLoadedDialog.data) return;
    navigate(`/rooms/${fileAlreadyLoadedDialog.data.existingRoom}`);
    setFileAlreadyLoadedDialog({ open: false, data: null, filePath: "" });
  };

  const handleCreateNewRoom = () => {
    if (!fileAlreadyLoadedDialog.data) return;
    createRoomMutation.mutate({
      sourceRoom: fileAlreadyLoadedDialog.data.existingRoom,
    });
  };

  const handleForceUpload = () => {
    if (!fileAlreadyLoadedDialog.data) return;
    
    const request: LoadFileRequest = {
      path: fileAlreadyLoadedDialog.filePath,
      force_upload: true,
    };

    setFileAlreadyLoadedDialog({ open: false, data: null, filePath: "" });
    loadFileMutation.mutate(request);
  };

  const handleBack = () => {
    navigate("/");
  };

  // Parse breadcrumbs from current path
  const pathParts = currentPath ? currentPath.split("/").filter(Boolean) : [];

  // Show error state
  if (error) {
    return (
      <Container maxWidth="lg" sx={{ mt: 4 }}>
        <Alert severity="error">
          {(error as any)?.response?.data?.error || "Failed to load directory"}
        </Alert>
        <Button onClick={() => navigate("/")} sx={{ mt: 2 }}>
          Back to Home
        </Button>
      </Container>
    );
  }

  return (
    <Box
      sx={{ flexGrow: 1 }}
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

      <AppBar position="static">
        <Toolbar>
          <IconButton
            edge="start"
            color="inherit"
            onClick={handleBack}
            sx={{ mr: 2 }}
          >
            <ArrowBackIcon />
          </IconButton>
          <Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
            File Browser
          </Typography>
          <Button
            color="inherit"
            startIcon={<UploadFileIcon />}
            onClick={handleFileUploadClick}
          >
            Upload File
          </Button>
        </Toolbar>
      </AppBar>

      <Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
        <Paper sx={{ p: 3 }}>
          {/* Breadcrumbs */}
          <Box sx={{ mb: 3, display: "flex", alignItems: "center", gap: 1 }}>
            <Tooltip title="Go to root">
              <IconButton size="small" onClick={handleGoToRoot}>
                <HomeIcon />
              </IconButton>
            </Tooltip>
            <Breadcrumbs>
              <Link
                component="button"
                variant="body1"
                onClick={handleGoToRoot}
                sx={{ cursor: "pointer" }}
              >
                Root
              </Link>
              {pathParts.map((part, index) => (
                <Link
                  key={index}
                  component="button"
                  variant="body1"
                  onClick={() => handleBreadcrumbClick(index)}
                  sx={{ cursor: "pointer" }}
                >
                  {part}
                </Link>
              ))}
            </Breadcrumbs>
          </Box>

          {/* Search Bar */}
          <Box sx={{ mb: 2 }}>
            <TextField
              fullWidth
              size="small"
              label="Search files and folders"
              placeholder="Filter by name (supports regex)"
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              variant="outlined"
            />
          </Box>

          <Divider sx={{ mb: 2 }} />

          {/* Loading state */}
          {isLoading && (
            <Box sx={{ display: "flex", justifyContent: "center", p: 4 }}>
              <CircularProgress />
            </Box>
          )}

          {/* Directory listing */}
          {!isLoading && directoryData && (
            <List>
              {/* Parent directory link */}
              {directoryData.parent !== null && (
                <>
                  <ListItem disablePadding>
                    <ListItemButton onClick={handleGoToParent}>
                      <ListItemIcon>
                        <FolderIcon />
                      </ListItemIcon>
                      <ListItemText primary=".." secondary="Parent directory" />
                    </ListItemButton>
                  </ListItem>
                  <Divider />
                </>
              )}

              {/* Items */}
              {directoryData.items.length === 0 && (
                <ListItem>
                  <ListItemText
                    primary="Empty directory"
                    secondary="No files or folders to display"
                  />
                </ListItem>
              )}

              {directoryData.items.map((item, index) => (
                <ListItem
                  key={index}
                  disablePadding
                  secondaryAction={
                    item.type === "file" ? (
                      item.alreadyLoaded ? (
                        <Tooltip title={`Already loaded in room: ${item.alreadyLoaded.room}`}>
                          <CheckCircleIcon color="primary" />
                        </Tooltip>
                      ) : item.supported ? (
                        <Tooltip title={item.format_info || "Supported file type"}>
                          <CheckCircleIcon color="success" />
                        </Tooltip>
                      ) : null
                    ) : null
                  }
                >
                  <ListItemButton onClick={() => handleItemClick(item)}>
                    <ListItemIcon>
                      {item.type === "directory" ? (
                        <FolderIcon color="primary" />
                      ) : (
                        <InsertDriveFileIcon color="action" />
                      )}
                    </ListItemIcon>
                    <ListItemText
                      primary={item.name}
                      secondary={
                        item.type === "file"
                          ? (() => {
                              const sizeStr = (item.size || 0) / 1024 > 1024
                                ? `${((item.size || 0) / 1024 / 1024).toFixed(2)} MB`
                                : `${((item.size || 0) / 1024).toFixed(2)} KB`;
                              const formatStr = item.format_info ? ` • ${item.format_info}` : "";
                              const loadedStr = item.alreadyLoaded
                                ? ` • Already loaded in '${item.alreadyLoaded.room}'`
                                : "";
                              return `${sizeStr}${formatStr}${loadedStr}`;
                            })()
                          : "Directory"
                      }
                    />
                  </ListItemButton>
                </ListItem>
              ))}
            </List>
          )}
        </Paper>
      </Container>

      {/* Load File Dialog */}
      <Dialog
        open={loadDialog.open}
        onClose={() => setLoadDialog({ open: false, file: null })}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>Load File into ZnDraw</DialogTitle>
        <DialogContent>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            File: {loadDialog.file?.name}
          </Typography>
          {loadDialog.file?.format_info && (
            <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
              Format: {loadDialog.file.format_info}
            </Typography>
          )}
          {!loadDialog.file?.format_info && (
            <Typography variant="body2" color="warning.main" sx={{ mb: 2 }}>
              Unknown format - will attempt to read via ASE
            </Typography>
          )}

          <TextField
            autoFocus
            margin="dense"
            label="Room Name (optional)"
            type="text"
            fullWidth
            variant="outlined"
            value={roomName}
            onChange={(e) => setRoomName(e.target.value)}
            helperText="Leave empty to auto-generate a room name from the filename"
            sx={{ mb: 2 }}
          />

          <Typography variant="body2" sx={{ mb: 1, fontWeight: 500 }}>
            Frame Selection (optional)
          </Typography>
          <Typography variant="caption" color="text.secondary" sx={{ mb: 2, display: "block" }}>
            Load specific frames from trajectory files. Leave empty to load all frames.
          </Typography>

          <Box sx={{ display: "flex", gap: 2 }}>
            <TextField
              size="small"
              label="Start"
              type="number"
              value={sliceParams.start}
              onChange={(e) =>
                setSliceParams({ ...sliceParams, start: e.target.value })
              }
              helperText="First frame"
              sx={{ flex: 1 }}
              slotProps={{ htmlInput: { min: 0 } }}
            />
            <TextField
              size="small"
              label="Stop"
              type="number"
              value={sliceParams.stop}
              onChange={(e) =>
                setSliceParams({ ...sliceParams, stop: e.target.value })
              }
              helperText="Last frame"
              sx={{ flex: 1 }}
              slotProps={{ htmlInput: { min: 0 } }}
            />
            <TextField
              size="small"
              label="Step"
              type="number"
              value={sliceParams.step}
              onChange={(e) =>
                setSliceParams({ ...sliceParams, step: e.target.value })
              }
              helperText="Interval"
              sx={{ flex: 1 }}
              slotProps={{ htmlInput: { min: 1 } }}
            />
          </Box>

          {(sliceParams.start || sliceParams.stop || sliceParams.step) && (
            <Typography
              variant="caption"
              color="primary.main"
              sx={{ mt: 1, display: "block" }}
            >
              Will load: [{sliceParams.start || "0"}:{sliceParams.stop || "end"}:
              {sliceParams.step || "1"}]
              {sliceParams.step && parseInt(sliceParams.step) > 1 &&
                ` (every ${sliceParams.step} frame${parseInt(sliceParams.step) > 1 ? "s" : ""})`}
            </Typography>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={() => setLoadDialog({ open: false, file: null })}>
            Cancel
          </Button>
          <Button
            onClick={handleLoadFile}
            variant="contained"
            disabled={loadFileMutation.isPending}
          >
            {loadFileMutation.isPending ? "Loading..." : "Load"}
          </Button>
        </DialogActions>
      </Dialog>

      {/* File Already Loaded Dialog */}
      <Dialog
        open={fileAlreadyLoadedDialog.open}
        onClose={() => setFileAlreadyLoadedDialog({ open: false, data: null, filePath: "" })}
        maxWidth="sm"
        fullWidth
      >
        <DialogTitle>File Already Loaded</DialogTitle>
        <DialogContent>
          <Typography variant="body1" sx={{ mb: 2 }}>
            {fileAlreadyLoadedDialog.data?.message}
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
            This exact file (same size and modification time) is already loaded in room{" "}
            <strong>{fileAlreadyLoadedDialog.data?.existingRoom}</strong>.
          </Typography>
          <Typography variant="body2" color="text.secondary">
            What would you like to do?
          </Typography>
        </DialogContent>
        <DialogActions sx={{ flexDirection: "column", gap: 1, alignItems: "stretch", p: 2 }}>
          <Button
            onClick={handleOpenExistingRoom}
            variant="contained"
            color="primary"
            fullWidth
          >
            Open Existing Room
          </Button>
          <Button
            onClick={handleCreateNewRoom}
            variant="contained"
            color="success"
            fullWidth
            disabled={createRoomMutation.isPending}
          >
            {createRoomMutation.isPending ? "Creating..." : "Create New Room (Reuse Storage - Fast!)"}
          </Button>
          <Button
            onClick={handleForceUpload}
            variant="outlined"
            color="warning"
            fullWidth
            disabled={loadFileMutation.isPending}
          >
            {loadFileMutation.isPending ? "Uploading..." : "Upload Anyway (Ignore Existing)"}
          </Button>
          <Button
            onClick={() => setFileAlreadyLoadedDialog({ open: false, data: null, filePath: "" })}
            variant="text"
            fullWidth
          >
            Cancel
          </Button>
        </DialogActions>
      </Dialog>
    </Box>
  );
}
