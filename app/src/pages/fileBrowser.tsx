import { useState, useEffect } from "react";
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
  Snackbar,
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
import {
  listDirectory,
  loadFile,
  DirectoryListResponse,
  FileItem,
  LoadFileRequest,
} from "../myapi/client";

/**
 * FileBrowser page allows browsing local filesystem and loading files into ZnDraw.
 */
export default function FileBrowserPage() {
  const navigate = useNavigate();
  const [currentPath, setCurrentPath] = useState<string>("");
  const [loadDialog, setLoadDialog] = useState<{
    open: boolean;
    file: FileItem | null;
  }>({ open: false, file: null });
  const [roomName, setRoomName] = useState<string>("");
  const [snackbar, setSnackbar] = useState<{
    open: boolean;
    message: string;
    severity: "success" | "error";
  }>({ open: false, message: "", severity: "success" });

  // Query for directory listing
  const {
    data: directoryData,
    isLoading,
    error,
    refetch,
  } = useQuery<DirectoryListResponse>({
    queryKey: ["directory", currentPath],
    queryFn: () => listDirectory(currentPath || undefined),
    retry: false,
  });

  // Mutation for loading files
  const loadFileMutation = useMutation({
    mutationFn: (request: LoadFileRequest) => loadFile(request),
    onSuccess: (data) => {
      setSnackbar({
        open: true,
        message: `File loading queued in room: ${data.room}`,
        severity: "success",
      });
      setLoadDialog({ open: false, file: null });

      // Navigate to the room with waitForCreation flag
      const userId = crypto.randomUUID();
      navigate(`/rooms/${data.room}/${userId}?waitForCreation=true`);
    },
    onError: (error: any) => {
      setSnackbar({
        open: true,
        message: error?.response?.data?.error || "Failed to load file",
        severity: "error",
      });
    },
  });

  const handleItemClick = (item: FileItem) => {
    if (item.type === "directory") {
      // Navigate to subdirectory
      const newPath = currentPath ? `${currentPath}/${item.name}` : item.name;
      setCurrentPath(newPath);
    } else if (item.supported) {
      // Open load dialog for supported files
      setLoadDialog({ open: true, file: item });
      setRoomName(""); // Reset room name
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
      setCurrentPath(directoryData.parent || "");
    }
  };

  const handleLoadFile = () => {
    if (!loadDialog.file) return;

    const filePath = currentPath
      ? `${currentPath}/${loadDialog.file.name}`
      : loadDialog.file.name;

    const request: LoadFileRequest = {
      path: filePath,
      room: roomName || undefined,
    };

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
    <Box sx={{ flexGrow: 1 }}>
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
                    item.type === "file" && item.supported ? (
                      <Tooltip title="Supported file type">
                        <CheckCircleIcon color="success" />
                      </Tooltip>
                    ) : null
                  }
                >
                  <ListItemButton
                    onClick={() => handleItemClick(item)}
                    disabled={item.type === "file" && !item.supported}
                  >
                    <ListItemIcon>
                      {item.type === "directory" ? (
                        <FolderIcon color="primary" />
                      ) : (
                        <InsertDriveFileIcon
                          color={item.supported ? "action" : "disabled"}
                        />
                      )}
                    </ListItemIcon>
                    <ListItemText
                      primary={item.name}
                      secondary={
                        item.type === "file"
                          ? `${(item.size || 0) / 1024 > 1024 ? `${((item.size || 0) / 1024 / 1024).toFixed(2)} MB` : `${((item.size || 0) / 1024).toFixed(2)} KB`}`
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
          <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
            File: {loadDialog.file?.name}
          </Typography>

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
          />
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
    </Box>
  );
}
