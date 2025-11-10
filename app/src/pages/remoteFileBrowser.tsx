import React, { useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { useQuery } from "@tanstack/react-query";
import {
  Box,
  Container,
  Typography,
  Paper,
  CircularProgress,
  Alert,
  Button,
  IconButton,
  AppBar,
  Toolbar,
  Divider,
} from "@mui/material";
import ArrowBackIcon from "@mui/icons-material/ArrowBack";
import CloudIcon from "@mui/icons-material/Cloud";
import {
  listFilesystems,
  listFilesystemFiles,
  FilesystemFileItem,
  LoadFilesystemFileRequest,
} from "../myapi/client";
import {
  FilesystemSelector,
  FileBreadcrumbs,
  FileList,
  LoadFileDialog,
} from "../components/filesystem";

/**
 * RemoteFileBrowser page allows browsing registered remote filesystems and loading files into ZnDraw.
 *
 * This page provides a UI for:
 * - Selecting from available remote filesystems
 * - Browsing directories and files
 * - Loading files into ZnDraw rooms with optional frame slicing
 */
export default function RemoteFileBrowserPage() {
  const navigate = useNavigate();
  const { roomId } = useParams<{ roomId: string }>();

  const [selectedFilesystem, setSelectedFilesystem] = useState<string>("");
  const [currentPath, setCurrentPath] = useState<string>("");
  const [loadDialog, setLoadDialog] = useState<{
    open: boolean;
    file: FilesystemFileItem | null;
  }>({ open: false, file: null });

  // Query for available filesystems
  const {
    data: filesystemsData,
    isLoading: isLoadingFilesystems,
    error: filesystemsError,
  } = useQuery({
    queryKey: ["filesystems", roomId],
    queryFn: () => listFilesystems(roomId!),
    enabled: !!roomId,
    retry: false,
  });

  // Auto-select first filesystem if only one is available
  React.useEffect(() => {
    if (
      filesystemsData &&
      filesystemsData.filesystems.length === 1 &&
      !selectedFilesystem
    ) {
      setSelectedFilesystem(filesystemsData.filesystems[0].name);
    }
  }, [filesystemsData, selectedFilesystem]);

  // Query for directory listing
  const {
    data: filesData,
    isLoading: isLoadingFiles,
    error: filesError,
  } = useQuery({
    queryKey: ["filesystemFiles", roomId, selectedFilesystem, currentPath],
    queryFn: () =>
      listFilesystemFiles(roomId!, selectedFilesystem, currentPath || undefined),
    enabled: !!roomId && !!selectedFilesystem,
    retry: false,
  });

  // Navigate to room immediately with load request state
  const handleLoadFile = (request: LoadFilesystemFileRequest & { targetRoom: string }) => {
    const targetRoom = request.targetRoom || roomId;

    // Close dialog immediately
    setLoadDialog({ open: false, file: null });

    // Navigate to target room with pending load state
    navigate(`/rooms/${targetRoom}`, {
      state: {
        pendingFilesystemLoad: {
          fsName: selectedFilesystem,
          request,
        },
      },
    });
  };

  // Navigation handlers
  const handleItemClick = (item: FilesystemFileItem) => {
    if (item.type === "directory") {
      const newPath = currentPath ? `${currentPath}/${item.name}` : item.name;
      setCurrentPath(newPath);
    } else {
      setLoadDialog({ open: true, file: item });
    }
  };

  const handleNavigate = (path: string) => {
    setCurrentPath(path);
  };

  const handleGoToRoot = () => {
    setCurrentPath("");
  };

  const handleGoToParent = () => {
    const pathParts = currentPath.split("/").filter(Boolean);
    if (pathParts.length > 0) {
      setCurrentPath(pathParts.slice(0, -1).join("/"));
    }
  };

  const handleBack = () => {
    if (roomId) {
      navigate(`/rooms/${roomId}`);
    } else {
      navigate("/");
    }
  };

  const handleFilesystemSelect = (name: string) => {
    setSelectedFilesystem(name);
    setCurrentPath(""); // Reset path when changing filesystem
  };

  // Show error state for filesystems
  if (filesystemsError) {
    return (
      <Container maxWidth="lg" sx={{ mt: 4 }}>
        <Alert severity="error">
          {(filesystemsError as any)?.response?.data?.error ||
            "Failed to load filesystems"}
        </Alert>
        <Button onClick={handleBack} sx={{ mt: 2 }}>
          Back
        </Button>
      </Container>
    );
  }

  // Show loading state
  if (isLoadingFilesystems) {
    return (
      <Container maxWidth="lg" sx={{ mt: 4 }}>
        <Box sx={{ display: "flex", justifyContent: "center", p: 4 }}>
          <CircularProgress />
        </Box>
      </Container>
    );
  }

  // Show no filesystems message
  if (filesystemsData && filesystemsData.filesystems.length === 0) {
    return (
      <Container maxWidth="lg" sx={{ mt: 4 }}>
        <Alert severity="info">
          No filesystems are currently registered for this room. Register a
          filesystem using the ZnDraw Python client.
        </Alert>
        <Button onClick={handleBack} sx={{ mt: 2 }}>
          Back
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
          <CloudIcon sx={{ mr: 1 }} />
          <Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
            Remote Filesystem Browser
          </Typography>
        </Toolbar>
      </AppBar>

      <Container maxWidth="lg" sx={{ mt: 4, mb: 4 }}>
        <Paper sx={{ p: 3 }}>
          {/* Filesystem Selector */}
          {filesystemsData && (
            <FilesystemSelector
              filesystems={filesystemsData.filesystems}
              selected={selectedFilesystem}
              onSelect={handleFilesystemSelect}
            />
          )}

          {/* Breadcrumbs */}
          {selectedFilesystem && (
            <>
              <FileBreadcrumbs
                currentPath={currentPath}
                onNavigate={handleNavigate}
                onGoToRoot={handleGoToRoot}
              />
              <Divider sx={{ mb: 2 }} />
            </>
          )}

          {/* File List */}
          <FileList
            files={filesData?.files || []}
            currentPath={currentPath}
            isLoading={isLoadingFiles}
            error={filesError as Error | null}
            onItemClick={handleItemClick}
            onGoToParent={handleGoToParent}
          />
        </Paper>
      </Container>

      {/* Load File Dialog */}
      <LoadFileDialog
        open={loadDialog.open}
        file={loadDialog.file}
        defaultRoomId={roomId || ""}
        isLoading={false}
        onClose={() => setLoadDialog({ open: false, file: null })}
        onLoad={handleLoadFile}
      />
    </Box>
  );
}
