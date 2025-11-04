import React from "react";
import {
  List,
  ListItem,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Divider,
  CircularProgress,
  Alert,
  Box,
} from "@mui/material";
import FolderIcon from "@mui/icons-material/Folder";
import InsertDriveFileIcon from "@mui/icons-material/InsertDriveFile";
import { FilesystemFileItem } from "../../myapi/client";

interface FileListProps {
  files: FilesystemFileItem[];
  currentPath: string;
  isLoading: boolean;
  error: Error | null;
  onItemClick: (item: FilesystemFileItem) => void;
  onGoToParent: () => void;
}

/**
 * List of files and directories in a filesystem.
 * Supports navigation and file selection.
 */
export function FileList({
  files,
  currentPath,
  isLoading,
  error,
  onItemClick,
  onGoToParent,
}: FileListProps) {
  if (isLoading) {
    return (
      <Box sx={{ display: "flex", justifyContent: "center", p: 4 }}>
        <CircularProgress />
      </Box>
    );
  }

  if (error) {
    return (
      <Alert severity="error" sx={{ mb: 2 }}>
        {(error as any)?.response?.data?.error || "Failed to load files"}
      </Alert>
    );
  }

  const formatFileSize = (size: number) => {
    return size / 1024 > 1024
      ? `${(size / 1024 / 1024).toFixed(2)} MB`
      : `${(size / 1024).toFixed(2)} KB`;
  };

  return (
    <List>
      {/* Parent directory link */}
      {currentPath && (
        <>
          <ListItem disablePadding>
            <ListItemButton onClick={onGoToParent}>
              <ListItemIcon>
                <FolderIcon />
              </ListItemIcon>
              <ListItemText primary=".." secondary="Parent directory" />
            </ListItemButton>
          </ListItem>
          <Divider />
        </>
      )}

      {/* Empty directory message */}
      {files.length === 0 && (
        <ListItem>
          <ListItemText primary="Empty directory" secondary="No files or folders to display" />
        </ListItem>
      )}

      {/* Items */}
      {files.map((item, index) => (
        <ListItem key={index} disablePadding>
          <ListItemButton onClick={() => onItemClick(item)}>
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
                      const sizeStr = formatFileSize(item.size);
                      const modifiedStr = item.modified
                        ? ` • Modified: ${new Date(item.modified * 1000).toLocaleDateString()}`
                        : "";
                      return `${sizeStr}${modifiedStr}`;
                    })()
                  : "Directory"
              }
            />
          </ListItemButton>
        </ListItem>
      ))}
    </List>
  );
}
