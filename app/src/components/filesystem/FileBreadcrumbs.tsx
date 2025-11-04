import React from "react";
import { Box, Breadcrumbs, Link, IconButton, Tooltip } from "@mui/material";
import HomeIcon from "@mui/icons-material/Home";

interface FileBreadcrumbsProps {
  currentPath: string;
  onNavigate: (path: string) => void;
  onGoToRoot: () => void;
}

/**
 * Breadcrumb navigation for filesystem paths.
 * Shows clickable path segments and a home button.
 */
export function FileBreadcrumbs({ currentPath, onNavigate, onGoToRoot }: FileBreadcrumbsProps) {
  const pathParts = currentPath ? currentPath.split("/").filter(Boolean) : [];

  return (
    <Box sx={{ mb: 3, display: "flex", alignItems: "center", gap: 1 }}>
      <Tooltip title="Go to root">
        <IconButton size="small" onClick={onGoToRoot}>
          <HomeIcon />
        </IconButton>
      </Tooltip>
      <Breadcrumbs>
        <Link component="button" variant="body1" onClick={onGoToRoot} sx={{ cursor: "pointer" }}>
          Root
        </Link>
        {pathParts.map((part, index) => {
          const path = pathParts.slice(0, index + 1).join("/");
          return (
            <Link
              key={index}
              component="button"
              variant="body1"
              onClick={() => onNavigate(path)}
              sx={{ cursor: "pointer" }}
            >
              {part}
            </Link>
          );
        })}
      </Breadcrumbs>
    </Box>
  );
}
