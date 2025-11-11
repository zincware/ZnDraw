import React from "react";
import { Box, Select, MenuItem, FormControl, InputLabel, Chip, Typography } from "@mui/material";
import { FilesystemInfo } from "../../myapi/client";

interface FilesystemSelectorProps {
  filesystems: FilesystemInfo[];
  selected: string;
  onSelect: (name: string) => void;
}

/**
 * Dropdown selector for choosing a filesystem.
 * Shows filesystem type and public status as chips.
 * Uses composite key "name:public" to uniquely identify filesystems across namespaces.
 */
export function FilesystemSelector({ filesystems, selected, onSelect }: FilesystemSelectorProps) {
  if (filesystems.length === 0) {
    return null;
  }

  // If only one filesystem, show info instead of dropdown
  if (filesystems.length === 1) {
    const fs = filesystems[0];
    return (
      <Box sx={{ mb: 2 }}>
        <Typography variant="body2" color="text.secondary">
          <strong>Filesystem:</strong> {fs.name} ({fs.fsType})
          {fs.public && " • Public"}
        </Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ mb: 3 }}>
      <FormControl fullWidth>
        <InputLabel>Select Filesystem</InputLabel>
        <Select
          value={selected}
          label="Select Filesystem"
          onChange={(e) => onSelect(e.target.value)}
        >
          {filesystems.map((fs) => {
            // Create composite key: "name:public"
            const compositeKey = `${fs.name}:${fs.public}`;
            return (
              <MenuItem key={fs.sessionId} value={compositeKey}>
                <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
                  <Typography>{fs.name}</Typography>
                  <Chip label={fs.fsType} size="small" color="primary" variant="outlined" />
                  {fs.public && <Chip label="Public" size="small" color="success" />}
                </Box>
              </MenuItem>
            );
          })}
        </Select>
      </FormControl>
    </Box>
  );
}
