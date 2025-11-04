import React, { useState } from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  TextField,
  Typography,
  Box,
} from "@mui/material";
import { FilesystemFileItem, LoadFilesystemFileRequest } from "../../myapi/client";

interface LoadFileDialogProps {
  open: boolean;
  file: FilesystemFileItem | null;
  defaultRoomId: string;
  isLoading: boolean;
  onClose: () => void;
  onLoad: (request: LoadFilesystemFileRequest & { targetRoom: string }) => void;
}

/**
 * Dialog for loading a file from filesystem with slice parameters.
 * Allows specifying target room and frame selection.
 */
export function LoadFileDialog({
  open,
  file,
  defaultRoomId,
  isLoading,
  onClose,
  onLoad,
}: LoadFileDialogProps) {
  const [targetRoom, setTargetRoom] = useState<string>(defaultRoomId);
  const [sliceParams, setSliceParams] = useState<{
    start: string;
    stop: string;
    step: string;
  }>({ start: "", stop: "", step: "" });

  // Reset state when dialog opens
  React.useEffect(() => {
    if (open) {
      setTargetRoom(defaultRoomId);
      setSliceParams({ start: "", stop: "", step: "" });
    }
  }, [open, defaultRoomId]);

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

  const handleLoad = () => {
    if (!file) return;

    const validationError = validateSliceParams();
    if (validationError) {
      alert(validationError);
      return;
    }

    const request: LoadFilesystemFileRequest & { targetRoom: string } = {
      path: file.path,
      targetRoom: targetRoom,
      ...(sliceParams.start && { start: parseInt(sliceParams.start) }),
      ...(sliceParams.stop && { stop: parseInt(sliceParams.stop) }),
      ...(sliceParams.step && { step: parseInt(sliceParams.step) }),
    };

    onLoad(request);
  };

  if (!file) return null;

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle>Load File from Remote Filesystem</DialogTitle>
      <DialogContent>
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          File: {file.name}
        </Typography>
        <Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
          Path: {file.path}
        </Typography>

        <TextField
          autoFocus
          margin="dense"
          label="Target Room"
          type="text"
          fullWidth
          variant="outlined"
          value={targetRoom}
          onChange={(e) => setTargetRoom(e.target.value)}
          helperText="Room to load the file into"
          sx={{ mb: 2 }}
        />

        <Typography variant="body2" sx={{ mb: 1, fontWeight: 500 }}>
          Frame Selection (optional)
        </Typography>
        <Typography
          variant="caption"
          color="text.secondary"
          sx={{ mb: 2, display: "block" }}
        >
          Load specific frames from trajectory files. Leave empty to load all frames.
        </Typography>

        <Box sx={{ display: "flex", gap: 2 }}>
          <TextField
            size="small"
            label="Start"
            type="number"
            value={sliceParams.start}
            onChange={(e) => setSliceParams({ ...sliceParams, start: e.target.value })}
            helperText="First frame"
            sx={{ flex: 1 }}
            slotProps={{ htmlInput: { min: 0 } }}
          />
          <TextField
            size="small"
            label="Stop"
            type="number"
            value={sliceParams.stop}
            onChange={(e) => setSliceParams({ ...sliceParams, stop: e.target.value })}
            helperText="Last frame"
            sx={{ flex: 1 }}
            slotProps={{ htmlInput: { min: 0 } }}
          />
          <TextField
            size="small"
            label="Step"
            type="number"
            value={sliceParams.step}
            onChange={(e) => setSliceParams({ ...sliceParams, step: e.target.value })}
            helperText="Interval"
            sx={{ flex: 1 }}
            slotProps={{ htmlInput: { min: 1 } }}
          />
        </Box>

        {(sliceParams.start || sliceParams.stop || sliceParams.step) && (
          <Typography variant="caption" color="primary.main" sx={{ mt: 1, display: "block" }}>
            Will load: [{sliceParams.start || "0"}:{sliceParams.stop || "end"}:
            {sliceParams.step || "1"}]
            {sliceParams.step &&
              parseInt(sliceParams.step) > 1 &&
              ` (every ${sliceParams.step} frame${
                parseInt(sliceParams.step) > 1 ? "s" : ""
              })`}
          </Typography>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Cancel</Button>
        <Button onClick={handleLoad} variant="contained" disabled={isLoading}>
          {isLoading ? "Loading..." : "Load"}
        </Button>
      </DialogActions>
    </Dialog>
  );
}
