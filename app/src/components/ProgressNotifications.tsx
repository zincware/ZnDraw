import { Paper, Box, Typography, CircularProgress, LinearProgress, IconButton, Fade } from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import { useState } from "react";
import { useAppStore } from "../store";
import { LAYOUT_CONSTANTS } from "../constants/layout";

export default function ProgressNotifications() {
  const progressTrackers = useAppStore((state) => state.progressTrackers);
  const roomId = useAppStore((state) => state.roomId);
  const [hiddenProgressIds, setHiddenProgressIds] = useState<Set<string>>(new Set());

  // Filter out hidden progress trackers and trackers from other rooms
  const progressList = Object.values(progressTrackers).filter(
    tracker => !hiddenProgressIds.has(tracker.progressId) && tracker.roomId === roomId
  );

  const handleHideProgress = (progressId: string) => {
    setHiddenProgressIds(prev => new Set(prev).add(progressId));
  };

  if (progressList.length === 0) {
    return null;
  }

  return (
    <Box
      sx={{
        position: 'fixed',
        bottom: LAYOUT_CONSTANTS.PROGRESSBAR_HEIGHT + 24,
        right: 24,
        zIndex: (theme) => theme.zIndex.drawer + 2,
        display: 'flex',
        flexDirection: 'column',
        gap: 1,
        maxWidth: 400,
      }}
    >
      {progressList.map((tracker) => (
        <Fade key={tracker.progressId} in={true}>
          <Paper
            elevation={8}
            sx={{
              p: 2,
              display: 'flex',
              flexDirection: 'column',
              gap: 1.5,
              minWidth: 320,
              position: 'relative',
            }}
          >
            {/* Close button */}
            <IconButton
              size="small"
              onClick={() => handleHideProgress(tracker.progressId)}
              sx={{
                position: 'absolute',
                top: 8,
                right: 8,
                opacity: 0.6,
                '&:hover': {
                  opacity: 1,
                },
              }}
            >
              <CloseIcon fontSize="small" />
            </IconButton>

            {/* Progress description */}
            <Box sx={{ pr: 4 }}>
              <Typography variant="body1" fontWeight={600}>
                {tracker.description}
              </Typography>
            </Box>

            {/* Progress indicator */}
            {tracker.progress !== null ? (
              <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.5 }}>
                <LinearProgress
                  variant="determinate"
                  value={tracker.progress}
                  sx={{ height: 6, borderRadius: 3 }}
                />
                <Typography variant="caption" color="text.secondary" textAlign="right">
                  {Math.round(tracker.progress)}%
                </Typography>
              </Box>
            ) : (
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
                <CircularProgress size={24} />
                <Typography variant="caption" color="text.secondary">
                  In progress...
                </Typography>
              </Box>
            )}
          </Paper>
        </Fade>
      ))}
    </Box>
  );
}
