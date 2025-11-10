import { Paper, Box, Typography, CircularProgress, LinearProgress, IconButton, Fade } from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import { useState } from "react";
import { useAppStore } from "../store";
import { LAYOUT_CONSTANTS } from "../constants/layout";

export default function TaskNotifications() {
  const tasks = useAppStore((state) => state.tasks);
  const roomId = useAppStore((state) => state.roomId);
  const [hiddenTaskIds, setHiddenTaskIds] = useState<Set<string>>(new Set());

  // Filter out hidden tasks and tasks from other rooms
  const taskList = Object.values(tasks).filter(
    task => !hiddenTaskIds.has(task.taskId) && task.roomId === roomId
  );

  const handleHideTask = (taskId: string) => {
    setHiddenTaskIds(prev => new Set(prev).add(taskId));
  };

  if (taskList.length === 0) {
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
      {taskList.map((task) => (
        <Fade key={task.taskId} in={true}>
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
              onClick={() => handleHideTask(task.taskId)}
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

            {/* Task description */}
            <Box sx={{ pr: 4 }}>
              <Typography variant="body1" fontWeight={600}>
                {task.description}
              </Typography>
            </Box>

            {/* Progress indicator */}
            {task.progress !== null ? (
              <Box sx={{ display: 'flex', flexDirection: 'column', gap: 0.5 }}>
                <LinearProgress
                  variant="determinate"
                  value={task.progress}
                  sx={{ height: 6, borderRadius: 3 }}
                />
                <Typography variant="caption" color="text.secondary" textAlign="right">
                  {Math.round(task.progress)}%
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
