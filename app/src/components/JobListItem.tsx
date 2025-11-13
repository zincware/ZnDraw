import { Box, Typography, Chip, Stack, Paper, Collapse, IconButton } from "@mui/material";
import { useState } from "react";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import PendingIcon from "@mui/icons-material/HourglassEmpty";
import AssignmentIcon from "@mui/icons-material/Assignment";
import PlayArrowIcon from "@mui/icons-material/PlayArrow";
import CheckCircleIcon from "@mui/icons-material/CheckCircle";
import ErrorIcon from "@mui/icons-material/Error";
import { type Job } from "../hooks/useSchemas";
import { getJobStatusColor, getJobStatusLabel, calculateJobDurations, formatDuration } from "../utils/jobUtils";

interface JobListItemProps {
  job: Job;
  showExtensionName?: boolean;
}

/**
 * Displays individual job as a slim one-liner with status, time, and ID
 */
export const JobListItem = ({ job, showExtensionName = false }: JobListItemProps) => {
  const [expanded, setExpanded] = useState(false);
  const durations = calculateJobDurations(job);
  const statusColor = getJobStatusColor(job.status);
  const statusLabel = getJobStatusLabel(job.status);

  const getStatusIcon = () => {
    switch (job.status) {
      case "pending":
        return <PendingIcon fontSize="small" />;
      case "assigned":
        return <AssignmentIcon fontSize="small" />;
      case "processing":
        return <PlayArrowIcon fontSize="small" />;
      case "completed":
        return <CheckCircleIcon fontSize="small" />;
      case "failed":
        return <ErrorIcon fontSize="small" />;
      default:
        return <PendingIcon fontSize="small" />;
    }
  };

  // Get time display: runtime for finished jobs, time since created for others
  const getTimeDisplay = () => {
    if (job.status === "completed" || job.status === "failed") {
      return durations.total ? formatDuration(durations.total) : "-";
    }

    // For pending/assigned/processing: show time since created
    const createdAt = new Date(job.created_at);
    const now = new Date();
    const elapsed = Math.floor((now.getTime() - createdAt.getTime()) / 1000);
    return formatDuration(elapsed);
  };

  return (
    <Paper
      elevation={1}
      sx={{
        p: 1.5,
        cursor: "pointer",
        "&:hover": {
          bgcolor: "action.hover"
        }
      }}
      onClick={() => setExpanded(!expanded)}
    >
      <Stack direction="row" alignItems="center" spacing={1.5}>
        <Chip
          icon={getStatusIcon()}
          label={statusLabel}
          color={statusColor}
          size="small"
          sx={{ minWidth: 100 }}
        />

        <Typography variant="body2" fontWeight="medium">
          {getTimeDisplay()}
        </Typography>

        <Typography variant="body2" color="text.secondary" sx={{ fontFamily: "monospace", fontSize: "0.85rem" }}>
          {job.id.substring(0, 8)}
        </Typography>

        <Box sx={{ flex: 1 }} />

        <IconButton
          size="small"
          onClick={(e) => {
            e.stopPropagation();
            setExpanded(!expanded);
          }}
        >
          {expanded ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
        </IconButton>
      </Stack>

      {/* Expanded details */}
      <Collapse in={expanded}>
        <Box sx={{ mt: 2, pt: 2, borderTop: 1, borderColor: "divider" }}>
          <Stack spacing={1}>
            <Box>
              <Typography variant="caption" color="text.secondary">
                Job ID:
              </Typography>
              <Typography variant="body2" fontFamily="monospace">
                {job.id}
              </Typography>
            </Box>

            {job.worker_id && (
              <Box>
                <Typography variant="caption" color="text.secondary">
                  Worker:
                </Typography>
                <Typography variant="body2" fontFamily="monospace">
                  {job.worker_id}
                </Typography>
              </Box>
            )}

            <Box>
              <Typography variant="caption" color="text.secondary">
                Provider:
              </Typography>
              <Typography variant="body2">
                {job.provider}
              </Typography>
            </Box>

            <Box>
              <Typography variant="caption" color="text.secondary">
                Created:
              </Typography>
              <Typography variant="body2">
                {new Date(job.created_at).toLocaleString()}
              </Typography>
            </Box>

            {job.assigned_at && (
              <Box>
                <Typography variant="caption" color="text.secondary">
                  Assigned:
                </Typography>
                <Typography variant="body2">
                  {new Date(job.assigned_at).toLocaleString()}
                </Typography>
              </Box>
            )}

            {job.started_at && (
              <Box>
                <Typography variant="caption" color="text.secondary">
                  Started:
                </Typography>
                <Typography variant="body2">
                  {new Date(job.started_at).toLocaleString()}
                </Typography>
              </Box>
            )}

            {job.completed_at && (
              <Box>
                <Typography variant="caption" color="text.secondary">
                  Completed:
                </Typography>
                <Typography variant="body2">
                  {new Date(job.completed_at).toLocaleString()}
                </Typography>
              </Box>
            )}

            {job.error && (
              <Box>
                <Typography variant="caption" color="error">
                  Error:
                </Typography>
                <Typography variant="body2" color="error">
                  {job.error}
                </Typography>
              </Box>
            )}

            {/* Timing information */}
            {(durations.waiting || durations.running) && (
              <Box sx={{ mt: 1 }}>
                <Typography variant="caption" color="text.secondary" display="block" sx={{ mb: 0.5 }}>
                  Timing:
                </Typography>
                <Stack direction="row" spacing={2}>
                  {durations.waiting && (
                    <Typography variant="body2">
                      <strong>Wait:</strong> {formatDuration(durations.waiting)}
                    </Typography>
                  )}
                  {durations.running && (
                    <Typography variant="body2">
                      <strong>Run:</strong> {formatDuration(durations.running)}
                    </Typography>
                  )}
                  {durations.total && (
                    <Typography variant="body2">
                      <strong>Total:</strong> {formatDuration(durations.total)}
                    </Typography>
                  )}
                </Stack>
              </Box>
            )}
          </Stack>
        </Box>
      </Collapse>
    </Paper>
  );
};
