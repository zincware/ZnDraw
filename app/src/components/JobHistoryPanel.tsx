import { Box, Typography, Stack, Divider, Button, Skeleton, IconButton, Collapse } from "@mui/material";
import { useMemo, useState } from "react";
import HistoryIcon from "@mui/icons-material/History";
import RefreshIcon from "@mui/icons-material/Refresh";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import { useJobs } from "../hooks/useSchemas";
import { JobListItem } from "./JobListItem";

interface JobHistoryPanelProps {
  roomId: string;
  category: string;
  extension: string | null;
  isPublic: boolean; // Whether this is a public/global extension
  maxJobs?: number;
}

/**
 * Shows recent job history for the current extension
 * - Displays last N jobs for the extension
 * - Shows state transitions with timestamps
 * - Expandable for details
 */
export const JobHistoryPanel = ({
  roomId,
  category,
  extension,
  isPublic,
  maxJobs = 10,
}: JobHistoryPanelProps) => {
  const { data: allJobs, isLoading, refetch } = useJobs(roomId);
  const [expanded, setExpanded] = useState(true);

  // Filter jobs for current extension and namespace (public/private) and limit to maxJobs
  const extensionJobs = useMemo(() => {
    if (!allJobs || !extension) {
      return [];
    }

    // Backend returns public as string "true"/"false"
    const publicStr = isPublic ? "true" : "false";

    const filtered = allJobs.filter(
      (job) =>
        job.category === category &&
        job.extension === extension &&
        job.public === publicStr // Filter by namespace
    );

    return filtered
      .sort((a, b) => {
        // Sort by created_at descending (newest first)
        return new Date(b.created_at).getTime() - new Date(a.created_at).getTime();
      })
      .slice(0, maxJobs);
  }, [allJobs, category, extension, isPublic, maxJobs]);

  if (!extension) {
    return null;
  }

  return (
    <Box sx={{ mt: 2 }}>
      <Divider sx={{ mb: 2 }} />

      <Stack direction="row" alignItems="center" justifyContent="space-between" sx={{ mb: 2 }}>
        <Stack direction="row" alignItems="center" spacing={1}>
          <IconButton
            size="small"
            onClick={() => setExpanded(!expanded)}
            aria-label={expanded ? "collapse recent jobs" : "expand recent jobs"}
          >
            {expanded ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
          </IconButton>
          <HistoryIcon fontSize="small" color="action" />
          <Typography variant="subtitle2" color="text.secondary">
            Recent Jobs
          </Typography>
        </Stack>

        <Button
          size="small"
          startIcon={<RefreshIcon />}
          onClick={() => refetch()}
          disabled={isLoading}
        >
          Refresh
        </Button>
      </Stack>

      <Collapse in={expanded}>
        {isLoading ? (
          <Stack spacing={1}>
            <Skeleton variant="rectangular" height={48} sx={{ borderRadius: 1 }} />
            <Skeleton variant="rectangular" height={48} sx={{ borderRadius: 1 }} />
            <Skeleton variant="rectangular" height={48} sx={{ borderRadius: 1 }} />
          </Stack>
        ) : extensionJobs.length === 0 ? (
          <Box sx={{ textAlign: "center", py: 2 }}>
            <Typography variant="body2" color="text.secondary">
              No jobs yet. Run an extension to see job history.
            </Typography>
          </Box>
        ) : (
          <Stack spacing={1}>
            {extensionJobs.map((job) => (
              <JobListItem key={job.id} job={job} showExtensionName={false} />
            ))}
          </Stack>
        )}
      </Collapse>
    </Box>
  );
};
