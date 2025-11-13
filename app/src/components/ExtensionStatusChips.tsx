import { Chip, Stack } from "@mui/material";
import WorkIcon from "@mui/icons-material/Work";
import QueueIcon from "@mui/icons-material/Queue";
import CloudIcon from "@mui/icons-material/Cloud";
import PlayArrowIcon from "@mui/icons-material/PlayArrow";
import CheckCircleIcon from "@mui/icons-material/CheckCircle";
import { type ExtensionMetadata } from "../hooks/useSchemas";
import { type ExtensionStats } from "../types/jobs";

interface ExtensionStatusChipsProps {
  metadata: ExtensionMetadata;
  stats?: ExtensionStats;
}

/**
 * Displays status chips for an extension showing:
 * - Provider type (server-side vs client workers)
 * - Worker availability (idle/busy)
 * - Number of pending tasks
 */
export const ExtensionStatusChips = ({
  metadata,
  stats,
}: ExtensionStatusChipsProps) => {
  return (
    <Stack direction="row" spacing={1} sx={{ mb: 2 }} flexWrap="wrap" gap={1}>
      {metadata.provider === "celery" ? (
        <Chip
          icon={<CloudIcon />}
          label="Server-side"
          color="primary"
          size="small"
        />
      ) : (
        <Chip
          icon={<WorkIcon />}
          label={`${metadata.provider} worker${metadata.provider !== 1 ? "s" : ""}`}
          color={metadata.provider > 0 ? "success" : "default"}
          size="small"
        />
      )}

      {/* Show worker availability from stats */}
      {stats && stats.idleWorkers > 0 && (
        <Chip
          icon={<CheckCircleIcon />}
          label={`${stats.idleWorkers} idle`}
          color="success"
          size="small"
        />
      )}

      {stats && stats.busyWorkers > 0 && (
        <Chip
          icon={<PlayArrowIcon />}
          label={`${stats.busyWorkers} busy`}
          color="info"
          size="small"
        />
      )}

      {stats && stats.pendingJobs > 0 && (
        <Chip
          icon={<QueueIcon />}
          label={`${stats.pendingJobs} pending`}
          color="warning"
          size="small"
        />
      )}
    </Stack>
  );
};
