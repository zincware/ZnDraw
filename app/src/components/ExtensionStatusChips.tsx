import { Chip, Stack } from "@mui/material";
import WorkIcon from "@mui/icons-material/Work";
import QueueIcon from "@mui/icons-material/Queue";
import CloudIcon from "@mui/icons-material/Cloud";
import PlayArrowIcon from "@mui/icons-material/PlayArrow";
import { type ExtensionMetadata } from "../hooks/useSchemas";

interface ExtensionStatusChipsProps {
  metadata: ExtensionMetadata;
}

/**
 * Displays status chips for an extension showing:
 * - Provider type (server-side vs client workers)
 * - Number of running tasks
 * - Number of queued tasks
 */
export const ExtensionStatusChips = ({
  metadata,
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
      {metadata.progressingWorkers > 0 && (
        <Chip
          icon={<PlayArrowIcon />}
          label={`${metadata.progressingWorkers} running`}
          color="info"
          size="small"
        />
      )}
      {metadata.queueLength > 0 && (
        <Chip
          icon={<QueueIcon />}
          label={`${metadata.queueLength} queued`}
          color="warning"
          size="small"
        />
      )}
    </Stack>
  );
};
