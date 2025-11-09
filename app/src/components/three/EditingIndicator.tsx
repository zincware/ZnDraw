import { Box, Chip, Typography } from "@mui/material";
import EditIcon from "@mui/icons-material/Edit";
import { useAppStore } from "../../store";

/**
 * EditingIndicator shows when the user is in geometry editing mode
 * Similar to DrawingIndicator but for transform controls editing
 */
export default function EditingIndicator() {
  // Use individual selectors to prevent unnecessary re-renders
  const mode = useAppStore((state) => state.mode);

  // Only show if editing mode is active
  if (mode !== 'editing') {
    return null;
  }

  return (
    <Box
      sx={{
        position: "fixed",
        bottom: 80,
        right: 16,
        zIndex: (theme) => theme.zIndex.tooltip,
      }}
    >
      <Chip
        icon={<EditIcon />}
        label={
          <Typography variant="body2">Editing Mode</Typography>
        }
        color="secondary"
        size="medium"
        sx={{
          "& .MuiChip-label": {
            px: 1,
          },
        }}
      />
    </Box>
  );
}
