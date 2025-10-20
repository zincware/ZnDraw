import { Box, Typography, Paper } from "@mui/material";
import CloudUploadIcon from "@mui/icons-material/CloudUpload";

interface DropOverlayProps {
  isDragging: boolean;
}

export default function DropOverlay({ isDragging }: DropOverlayProps) {
  if (!isDragging) return null;

  return (
    <Box
      sx={{
        position: "absolute",
        top: 0,
        left: 0,
        right: 0,
        bottom: 0,
        bgcolor: "rgba(0, 0, 0, 0.7)",
        backdropFilter: "blur(4px)",
        zIndex: 9999,
        display: "flex",
        alignItems: "center",
        justifyContent: "center",
        pointerEvents: "none",
      }}
    >
      <Paper
        elevation={8}
        sx={{
          p: 4,
          display: "flex",
          flexDirection: "column",
          alignItems: "center",
          gap: 2,
          bgcolor: "background.paper",
          border: 2,
          borderColor: "primary.main",
          borderStyle: "dashed",
          borderRadius: 2,
        }}
      >
        <CloudUploadIcon sx={{ fontSize: 64, color: "primary.main" }} />
        <Typography variant="h4" component="div">
          Drop file to upload
        </Typography>
        <Typography variant="body2" color="text.secondary">
          Supported formats: XYZ, PDB, CIF, H5MD, and more
        </Typography>
      </Paper>
    </Box>
  );
}
