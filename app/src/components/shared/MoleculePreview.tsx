// src/components/shared/MoleculePreview.tsx
// Reusable component for displaying molecule structure images

import { Box, CircularProgress, Alert, Tooltip, Typography } from "@mui/material";
import { useMoleculeImage } from "../../hooks/useMoleculeImage";

interface MoleculePreviewProps {
  /** SMILES notation string */
  smiles: string;
  /** Optional click handler */
  onClick?: () => void;
  /** Size variant */
  size?: "small" | "medium" | "large";
  /** Custom alt text for the image */
  alt?: string;
  /** Show error as tooltip instead of alert */
  errorAsTooltip?: boolean;
  /** Throttle image fetching (useful for real-time editing) */
  throttleMs?: number;
}

const SIZE_CONFIG = {
  small: {
    height: 60,
    spinnerSize: 24,
  },
  medium: {
    height: 150,
    spinnerSize: 40,
  },
  large: {
    height: 300,
    spinnerSize: 60,
  },
};

/**
 * Displays a molecule structure image from SMILES notation.
 * Handles loading, error states, and optional click interactions.
 */
export const MoleculePreview: React.FC<MoleculePreviewProps> = ({
  smiles,
  onClick,
  size = "medium",
  alt = "Molecule structure",
  errorAsTooltip = false,
  throttleMs = 0,
}) => {
  const { image, loading, error } = useMoleculeImage(smiles, throttleMs);
  const config = SIZE_CONFIG[size];

  // No SMILES provided
  if (!smiles || smiles.trim() === "") {
    return null;
  }

  // Loading state
  if (loading) {
    return (
      <Box
        sx={{
          display: "flex",
          alignItems: "center",
          justifyContent: "center",
          height: config.height,
        }}
      >
        <CircularProgress size={config.spinnerSize} />
      </Box>
    );
  }

  // Error state
  if (error) {
    if (errorAsTooltip) {
      return (
        <Tooltip title={error}>
          <Typography variant="caption" color="error">
            Error
          </Typography>
        </Tooltip>
      );
    }
    return (
      <Alert severity="error" sx={{ width: "100%" }}>
        {error}
      </Alert>
    );
  }

  // Success state - show image
  if (image) {
    return (
      <Box
        component="img"
        src={image}
        alt={alt}
        onClick={onClick}
        sx={{
          width: "100%",
          height: config.height,
          objectFit: "contain",
          backgroundColor: "white",
          cursor: onClick ? "pointer" : "default",
          ...(onClick && {
            "&:hover": {
              opacity: 0.8,
            },
          }),
        }}
      />
    );
  }

  return null;
};

export default MoleculePreview;
