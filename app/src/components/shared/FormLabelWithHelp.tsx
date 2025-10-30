// src/components/shared/FormLabelWithHelp.tsx
// Reusable form label with optional help tooltip

import { Box, FormLabel, IconButton, Tooltip } from "@mui/material";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";

interface FormLabelWithHelpProps {
  /** Label text */
  label: string;
  /** Whether the field is required */
  required?: boolean;
  /** Help text to display in tooltip */
  helpText?: string;
}

/**
 * Form label with an optional help icon tooltip.
 * Provides consistent styling for form field labels.
 */
export const FormLabelWithHelp: React.FC<FormLabelWithHelpProps> = ({
  label,
  required = false,
  helpText,
}) => {
  return (
    <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
      <FormLabel required={required}>{label}</FormLabel>
      {helpText && (
        <Tooltip title={helpText} placement="top">
          <IconButton size="small" sx={{ padding: 0.5 }}>
            <HelpOutlineIcon fontSize="small" />
          </IconButton>
        </Tooltip>
      )}
    </Box>
  );
};

export default FormLabelWithHelp;
