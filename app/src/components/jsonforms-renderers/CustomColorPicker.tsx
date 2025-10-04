// src/components/jsonforms-renderers/CustomColorPicker.tsx

import { withJsonFormsControlProps } from "@jsonforms/react";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import { Box, FormLabel, IconButton, Tooltip } from "@mui/material";
import { rankWith, schemaMatches, ControlProps } from "@jsonforms/core";

// Define the component as a pure UI component.
// It receives all its props from the HOC.
const CustomColorPicker = ({
  data,
  handleChange,
  path,
  label,
  schema,
  required,
}: ControlProps) => {
  // The initial value logic is correct.
  const value = data ?? schema.default ?? "#000000";

  return (
    <Box sx={{ marginBottom: 2 }}>
      <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
        <FormLabel required={required}>{label}</FormLabel>
        {schema.description && (
          <Tooltip title={schema.description} placement="top">
            <IconButton size="small" sx={{ padding: 0.5 }}>
              <HelpOutlineIcon fontSize="small" />
            </IconButton>
          </Tooltip>
        )}
      </Box>
      <Box sx={{ display: "flex", alignItems: "center", gap: 2, marginTop: 1 }}>
        <input
          type="color"
          value={value}
          // This will now work because handleChange is correctly wired by the HOC.
          onChange={(e) => handleChange(path, e.target.value)}
          style={{
            width: "40px",
            height: "40px",
            border: "none",
            background: "none",
            cursor: "pointer",
            // A small visual improvement to show a border around the color swatch
            borderRadius: "4px",
          }}
        />
        <code>{value}</code>
      </Box>
    </Box>
  );
};

// The tester remains the same. It's correct.
export const customColorPickerTester = rankWith(
  5, // High rank ensures it's chosen over default string renderers
  schemaMatches((schema) => {
    return schema?.format === "color" && schema.type === "string";
  }),
);

// Export the component wrapped in the HOC. This is the crucial step.
export default withJsonFormsControlProps(CustomColorPicker);
