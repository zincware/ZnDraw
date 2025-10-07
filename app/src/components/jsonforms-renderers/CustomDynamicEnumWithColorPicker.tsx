import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, and, ControlProps } from "@jsonforms/core";
import { Box, Autocomplete, TextField } from "@mui/material";

const CustomDynamicEnumWithColorPicker = ({
  data,
  handleChange,
  path,
  label,
  schema,
  required,
  errors,
}: ControlProps) => {
  const enumValues = schema.enum || [];
  const value = data ?? schema.default ?? "";

  // Determine if value looks like a hex color for the color picker
  const isHexColor =
    typeof value === "string" && value.startsWith("#") && value.length === 7;
  const colorValue = isHexColor ? value : "#000000";

  return (
    <Box sx={{ marginBottom: 2 }}>
      <Box sx={{ display: "flex", gap: 1, alignItems: "flex-start" }}>
        <Autocomplete
          freeSolo
          options={enumValues}
          value={value}
          inputValue={value}
          onChange={(_, newValue) => {
            // newValue can be null, a string from options, or typed value
            handleChange(path, newValue || "");
          }}
          onInputChange={(_, newInputValue, reason) => {
            // Only update on input change, not on blur/reset
            if (reason === "input") {
              handleChange(path, newInputValue);
            }
          }}
          renderInput={(params) => (
            <TextField
              {...params}
              label={label}
              required={required}
              error={!!errors}
              helperText={errors}
            />
          )}
          fullWidth
        />
        <input
          type="color"
          value={colorValue}
          onChange={(e) => handleChange(path, e.target.value)}
          style={{
            width: "50px",
            height: "40px",
            marginTop: "8px",
            border: "1px solid rgba(0, 0, 0, 0.23)",
            borderRadius: "4px",
            cursor: "pointer",
          }}
          title="Pick color"
        />
      </Box>
    </Box>
  );
};

// High priority tester - matches fields with both x-dynamic-enum AND x-color-picker
export const customDynamicEnumWithColorPickerTester = rankWith(
  10, // Higher priority than standard renderers
  and(
    schemaMatches((schema) => schema?.hasOwnProperty("x-dynamic-enum")),
    schemaMatches((schema) => schema?.["x-color-picker"] === true)
  )
);

export default withJsonFormsControlProps(CustomDynamicEnumWithColorPicker);
