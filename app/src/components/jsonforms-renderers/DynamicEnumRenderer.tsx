import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, and, uiTypeIs, ControlProps } from "@jsonforms/core";
import { Box, Autocomplete, TextField, Typography } from "@mui/material";

/**
 * StaticValueDisplay component for displaying arrays and numbers
 * Uses a disabled TextField to ensure JSONForms recognizes this as a form control
 * and preserves the data in the form submission
 */
const StaticValueDisplay = ({ value, label, required, errors }: { 
  value: any; 
  label: string;
  required?: boolean;
  errors?: string;
}) => {
  const displayValue = Array.isArray(value) 
    ? JSON.stringify(value)
    : String(value);
  
  return (
    <TextField
      fullWidth
      label={label}
      value={displayValue}
      disabled
      required={required}
      error={!!errors}
      helperText={errors || "Static value (not editable)"}
      sx={{ 
        marginBottom: 2,
        "& .MuiInputBase-input.Mui-disabled": {
          WebkitTextFillColor: "rgba(0, 0, 0, 0.6)",
          fontFamily: "monospace",
        }
      }}
    />
  );
};

/**
 * DynamicEnumRenderer - A composable renderer for dynamic enum fields
 * 
 * Features controlled by x-features array:
 * - "dynamic-atom-props": Populate dropdown from metadata keys
 * - "free-solo": Allow custom text input
 * - "color-picker": Add color picker UI alongside autocomplete
 */
const DynamicEnumRenderer = ({
  data,
  handleChange,
  path,
  label,
  schema,
  required,
  errors,
  visible,
}: ControlProps) => {
  // Extract features from schema (with type assertion for custom properties)
  const features = (schema as any)["x-features"] || [];
  const hasColorPicker = features.includes("color-picker");
  const hasFreeSolo = features.includes("free-solo");
  const hasDynamicProps = features.includes("dynamic-atom-props");

  // Get options from injected enum or empty array
  const options = schema.enum || [];
  
  // Detect if current value is static (array/number)
  const isStaticValue = Array.isArray(data) || typeof data === "number";
  
  // If value is static, show read-only disabled TextField
  // Using TextField instead of custom display ensures JSONForms recognizes this
  // as a valid form control and preserves the data in form submissions
  if (isStaticValue) {
    return <StaticValueDisplay value={data} label={label || ""} required={required} errors={errors} />;
  }
  
  // If field is hidden, don't render but preserve data
  if (!visible) {
    return null;
  }

  // Get current value with fallback to default
  const value = data ?? schema.default ?? "";

  // Determine color picker value (if hex color)
  const isHexColor =
    typeof value === "string" && value.startsWith("#") && value.length === 7;
  const colorValue = isHexColor ? value : "#000000";

  return (
    <Box sx={{ marginBottom: 2 }}>
      <Box sx={{ display: "flex", gap: 1, alignItems: "flex-start" }}>
        <Autocomplete
          freeSolo={hasFreeSolo}
          options={options}
          value={value}
          inputValue={typeof value === "string" ? value : ""}
          onChange={(_, newValue) => {
            // newValue can be null, a string from options, or typed value
            handleChange(path, newValue || "");
          }}
          onInputChange={(_, newInputValue, reason) => {
            // Only update on input change when freeSolo is enabled
            if (reason === "input" && hasFreeSolo) {
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
        
        {hasColorPicker && (
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
        )}
      </Box>
    </Box>
  );
};

/**
 * Tester function for DynamicEnumRenderer
 * Matches fields with x-custom-type="dynamic-enum"
 * Priority 10 to override default renderers
 * 
 * Note: This renderer handles both string values (editable) and static values
 * (arrays/numbers - read-only display). For static values, it renders a 
 * display-only component but doesn't modify the data, ensuring it's preserved.
 */
export const dynamicEnumTester = rankWith(
  10, // High priority
  and(
    schemaMatches((schema) => (schema as any)["x-custom-type"] === "dynamic-enum"),
    uiTypeIs("Control")
  )
);

export default withJsonFormsControlProps(DynamicEnumRenderer);
