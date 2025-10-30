import { Box, Autocomplete, TextField, Tooltip, IconButton } from "@mui/material";
import FilterAltIcon from "@mui/icons-material/FilterAlt";
import ArrayFieldToolbar from "./ArrayFieldToolbar";

interface StringEnumControlProps {
  value: string;
  label: string;
  required?: boolean;
  errors?: string;
  options: string[];
  hasFreeSolo: boolean;
  hasColorPicker: boolean;
  hasEditableArray: boolean;
  hasTransform: boolean;
  onChange: (newValue: string) => void;
  onOpenArrayEditor?: () => void;
  onCreateTransform?: () => void;
}

/**
 * StringEnumControl component for handling string dropdown with optional features.
 *
 * Features:
 * - Autocomplete dropdown with optional free-solo
 * - Optional color picker for hex colors
 * - Optional array editor button
 * - Optional transform mode button
 *
 * Follows SRP by only handling string value selection and UI controls.
 */
export default function StringEnumControl({
  value,
  label,
  required,
  errors,
  options,
  hasFreeSolo,
  hasColorPicker,
  hasEditableArray,
  hasTransform,
  onChange,
  onOpenArrayEditor,
  onCreateTransform,
}: StringEnumControlProps) {
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
          getOptionLabel={(option) => {
            // Handle non-string options gracefully
            if (typeof option === "string") return option;
            if (Array.isArray(option)) return JSON.stringify(option);
            return String(option);
          }}
          onChange={(_, newValue) => {
            // newValue can be null, a string from options, or typed value
            onChange(newValue || "");
          }}
          onInputChange={(_, newInputValue, reason) => {
            // Only update on input change when freeSolo is enabled
            if (reason === "input" && hasFreeSolo) {
              onChange(newInputValue);
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
            onChange={(e) => onChange(e.target.value)}
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

        {hasEditableArray && (
          <Box sx={{ marginTop: "8px" }}>
            <ArrayFieldToolbar
              valueType="string"
              onEditArray={onOpenArrayEditor}
            />
          </Box>
        )}

        {hasTransform && (
          <Tooltip title="Enable transform mode to filter data based on frame data">
            <IconButton
              size="small"
              onClick={onCreateTransform}
              sx={{ marginTop: "8px" }}
              color="primary"
            >
              <FilterAltIcon />
            </IconButton>
          </Tooltip>
        )}
      </Box>
    </Box>
  );
}
