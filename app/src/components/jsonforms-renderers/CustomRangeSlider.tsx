import { withJsonFormsControlProps } from "@jsonforms/react";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import {
  Box,
  FormLabel,
  IconButton,
  Slider,
  Tooltip,
  Typography,
} from "@mui/material";
import { ControlProps, rankWith, schemaMatches } from "@jsonforms/core";

// The component itself, now typed with ControlProps
const CustomRangeSlider = ({
  data,
  handleChange,
  path,
  label,
  schema,
  required,
}: ControlProps) => {
  // Your logic for defaults and values is solid.
  // We can add a Number() cast for extra safety in case the data is a string.
  const value = Number(data ?? schema.default ?? schema.minimum ?? 0);

  // The step property is correctly read from the schema.
  const step = schema.step || 0.1;

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
      <Box sx={{ display: "flex", alignItems: "center", gap: 2, paddingX: 1 }}>
        <Slider
          value={value}
          // The MUI Slider's onChange provides (event, newValue).
          // This correctly passes only the newValue to jsonforms.
          onChange={(_event, newValue) => handleChange(path, newValue)}
          min={schema.minimum}
          max={schema.maximum}
          step={step}
          aria-labelledby="input-slider"
        />
        <Typography sx={{ minWidth: 35 }}>{value}</Typography>
      </Box>
    </Box>
  );
};

// The tester logic is correct and now co-located with its component.
export const customRangeSliderTester = rankWith(
  5, // High rank to ensure it's chosen for schemas with format: 'range'
  schemaMatches((schema) => {
    // This will match schemas like: { "type": "number", "format": "range" }
    return (
      schema?.format === "range" &&
      (schema.type === "number" || schema.type === "integer")
    );
  }),
);

// The default export is the component wrapped in the HOC.
export default withJsonFormsControlProps(CustomRangeSlider);
