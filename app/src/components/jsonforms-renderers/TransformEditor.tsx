import { Box, TextField, Typography, Tooltip, IconButton, Paper } from "@mui/material";
import FilterAltIcon from "@mui/icons-material/FilterAlt";
import CloseIcon from "@mui/icons-material/Close";
import { Transform } from "../../utils/transformProcessor";

interface TransformEditorProps {
  value: Transform;
  label: string;
  required?: boolean;
  onChange: (newValue: Transform) => void;
  onClear: () => void;
}

/**
 * TransformEditor component - inline editor for transform objects.
 *
 * Provides a focused UI for editing InArrayTransform objects with fields for:
 * - source: Frame data key containing indices
 * - path: Dot-separated path to extract indices
 * - filter: Frame data key to filter
 *
 * Follows SRP by only handling transform object editing.
 */
export default function TransformEditor({
  value,
  label,
  required,
  onChange,
  onClear,
}: TransformEditorProps) {
  const transform = value || { type: "in_array", source: "", path: "", filter: "" };

  const updateField = (field: string, newValue: string) => {
    onChange({
      ...transform,
      [field]: newValue,
    } as Transform);
  };

  return (
    <Box sx={{ marginBottom: 2 }}>
      <Box sx={{ display: "flex", alignItems: "center", justifyContent: "space-between", marginBottom: 1 }}>
        <Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
          <FilterAltIcon fontSize="small" color="primary" />
          <Typography variant="subtitle2" color="primary">
            {label}{required && ' *'} (Transform)
          </Typography>
        </Box>
        <Tooltip title="Remove transform and use simple key">
          <IconButton size="small" onClick={onClear}>
            <CloseIcon fontSize="small" />
          </IconButton>
        </Tooltip>
      </Box>

      <Paper elevation={1} sx={{ padding: 2, backgroundColor: "#f9f9f9" }}>
        <Box sx={{ display: "flex", flexDirection: "column", gap: 1.5 }}>
          <TextField
            fullWidth
            size="small"
            label="Source"
            value={transform.source || ""}
            onChange={(e) => updateField("source", e.target.value)}
            placeholder="e.g., constraints"
            helperText="Frame data key containing indices"
          />
          <TextField
            fullWidth
            size="small"
            label="Path"
            value={transform.path || ""}
            onChange={(e) => updateField("path", e.target.value)}
            placeholder="e.g., 0.kwargs.indices"
            helperText="Dot-separated path to extract indices"
          />
          <TextField
            fullWidth
            size="small"
            label="Filter"
            value={transform.filter || ""}
            onChange={(e) => updateField("filter", e.target.value)}
            placeholder="e.g., arrays.positions"
            helperText="Frame data key to filter"
          />
        </Box>
      </Paper>
    </Box>
  );
}
