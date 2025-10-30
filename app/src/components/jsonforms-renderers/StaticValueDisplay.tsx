import { Box, TextField, Typography, Chip, Tooltip } from "@mui/material";
import TableViewIcon from "@mui/icons-material/TableView";
import ArrayFieldToolbar from "./ArrayFieldToolbar";
import { ArrayFieldType, getArrayShapeLabel, getArrayPreview } from "../../utils/arrayEditor";

interface StaticValueDisplayProps {
  value: any;
  label: string;
  required?: boolean;
  errors?: string;
  onEdit?: () => void;
  onClear?: () => void;
  fieldType?: ArrayFieldType;
  onChange?: (newValue: any) => void;
}

/**
 * ArrayShapeChip component - displays array shape as an interactive chip.
 */
function ArrayShapeChip({
  value,
  fieldType,
  label,
  onEdit,
  onDelete,
}: {
  value: (string | number)[] | (string | number)[][];
  fieldType: ArrayFieldType;
  label: string;
  onEdit: () => void;
  onDelete?: () => void;
}) {
  const shapeLabel = getArrayShapeLabel(value, fieldType);
  const preview = getArrayPreview(value, 2);

  return (
    <Tooltip title={preview} arrow>
      <Chip
        icon={<TableViewIcon />}
        label={shapeLabel}
        onClick={onEdit}
        onDelete={onDelete}
        color="primary"
        variant="outlined"
        sx={{
          cursor: 'pointer',
          fontSize: '0.875rem',
          height: '32px',
          '&:hover': {
            backgroundColor: 'primary.light',
            boxShadow: 1,
          },
        }}
      />
    </Tooltip>
  );
}

/**
 * StaticValueDisplay component for displaying arrays and numbers.
 *
 * Different rendering based on value type:
 * - Single number: editable TextField
 * - Single color: inline color picker
 * - Array: clickable chip with shape info
 *
 * Follows SRP by only handling static value display and editing.
 */
export default function StaticValueDisplay({
  value,
  label,
  required,
  errors,
  onEdit,
  onClear,
  fieldType,
  onChange,
}: StaticValueDisplayProps) {
  // Case 1: Single number - show editable TextField
  if (typeof value === 'number') {
    return (
      <Box sx={{ display: 'flex', alignItems: 'flex-start', gap: 1, marginBottom: 2 }}>
        <TextField
          fullWidth
          type="number"
          label={label}
          value={value}
          onChange={(e) => onChange?.(parseFloat(e.target.value) || 0)}
          required={required}
          error={!!errors}
          helperText={errors}
        />
        <ArrayFieldToolbar
          valueType="number"
          onEditArray={onEdit}
          onClear={onClear}
        />
      </Box>
    );
  }

  // Case 2: Single-element hex color array - show inline color picker
  const isSingleColorArray =
    fieldType === 'color' &&
    Array.isArray(value) &&
    value.length === 1 &&
    typeof value[0] === 'string' &&
    /^#[0-9A-Fa-f]{6}$/.test(value[0]);

  if (isSingleColorArray) {
    const colorValue = value[0] as string;
    return (
      <Box sx={{ display: 'flex', alignItems: 'flex-start', gap: 1, marginBottom: 2 }}>
        <Box sx={{ flex: 1 }}>
          <Typography variant="caption" sx={{ display: 'block', mb: 0.5, color: 'text.secondary' }}>
            {label}{required && ' *'}
          </Typography>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
            <input
              type="color"
              value={colorValue}
              onChange={(e) => onChange?.([e.target.value])}
              style={{
                width: '40px',
                height: '40px',
                border: '1px solid rgba(0, 0, 0, 0.23)',
                borderRadius: '4px',
                cursor: 'pointer',
              }}
              title="Pick color"
            />
            <Typography variant="body2" sx={{ fontFamily: 'monospace' }}>
              {colorValue}
            </Typography>
          </Box>
          {errors && (
            <Typography variant="caption" color="error" sx={{ display: 'block', mt: 0.5 }}>
              {errors}
            </Typography>
          )}
        </Box>
        <ArrayFieldToolbar
          valueType="array"
          onEditArray={onEdit}
          onClear={onClear}
        />
      </Box>
    );
  }

  // Case 3: Array (number[] or number[][]) - show shape chip
  if (Array.isArray(value) && fieldType) {
    return (
      <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1, marginBottom: 2 }}>
        <Typography variant="caption" sx={{ color: 'text.secondary' }}>
          {label}{required && ' *'}
        </Typography>
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <ArrayShapeChip
            value={value}
            fieldType={fieldType}
            label={label}
            onEdit={onEdit!}
            onDelete={onClear}
          />
          {errors && (
            <Typography variant="caption" color="error">
              {errors}
            </Typography>
          )}
        </Box>
      </Box>
    );
  }

  // Fallback: shouldn't reach here
  return null;
}
