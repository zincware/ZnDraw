import { Box, IconButton, Tooltip } from '@mui/material';
import TableChartIcon from '@mui/icons-material/TableChart';
import ClearIcon from '@mui/icons-material/Clear';

interface ArrayFieldToolbarProps {
  /** Current value type */
  valueType: 'string' | 'number' | 'array';
  /** Whether the field supports array editing */
  supportsArrayEdit?: boolean;
  /** Callback when user wants to edit the array in table view */
  onEditArray?: () => void;
  /** Callback when user wants to clear the value */
  onClear?: () => void;
}

/**
 * Toolbar with action buttons for array field editing
 */
export default function ArrayFieldToolbar({
  valueType,
  supportsArrayEdit = true,
  onEditArray,
  onClear,
}: ArrayFieldToolbarProps) {
  return (
    <Box
      sx={{
        display: 'flex',
        gap: 0.5,
        alignItems: 'center',
        ml: 1,
      }}
    >
      {/* Edit Array Button - shown for all value types */}
      {supportsArrayEdit && onEditArray && (
        <Tooltip title={valueType === 'string' ? 'Edit as table' : 'Edit in table'}>
          <IconButton size="small" onClick={onEditArray} color="primary">
            <TableChartIcon fontSize="small" />
          </IconButton>
        </Tooltip>
      )}

      {/* Clear Button - shown when value is array or number */}
      {(valueType === 'array' || valueType === 'number') && onClear && (
        <Tooltip title="Clear and switch back to dropdown">
          <IconButton size="small" onClick={onClear} color="error">
            <ClearIcon fontSize="small" />
          </IconButton>
        </Tooltip>
      )}
    </Box>
  );
}
