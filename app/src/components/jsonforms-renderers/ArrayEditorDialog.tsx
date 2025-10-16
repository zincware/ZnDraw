import { useState, useMemo, useCallback, useEffect } from 'react';
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  Box,
  Autocomplete,
  TextField,
  Alert,
  Tooltip,
  Typography,
  Chip,
  CircularProgress,
} from '@mui/material';
import {
  DataGrid,
  GridColDef,
  GridRowsProp,
  GridRowModel,
  GridRowSelectionModel,
} from '@mui/x-data-grid';
import AddIcon from '@mui/icons-material/Add';
import DeleteIcon from '@mui/icons-material/Delete';
import ContentPasteIcon from '@mui/icons-material/ContentPaste';
import {
  ArrayFieldType,
  getFieldTypeConfig,
  normalizeToArray,
  denormalizeFromArray,
  validateArrayData,
  createDefaultRow,
  parseClipboardData,
  getDefaultArrayValue,
} from '../../utils/arrayEditor';

interface ArrayEditorDialogProps {
  open: boolean;
  value: string | number | (string | number)[] | (string | number)[][];
  fieldType: ArrayFieldType;
  fieldLabel: string;
  enumOptions: string[];
  onSave: (newValue: string | number | (string | number)[] | (string | number)[][]) => void;
  onCancel: () => void;
  /** Optional callback to fetch data from server (for string values) */
  onFetchFromServer?: (key: string) => Promise<(string | number)[] | (string | number)[][]>;
  /** Optional instance count (e.g., from position field - ground truth for particle count) */
  instanceCount?: number;
}

interface RowData {
  id: number;
  [key: string]: string | number;
}

export default function ArrayEditorDialog({
  open,
  value,
  fieldType,
  fieldLabel,
  enumOptions,
  onSave,
  onCancel,
  onFetchFromServer,
  instanceCount,
}: ArrayEditorDialogProps) {
  const config = getFieldTypeConfig(fieldType);

  // Store the original string value if applicable (updated when dialog opens)
  const [originalStringValue, setOriginalStringValue] = useState<string | null>(() =>
    typeof value === 'string' ? value : null
  );

  const [arrayData, setArrayData] = useState<(string | number)[][]>(() => {
    // If value is string, use default values
    if (typeof value === 'string') {
      return getDefaultArrayValue(fieldType);
    }
    return normalizeToArray(value, fieldType);
  });

  const [selectedRows, setSelectedRows] = useState<GridRowSelectionModel>({
    type: 'include',
    ids: new Set(),
  });
  const [showDropdownSwitch, setShowDropdownSwitch] = useState(false);
  const [selectedDropdownValue, setSelectedDropdownValue] = useState<string>('');
  const [validationErrors, setValidationErrors] = useState<string[]>([]);
  const [isLoadingFromServer, setIsLoadingFromServer] = useState(false);

  // Reset state when dialog opens with new value
  useEffect(() => {
    if (open) {
      // Update originalStringValue when dialog opens
      setOriginalStringValue(typeof value === 'string' ? value : null);

      // If value is string, use defaults; otherwise normalize
      const initialData = typeof value === 'string'
        ? getDefaultArrayValue(fieldType)
        : normalizeToArray(value, fieldType);

      setArrayData(initialData);
      setSelectedRows({ type: 'include', ids: new Set() });
      setShowDropdownSwitch(false);
      setSelectedDropdownValue('');
      setValidationErrors([]);
      setIsLoadingFromServer(false);
    }
  }, [open, value, fieldType]);

  // Check if we're in single-value mode (only one row for single-value supporting fields)
  const isSingleValueMode = config.supportsSingleValue && arrayData.length === 1;

  // Validate on data change
  useEffect(() => {
    const validation = validateArrayData(arrayData, fieldType);
    setValidationErrors(validation.errors);
  }, [arrayData, fieldType]);

  // Convert array data to DataGrid rows
  const rows: GridRowsProp<RowData> = useMemo(() => {
    return arrayData.map((row, idx) => {
      const rowData: RowData = { id: idx };
      row.forEach((val, colIdx) => {
        rowData[`col${colIdx}`] = val;
      });
      return rowData;
    });
  }, [arrayData]);

  // Define columns
  const columns: GridColDef<RowData>[] = useMemo(() => {
    const cols: GridColDef<RowData>[] = [
      {
        field: 'id',
        headerName: '#',
        width: 60,
        editable: false,
        valueGetter: (value, row) => row.id + 1,
      },
    ];

    for (let i = 0; i < config.dimensions; i++) {
      const column: GridColDef<RowData> = {
        field: `col${i}`,
        headerName: config.columnLabels[i] || `Col ${i + 1}`,
        width: config.isStringType && fieldType === 'color' ? 200 : 120,
        editable: true,
        type: config.isStringType ? 'string' : 'number',
        valueGetter: (value, row) => row[`col${i}`],
      };

      // Add custom color picker cell renderer for color fields
      if (config.isStringType && fieldType === 'color') {
        column.renderCell = (params) => {
          const cellValue = params.value as string;
          return (
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, width: '100%' }}>
              <input
                type="color"
                value={cellValue || '#000000'}
                onChange={(e) => {
                  const newValue = e.target.value;
                  const rowIndex = params.row.id;
                  const newArrayData = [...arrayData];
                  newArrayData[rowIndex] = [newValue];
                  setArrayData(newArrayData);
                }}
                style={{
                  width: '32px',
                  height: '32px',
                  border: '1px solid rgba(0, 0, 0, 0.23)',
                  borderRadius: '4px',
                  cursor: 'pointer',
                }}
              />
              <Typography variant="body2" sx={{ fontFamily: 'monospace' }}>
                {cellValue}
              </Typography>
            </Box>
          );
        };
        column.editable = false; // Disable text editing, use color picker only
      }

      cols.push(column);
    }

    return cols;
  }, [config, fieldType, arrayData]);

  // Handle cell edit
  const processRowUpdate = useCallback(
    (newRow: GridRowModel<RowData>) => {
      const rowIndex = newRow.id;
      const newArrayData = [...arrayData];
      const updatedRow: (string | number)[] = [];

      for (let i = 0; i < config.dimensions; i++) {
        const value = newRow[`col${i}`];
        if (config.isStringType) {
          updatedRow.push(typeof value === 'string' ? value : (config.defaultValue as string));
        } else {
          updatedRow.push(typeof value === 'number' ? value : (config.defaultValue as number));
        }
      }

      newArrayData[rowIndex] = updatedRow;
      setArrayData(newArrayData);

      return newRow;
    },
    [arrayData, config]
  );

  // Add new row
  const handleAddRow = useCallback(() => {
    const newRow = createDefaultRow(fieldType, arrayData);
    setArrayData([...arrayData, newRow]);
  }, [arrayData, fieldType]);

  // Delete selected rows
  const handleDeleteSelected = useCallback(() => {
    if (selectedRows.ids.size === 0) return;

    // Filter out selected rows
    const newArrayData = arrayData.filter((_, idx) => !selectedRows.ids.has(idx));

    // If we delete all rows, add one default row
    if (newArrayData.length === 0) {
      newArrayData.push(createDefaultRow(fieldType));
    }

    setArrayData(newArrayData);
    setSelectedRows({ type: 'include', ids: new Set() });
  }, [arrayData, selectedRows, fieldType]);

  // Toggle single value mode
  const handleToggleSingleValue = useCallback(() => {
    if (isSingleValueMode) {
      // Switching from single to per-instance: this should not happen for position/direction
      // For other fields, duplicate the row based on a reasonable default (or stay as-is)
      if (fieldType === 'position' || fieldType === 'direction') {
        // These fields should never have a toggle button, but guard just in case
        return;
      }
      // For other per-instance fields, expand to match instance count or use sensible default
      // Use instanceCount if available (ground truth from position field), otherwise expand by 1
      const targetCount = instanceCount ?? Math.max(2, arrayData.length + 1);
      const expandedRows = Array(targetCount).fill(arrayData[0]).map(row => [...row]);
      setArrayData(expandedRows);
    } else {
      // Switching to single value: keep only first row
      setArrayData([arrayData[0]]);
    }
  }, [arrayData, isSingleValueMode, fieldType, instanceCount]);

  // Handle paste
  const handlePaste = useCallback(async () => {
    try {
      const text = await navigator.clipboard.readText();
      const parsed = parseClipboardData(text);

      if (!parsed) {
        setValidationErrors(['Could not parse clipboard data. Expected CSV, TSV, or JSON format.']);
        return;
      }

      // Validate dimensions
      const invalidRows = parsed.filter(row => row.length !== config.dimensions);
      if (invalidRows.length > 0) {
        setValidationErrors([
          `Pasted data has invalid dimensions. Expected ${config.dimensions} columns per row.`,
        ]);
        return;
      }

      setArrayData(parsed);
      setValidationErrors([]);
    } catch (error) {
      setValidationErrors(['Failed to read clipboard. Please check browser permissions.']);
    }
  }, [config.dimensions]);

  // Handle save
  const handleSave = useCallback(() => {
    // Validate before saving
    const validation = validateArrayData(arrayData, fieldType);
    if (!validation.valid) {
      setValidationErrors(validation.errors);
      return;
    }

    const denormalized = denormalizeFromArray(arrayData, fieldType);
    onSave(denormalized);
  }, [arrayData, fieldType, onSave]);

  // Handle switch to dropdown
  const handleSwitchToDropdown = useCallback(() => {
    if (!selectedDropdownValue) return;
    onSave(selectedDropdownValue);
  }, [selectedDropdownValue, onSave]);

  // Handle load from server
  const handleLoadFromServer = useCallback(async () => {
    if (!originalStringValue || !onFetchFromServer) return;

    setIsLoadingFromServer(true);
    setValidationErrors([]);

    try {
      const fetchedData = await onFetchFromServer(originalStringValue);
      const normalizedData = normalizeToArray(fetchedData, fieldType);
      setArrayData(normalizedData);
    } catch (error) {
      setValidationErrors([
        `Failed to load data from server: ${error instanceof Error ? error.message : 'Unknown error'}`,
      ]);
    } finally {
      setIsLoadingFromServer(false);
    }
  }, [originalStringValue, onFetchFromServer, fieldType]);

  return (
    <Dialog
      open={open}
      onClose={onCancel}
      maxWidth="md"
      fullWidth
      PaperProps={{ sx: { height: '80vh' } }}
    >
      <DialogTitle>
        <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', gap: 2 }}>
          <Typography variant="h6">Edit {fieldLabel}</Typography>
          <Box sx={{ display: 'flex', gap: 1 }}>
            {/* Show Load from Server button if we have a fetchable string value */}
            {originalStringValue && onFetchFromServer && !showDropdownSwitch && (
              <Button
                size="small"
                variant="outlined"
                onClick={handleLoadFromServer}
                disabled={isLoadingFromServer}
                startIcon={isLoadingFromServer ? <CircularProgress size={16} /> : undefined}
              >
                {isLoadingFromServer ? 'Loading...' : `Load "${originalStringValue}"`}
              </Button>
            )}
            <Button
              size="small"
              onClick={() => setShowDropdownSwitch(!showDropdownSwitch)}
            >
              {showDropdownSwitch ? 'Edit Table' : 'Switch to Dropdown'}
            </Button>
          </Box>
        </Box>
      </DialogTitle>

      <DialogContent dividers>
        {validationErrors.length > 0 && (
          <Alert severity="error" sx={{ mb: 2 }} onClose={() => setValidationErrors([])}>
            {validationErrors.map((err, idx) => (
              <div key={idx}>{err}</div>
            ))}
          </Alert>
        )}

        {showDropdownSwitch ? (
          <Box sx={{ p: 2 }}>
            <Typography variant="body2" sx={{ mb: 2 }}>
              Select a dynamic value from the dropdown to replace the current array data:
            </Typography>
            <Autocomplete
              options={enumOptions}
              value={selectedDropdownValue}
              onChange={(_, newValue) => setSelectedDropdownValue(newValue || '')}
              renderInput={(params) => (
                <TextField
                  {...params}
                  label="Select Value"
                  placeholder="e.g., arrays.positions"
                />
              )}
              freeSolo
              fullWidth
            />
            {config.supportsSingleValue && (
              <Alert severity="info" sx={{ mt: 2 }}>
                Note: You can also enter a single numeric value instead of selecting from the dropdown.
              </Alert>
            )}
          </Box>
        ) : (
          <Box sx={{ display: 'flex', flexDirection: 'column', gap: 2, height: 'calc(80vh - 200px)' }}>
            <Box sx={{ flex: 1, minHeight: 0 }}>
              <DataGrid
                key={arrayData.length} // Force re-render when row count changes
                rows={rows}
                columns={columns}
                processRowUpdate={processRowUpdate}
                onProcessRowUpdateError={(error) => {
                  console.error('Row update error:', error);
                  setValidationErrors(['Failed to update row. Please check your input.']);
                }}
                checkboxSelection
                rowSelectionModel={selectedRows}
                onRowSelectionModelChange={(newSelection) => {
                  setSelectedRows(newSelection);
                }}
                sx={{
                  '& .MuiDataGrid-cell:focus': {
                    outline: 'none',
                  },
                  '& .MuiDataGrid-cell:focus-within': {
                    outline: '2px solid #1976d2',
                  },
                }}
              />
            </Box>

            {/* Action buttons below the table */}
            <Box sx={{ display: 'flex', gap: 1, alignItems: 'center', flexWrap: 'wrap' }}>
              <Button
                size="small"
                startIcon={<AddIcon />}
                onClick={handleAddRow}
                variant="outlined"
              >
                Add Row
              </Button>
              {selectedRows.ids.size > 0 && (
                <Button
                  size="small"
                  startIcon={<DeleteIcon />}
                  onClick={handleDeleteSelected}
                  color="error"
                  variant="outlined"
                >
                  Delete Selected ({selectedRows.ids.size})
                </Button>
              )}
              <Tooltip title="Paste data from clipboard (CSV, TSV, or JSON)">
                <Button
                  size="small"
                  startIcon={<ContentPasteIcon />}
                  onClick={handlePaste}
                  variant="outlined"
                >
                  Paste
                </Button>
              </Tooltip>
              {config.supportsSingleValue && fieldType !== 'position' && fieldType !== 'direction' && (
                <Box sx={{ ml: 'auto' }}>
                  <Chip
                    label={isSingleValueMode ? "Single Value (applies to all)" : "Per Instance"}
                    color={isSingleValueMode ? "primary" : "default"}
                    size="small"
                    onClick={handleToggleSingleValue}
                    sx={{ cursor: 'pointer' }}
                  />
                </Box>
              )}
            </Box>
          </Box>
        )}
      </DialogContent>

      <DialogActions sx={{ p: 2 }}>
        <Button onClick={onCancel}>Cancel</Button>
        {showDropdownSwitch ? (
          <Button
            onClick={handleSwitchToDropdown}
            variant="contained"
            disabled={!selectedDropdownValue}
          >
            Switch to Dropdown
          </Button>
        ) : (
          <Button
            onClick={handleSave}
            variant="contained"
            disabled={validationErrors.length > 0}
          >
            Save ({arrayData.length} {arrayData.length === 1 ? 'row' : 'rows'})
          </Button>
        )}
      </DialogActions>
    </Dialog>
  );
}
