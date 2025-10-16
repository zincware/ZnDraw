import { useState, useCallback } from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, and, uiTypeIs, ControlProps } from "@jsonforms/core";
import { Box, Autocomplete, TextField, Typography, Chip, Tooltip } from "@mui/material";
import TableViewIcon from "@mui/icons-material/TableView";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import ArrayFieldToolbar from "./ArrayFieldToolbar";
import ArrayEditorDialog from "./ArrayEditorDialog";
import { inferFieldType, ArrayFieldType, getArrayShapeLabel, getArrayPreview } from "../../utils/arrayEditor";
import { useAppStore } from "../../store";
import { getFrames } from "../../myapi/client";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";

/**
 * ArrayShapeChip component - displays array shape as an interactive chip
 */
const ArrayShapeChip = ({
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
}) => {
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
};

/**
 * StaticValueDisplay component for displaying arrays and numbers
 * Different rendering based on value type:
 * - Single number: editable TextField
 * - Single color: inline color picker
 * - Array: clickable chip with shape info
 */
const StaticValueDisplay = ({
  value,
  label,
  required,
  errors,
  onEdit,
  onClear,
  fieldType,
  onChange,
}: {
  value: any;
  label: string;
  required?: boolean;
  errors?: string;
  onEdit?: () => void;
  onClear?: () => void;
  fieldType?: ArrayFieldType;
  onChange?: (newValue: any) => void;
}) => {
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
};

/**
 * DynamicEnumRenderer - A composable renderer for dynamic enum fields
 *
 * Features controlled by x-features array:
 * - "dynamic-atom-props": Populate dropdown from metadata keys
 * - "free-solo": Allow custom text input
 * - "color-picker": Add color picker UI alongside autocomplete
 * - "editable-array": Enable array editing via table (automatically enabled for arrays)
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
  rootData,
}: ControlProps & { rootData?: any }) => {
  // Extract features from schema (with type assertion for custom properties)
  const features = (schema as any)["x-features"] || [];
  const hasColorPicker = features.includes("color-picker");
  const hasFreeSolo = features.includes("free-solo");
  const hasDynamicProps = features.includes("dynamic-atom-props");

  // Get options from injected enum or empty array
  const options = schema.enum || [];

  // Detect if current value is static (array/number)
  const isStaticValue = Array.isArray(data) || typeof data === "number";

  // State for array editor dialog
  const [arrayEditorOpen, setArrayEditorOpen] = useState(false);

  // Get app store for fetching frame data
  const { roomId, currentFrame, clientId } = useAppStore();

  // Infer field type from path and schema
  const fieldType: ArrayFieldType = inferFieldType(path, schema);

  // Compute instance count from position data (ground truth for instance count)
  const instanceCount = (() => {
    const position = rootData?.position;
    if (position) {
      // If position is a 2D array of tuples or numbers
      if (Array.isArray(position) && Array.isArray(position[0])) {
        return position.length; // Number of instances
      }
      // If position is a single tuple/array
      if (Array.isArray(position) && position.length > 0 && !Array.isArray(position[0])) {
        return 1; // Single position = 1 instance
      }
    }
    return undefined;
  })();

  // Determine if the current value is a fetchable string (not a hex color)
  const isFetchableString = typeof data === 'string' && shouldFetchAsFrameData(data);

  // React-query for fetching server data (only enabled when needed)
  const { data: fetchedServerData, refetch: refetchServerData, isFetching: isFetchingServerData } = useQuery({
    queryKey: ["frame", roomId, currentFrame, data],
    queryFn: ({ signal }: { signal: AbortSignal }) =>
      getFrames(roomId!, currentFrame, [data as string], signal),
    enabled: false, // Manual fetching only
    placeholderData: keepPreviousData,
  });

  // Handle opening array editor
  const handleOpenArrayEditor = useCallback(() => {
    setArrayEditorOpen(true);
  }, []);

  // Handle saving from array editor
  const handleSaveArrayEditor = useCallback((newValue: string | number | number[] | number[][]) => {
    handleChange(path, newValue);
    setArrayEditorOpen(false);
  }, [handleChange, path]);

  // Handle clearing static value back to empty string (dropdown mode)
  const handleClearValue = useCallback(() => {
    handleChange(path, "");
  }, [handleChange, path]);

  // Callback for loading from server (triggers refetch)
  const handleLoadFromServer = useCallback(async (): Promise<number[] | number[][]> => {
    const result = await refetchServerData();
    if (result.data && typeof data === 'string') {
      const fetchedData = result.data[data as string];
      if (fetchedData) {
        // Convert TypedArray to regular array
        const converted = Array.from(fetchedData as ArrayLike<number>);
        return converted as number[] | number[][];
      }
    }
    throw new Error('No data returned from server');
  }, [refetchServerData, data]);

  // If value is static, show read-only disabled TextField with edit button
  // Using TextField instead of custom display ensures JSONForms recognizes this
  // as a valid form control and preserves the data in form submissions
  if (isStaticValue) {
    return (
      <>
        <StaticValueDisplay
          value={data}
          label={label || ""}
          required={required}
          errors={errors}
          onEdit={handleOpenArrayEditor}
          onClear={handleClearValue}
          fieldType={fieldType}
          onChange={(newValue) => handleChange(path, newValue)}
        />
        <ArrayEditorDialog
          open={arrayEditorOpen}
          value={data}
          fieldType={fieldType}
          fieldLabel={label || path}
          enumOptions={options as string[]}
          onSave={handleSaveArrayEditor}
          onCancel={() => setArrayEditorOpen(false)}
          onFetchFromServer={isFetchableString ? handleLoadFromServer : undefined}
          instanceCount={instanceCount}
        />
      </>
    );
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
    <>
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

          {/* Add edit button for string values to open table with defaults */}
          <ArrayFieldToolbar
            valueType="string"
            onEditArray={handleOpenArrayEditor}
          />
        </Box>
      </Box>

      {/* Array editor dialog (available for all value types) */}
      <ArrayEditorDialog
        open={arrayEditorOpen}
        value={data}
        fieldType={fieldType}
        fieldLabel={label || path}
        enumOptions={options as string[]}
        onSave={handleSaveArrayEditor}
        onCancel={() => setArrayEditorOpen(false)}
        onFetchFromServer={isFetchableString ? handleLoadFromServer : undefined}
        instanceCount={instanceCount}
      />
    </>
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
