import { useState, useCallback } from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, and, uiTypeIs, ControlProps } from "@jsonforms/core";
import { useQuery, keepPreviousData } from "@tanstack/react-query";
import { useAppStore } from "../../store";
import { getFrames } from "../../myapi/client";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import { isTransform } from "../../utils/transformProcessor";
import { inferFieldType, ArrayFieldType } from "../../utils/arrayEditor";
import ArrayEditorDialog from "./ArrayEditorDialog";
import TransformEditor from "./TransformEditor";
import StaticValueDisplay from "./StaticValueDisplay";
import StringEnumControl from "./StringEnumControl";

/**
 * DynamicEnumRenderer - A composable renderer for dynamic enum fields
 *
 * This is now a lightweight orchestrator that delegates to focused sub-components:
 * - TransformEditor: Handles transform object editing
 * - StaticValueDisplay: Handles arrays and numbers display
 * - StringEnumControl: Handles string dropdown with optional features
 * - ArrayEditorDialog: Handles array editing in a dialog
 *
 * Features controlled by x-features array:
 * - "dynamic-atom-props": Populate dropdown from metadata keys
 * - "free-solo": Allow custom text input
 * - "color-picker": Add color picker UI alongside autocomplete
 * - "editable-array": Enable array editing via table
 * - "transform": Enable transform mode to filter data based on frame data
 *
 * Transform mode allows creating InArrayTransform objects that filter geometry
 * data (positions, radii, colors) based on indices extracted from frame data
 * (e.g., constraint indices).
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
  const hasEditableArray = features.includes("editable-array");
  const hasTransform = features.includes("transform");

  // Get options from injected enum or empty array
  const options = schema.enum || [];

  // Detect value types
  const isTransformValue = isTransform(data);
  const isStaticValue = Array.isArray(data) || typeof data === "number";

  // State for array editor dialog
  const [arrayEditorOpen, setArrayEditorOpen] = useState(false);

  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const currentFrame = useAppStore((state) => state.currentFrame);
  const userName = useAppStore((state) => state.userName);

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
  const handleSaveArrayEditor = useCallback((newValue: string | number | (string | number)[] | (string | number)[][]) => {
    handleChange(path, newValue);
    setArrayEditorOpen(false);
  }, [handleChange, path]);

  // Handle clearing static value back to empty string (dropdown mode)
  const handleClearValue = useCallback(() => {
    handleChange(path, "");
  }, [handleChange, path]);

  // Handle creating a new transform
  const handleCreateTransform = useCallback(() => {
    handleChange(path, {
      type: "in_array",
      source: "",
      path: "",
      filter: "",
    });
  }, [handleChange, path]);

  // Handle clearing transform back to empty string
  const handleClearTransform = useCallback(() => {
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

  // If field is hidden, don't render but preserve data
  if (!visible) {
    return null;
  }

  // --- Component Delegation ---

  // If value is a transform, delegate to TransformEditor
  if (isTransformValue && hasTransform) {
    return (
      <TransformEditor
        value={data}
        label={label || ""}
        required={required}
        onChange={(newValue) => handleChange(path, newValue)}
        onClear={handleClearTransform}
      />
    );
  }

  // If value is static (array/number), delegate to StaticValueDisplay
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

  // Otherwise, delegate to StringEnumControl
  const value = data ?? schema.default ?? "";

  return (
    <>
      <StringEnumControl
        value={value}
        label={label || ""}
        required={required}
        errors={errors}
        options={options as string[]}
        hasFreeSolo={hasFreeSolo}
        hasColorPicker={hasColorPicker}
        hasEditableArray={hasEditableArray}
        hasTransform={hasTransform}
        onChange={(newValue) => handleChange(path, newValue)}
        onOpenArrayEditor={handleOpenArrayEditor}
        onCreateTransform={handleCreateTransform}
      />

      {/* Array editor dialog (only available if editable-array feature is enabled) */}
      {hasEditableArray && (
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
      )}
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
