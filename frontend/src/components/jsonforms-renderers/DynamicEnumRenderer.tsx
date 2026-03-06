import {
	type ControlProps,
	and,
	rankWith,
	schemaMatches,
	uiTypeIs,
} from "@jsonforms/core";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { keepPreviousData, useQuery } from "@tanstack/react-query";
import { useCallback, useEffect, useState } from "react";
import { useAvailableProperties } from "../../hooks/usePropertyInspector";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import { shouldFetchAsFrameData } from "../../utils/colorUtils";
import { isTransform } from "../../utils/transformProcessor";
import ArrayEditorDialog from "./ArrayEditorDialog";
import StaticValueDisplay from "./StaticValueDisplay";
import StringEnumControl from "./StringEnumControl";
import TransformEditor from "./TransformEditor";

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
}: ControlProps) => {
	// Extract features from schema (with type assertion for custom properties)
	const features = (schema as any)["x-features"] || [];
	const hasColorPicker = features.includes("color-picker");
	const hasFreeSolo = features.includes("free-solo");
	const hasEditableArray = features.includes("editable-array");
	const hasTransform = features.includes("transform");
	const hasDynamicAtomProps = features.includes("dynamic-atom-props");
	const hasStep = features.includes("step");

	// Detect value types
	const isTransformValue = isTransform(data);
	const isStaticValue = Array.isArray(data) || typeof data === "number";

	// State for array editor dialog
	const [arrayEditorOpen, setArrayEditorOpen] = useState(false);

	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);
	const currentFrame = useAppStore((state) => state.currentFrame);
	const particleCount = useAppStore((state) => state.particleCount);
	const showSnackbar = useAppStore((state) => state.showSnackbar);

	// Fetch available properties using the same query key as PropertyInspector.
	// Enabled only for fields with the dynamic-atom-props feature.
	const {
		data: atomProperties,
		isLoading: isLoadingAtomProps,
		isError: isAtomPropsError,
	} = useAvailableProperties(
		roomId || undefined,
		currentFrame,
		particleCount,
		hasDynamicAtomProps,
	);

	useEffect(() => {
		if (isAtomPropsError) {
			showSnackbar("Failed to load atom properties", "error");
		}
	}, [isAtomPropsError, showSnackbar]);

	// For dynamic-atom-props fields: options come from the same source as PropertyInspector.
	// For other fields: options come from the schema enum (populated by injectDynamicEnums).
	const options: string[] = hasDynamicAtomProps
		? [
				...(hasStep ? ["step"] : []),
				...(atomProperties?.perParticle.map((p) => p.key) ?? []),
				...(atomProperties?.global.map((g) => g.key) ?? []),
			]
		: ((schema.enum as string[]) ?? []);

	// Determine if the current value is a fetchable string (not a hex color)
	const isFetchableString =
		typeof data === "string" && shouldFetchAsFrameData(data);

	// React-query for fetching server data (only enabled when needed)
	const { refetch: refetchServerData } = useQuery({
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
	const handleSaveArrayEditor = useCallback(
		(
			newValue: string | number | (string | number)[] | (string | number)[][],
		) => {
			handleChange(path, newValue);
			setArrayEditorOpen(false);
		},
		[handleChange, path],
	);

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
	const handleLoadFromServer = useCallback(async (): Promise<
		number[] | number[][]
	> => {
		const result = await refetchServerData();
		if (result.data && typeof data === "string") {
			const fetchedData = result.data[data as string];
			if (fetchedData) {
				// Convert TypedArray to regular array
				const converted = Array.from(fetchedData as ArrayLike<number>);
				return converted as number[] | number[][];
			}
		}
		throw new Error("No data returned from server");
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
					schema={schema}
					onChange={(newValue) => handleChange(path, newValue)}
				/>
				<ArrayEditorDialog
					open={arrayEditorOpen}
					value={data}
					schema={schema}
					fieldLabel={label || path}
					enumOptions={options as string[]}
					onSave={handleSaveArrayEditor}
					onCancel={() => setArrayEditorOpen(false)}
					onFetchFromServer={
						isFetchableString ? handleLoadFromServer : undefined
					}
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
				isLoading={hasDynamicAtomProps && isLoadingAtomProps}
				onChange={(newValue) => handleChange(path, newValue)}
				onOpenArrayEditor={handleOpenArrayEditor}
				onCreateTransform={handleCreateTransform}
			/>

			{/* Array editor dialog (only available if editable-array feature is enabled) */}
			{hasEditableArray && (
				<ArrayEditorDialog
					open={arrayEditorOpen}
					value={data}
					schema={schema}
					fieldLabel={label || path}
					enumOptions={options as string[]}
					onSave={handleSaveArrayEditor}
					onCancel={() => setArrayEditorOpen(false)}
					onFetchFromServer={
						isFetchableString ? handleLoadFromServer : undefined
					}
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
		schemaMatches(
			(schema) => (schema as any)["x-custom-type"] === "dynamic-enum",
		),
		uiTypeIs("Control"),
	),
);

export default withJsonFormsControlProps(DynamicEnumRenderer);
