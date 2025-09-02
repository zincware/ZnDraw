import type { ErrorObject, JsonSchema, UISchemaElement } from "@jsonforms/core";
import {
	materialCells,
	materialRenderers,
} from "@jsonforms/material-renderers";
import { JsonForms } from "@jsonforms/react";
import type React from "react";
import { useCallback, useMemo } from "react";
import CustomColorPicker from "./jsonforms-renderers/CustomColorPicker";
// Import your custom renderers and testers
import CustomRangeSlider from "./jsonforms-renderers/CustomRangeSlider";
import CustomSmilesRenderer from "./jsonforms-renderers/CustomSmilesRenderer";
import { customColorPickerTester } from "./jsonforms-renderers/customColorPickerTester";
import { customRangeSliderTester } from "./jsonforms-renderers/customRangeSliderTester";
import { customSmilesRendererTester } from "./jsonforms-renderers/customSmilesRendererTester";

// TODO: do we need to memoize this?
const customRenderers = [
	...materialRenderers,
	{ tester: customRangeSliderTester, renderer: CustomRangeSlider },
	{ tester: customColorPickerTester, renderer: CustomColorPicker },
	{ tester: customSmilesRendererTester, renderer: CustomSmilesRenderer },
];

interface JSONFormsState {
	data: unknown;
	errors?: ErrorObject[];
}

interface Schema {
	data?: JsonSchema;
	ui?: UISchemaElement;
	[key: string]: unknown;
}

interface JSONFormsEditorProps {
	schema: Schema;
	data?: unknown;
	onChange: (data: unknown) => void;
	onValidationChange?: (errors: ErrorObject[]) => void;
}

export const JSONFormsEditor: React.FC<JSONFormsEditorProps> = ({
	schema,
	data,
	onChange,
	onValidationChange,
}) => {
	// Use Material-UI renderers
	const cells = useMemo(() => materialCells, []);

	const handleChange = useCallback(
		(state: JSONFormsState) => {
			console.log("JSONForms data changed:", state.data);
			onChange(state.data);

			if (onValidationChange) {
				onValidationChange(state.errors || []);
			}
		},
		[onChange, onValidationChange],
	);

	// Memoize the schema to prevent unnecessary re-renders
	const memoizedSchema = useMemo(() => {
		console.log("JSONForms schema updated:", schema);
		return schema?.data;
	}, [schema]);

	const memoizedUiSchema = useMemo(() => {
		console.log("JSONForms uischema updated:", schema);
		const ui = schema?.ui;
		if (Object.keys(ui || {}).length === 0) {
			return undefined;
		}
		return ui;
	}, [schema]);

	// Memoize the data to prevent unnecessary re-renders
	const memoizedData = useMemo(() => {
		console.log("JSONForms data updated:", data);
		return data || {};
	}, [data]);

	return (
		<div className="jsonforms-editor">
			<JsonForms
				schema={memoizedSchema}
				uischema={memoizedUiSchema}
				data={memoizedData}
				renderers={customRenderers}
				cells={cells}
				onChange={handleChange}
			/>
		</div>
	);
};

export default JSONFormsEditor;
