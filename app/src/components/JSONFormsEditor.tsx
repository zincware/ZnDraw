import {
	materialCells,
	materialRenderers,
} from "@jsonforms/material-renderers";
import { JsonForms } from "@jsonforms/react";
import type React from "react";
import { useCallback, useMemo } from "react";
// Import your custom renderers and testers
import CustomRangeSlider from './jsonforms-renderers/CustomRangeSlider';
import { customRangeSliderTester } from './jsonforms-renderers/customRangeSliderTester';
import CustomColorPicker from './jsonforms-renderers/CustomColorPicker';
import { customColorPickerTester } from './jsonforms-renderers/customColorPickerTester';


// TODO: do we need to memoize this?
const customRenderers = [
  ...materialRenderers,
  { tester: customRangeSliderTester, renderer: CustomRangeSlider },
  { tester: customColorPickerTester, renderer: CustomColorPicker }
];
interface JSONFormsEditorProps {
	schema: any;
	data?: any;
	onChange: (data: any) => void;
	onValidationChange?: (errors: any[]) => void;
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
		(state: any) => {
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
		return schema?.ui;
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

