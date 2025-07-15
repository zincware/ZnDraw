import {
	materialCells,
	materialRenderers,
} from "@jsonforms/material-renderers";
import { JsonForms } from "@jsonforms/react";
import type React from "react";
import { useCallback, useMemo } from "react";

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
	const renderers = useMemo(() => materialRenderers, []);
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
		return schema;
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
				data={memoizedData}
				renderers={renderers}
				cells={cells}
				onChange={handleChange}
			/>
		</div>
	);
};

export default JSONFormsEditor;
