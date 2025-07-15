import React, { useMemo, useCallback } from "react";
import { JsonForms } from "@jsonforms/react";
import { vanillaRenderers, vanillaCells } from "@jsonforms/vanilla-renderers";

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
	// Use vanilla renderers for bootstrap compatibility
	const renderers = useMemo(() => vanillaRenderers, []);
	const cells = useMemo(() => vanillaCells, []);

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
