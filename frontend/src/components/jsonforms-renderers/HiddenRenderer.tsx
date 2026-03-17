import { rankWith, schemaMatches } from "@jsonforms/core";
import { withJsonFormsControlProps } from "@jsonforms/react";

/**
 * Renderer that hides fields marked with x-hidden: true in the schema.
 * The field still exists in the form data, it's just not shown.
 */
const HiddenRenderer = () => null;

export const hiddenRendererTester = rankWith(
	100, // Highest priority — must override all other renderers
	schemaMatches((schema) => {
		const ext = schema as Record<string, unknown>;
		return ext?.["x-hidden"] === true;
	}),
);

export default withJsonFormsControlProps(HiddenRenderer);
