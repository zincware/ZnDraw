import { rankWith, schemaMatches } from "@jsonforms/core";

export const customSmilesRendererTester = rankWith(
	5, // High rank to override default string renderer
	schemaMatches((schema) => {
		return schema?.format === "smiles" && schema.type === "string";
	}),
);