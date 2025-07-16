// src/testers/customColorPickerTester.js
import { rankWith, schemaMatches } from "@jsonforms/core";

export const customColorPickerTester = rankWith(
	5, // A high rank, same as the other custom one
	schemaMatches((schema) => {
		return schema?.format === "color" && schema.type === "string";
	}),
);
