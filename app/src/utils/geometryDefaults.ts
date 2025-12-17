import { merge } from "lodash";

/**
 * Utility functions for handling geometry defaults from Pydantic models.
 * Follows Single Responsibility Principle - Pydantic is the single source of truth for defaults.
 */

/**
 * Extracts default values from a JSON schema.
 * Traverses the schema's properties and collects the 'default' value for each field.
 *
 * @param schema - The JSON schema object (from Pydantic model)
 * @returns An object with field names as keys and their default values
 */
function extractDefaultsFromSchema(
	schema: Record<string, any>,
): Record<string, any> {
	const defaults: Record<string, any> = {};

	if (!schema?.properties) {
		return defaults;
	}

	for (const [key, propSchema] of Object.entries(schema.properties)) {
		const prop = propSchema as Record<string, any>;
		if (prop.default !== undefined) {
			defaults[key] = prop.default;
		} else if (prop.items?.default !== undefined) {
			// Handle array schema defaults (e.g., list[Vec3] with item defaults)
			defaults[key] = prop.items.default;
		} else if (prop.properties) {
			// Handle nested objects (like InteractionSettings)
			defaults[key] = extractDefaultsFromSchema(prop);
		}
	}

	return defaults;
}

/**
 * Merges geometry data with defaults from Pydantic models.
 *
 * This ensures that all required fields have values, even if they weren't sent by the API.
 * Pydantic models are the single source of truth for defaults.
 *
 * @param data - The geometry data from the server (may be partial)
 * @param type - The geometry type (e.g., "Arrow", "Box", "Sphere")
 * @param schemas - The geometry schemas object from Zustand store (JSON schemas)
 * @returns Complete geometry data with all defaults applied
 *
 * @example
 * ```ts
 * const arrowData = getGeometryWithDefaults(
 *   { position: "arrays.positions", direction: "calc.forces" },
 *   "Arrow",
 *   geometrySchemas
 * );
 * // Result will include opacity, selecting, hovering from Arrow defaults
 * ```
 */
export function getGeometryWithDefaults<T extends Record<string, any>>(
	data: Partial<T>,
	type: string,
	schemas: Record<string, any>,
): T {
	// If schemas aren't loaded yet (during app initialization), use data as-is
	if (!schemas || Object.keys(schemas).length === 0) {
		return data as T;
	}

	const schemaForType = schemas[type];

	// If this specific geometry type has no schema, use data as-is
	if (!schemaForType) {
		return data as T;
	}

	// Extract default values from the JSON schema
	const defaultsForType = extractDefaultsFromSchema(schemaForType);

	// Deep merge: defaults first, then override with actual data
	// lodash merge mutates the first argument, so pass empty object
	return merge({}, defaultsForType, data) as T;
}
