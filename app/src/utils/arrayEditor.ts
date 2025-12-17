/**
 * Utility functions for array editing in geometry forms.
 * All configuration is derived from JSON schema - no hardcoded field types.
 */

/**
 * Information extracted from JSON schema for array fields
 */
export interface ArrayFieldInfo {
	/** Number of values per row (e.g., 3 for Vec3, 2 for Vec2, 1 for scalar) */
	dimensions: number;
	/** Whether items are strings (e.g., hex colors) vs numbers */
	isStringType: boolean;
	/** Minimum value for numeric items (from schema) */
	itemMinimum?: number;
	/** Maximum value for numeric items (from schema) */
	itemMaximum?: number;
}

/**
 * Extract array field information from JSON schema.
 * Parses the anyOf structure to find the array type and its item structure.
 *
 * Schema patterns:
 * - Vec3: items: {prefixItems: [{type: "number"}, ...], minItems: 3, maxItems: 3}
 * - Vec2: items: {prefixItems: [{type: "number"}, ...], minItems: 2, maxItems: 2}
 * - list[float]: items: {type: "number"}
 * - list[str]: items: {type: "string"}
 */
export function getArrayFieldInfo(schema: any): ArrayFieldInfo | null {
	if (!schema) return null;

	// Find the array type option in anyOf
	const arrayOption = schema.anyOf?.find(
		(opt: any) => opt.type === "array" && opt.items,
	);
	if (!arrayOption) return null;

	const items = arrayOption.items;
	if (!items) return null;

	// Check if items is a tuple (has prefixItems - Vec2, Vec3, etc.)
	if (items.prefixItems && items.prefixItems.length > 0) {
		const firstItem = items.prefixItems[0];
		return {
			dimensions: items.prefixItems.length,
			isStringType: firstItem?.type === "string",
			itemMinimum: firstItem?.minimum,
			itemMaximum: firstItem?.maximum,
		};
	}

	// Check for minItems/maxItems (alternative tuple definition)
	if (items.minItems !== undefined && items.minItems > 1) {
		return {
			dimensions: items.minItems,
			isStringType: false,
			itemMinimum: undefined,
			itemMaximum: undefined,
		};
	}

	// Simple array (list[float] or list[str]) - 1 dimension
	return {
		dimensions: 1,
		isStringType: items.type === "string",
		itemMinimum: items.minimum,
		itemMaximum: items.maximum,
	};
}

/**
 * Get column labels for the data grid.
 * Uses simple defaults based on dimensions.
 */
export function getColumnLabels(dimensions: number): string[] {
	if (dimensions === 1) return ["Value"];
	if (dimensions === 2) return ["X", "Y"];
	if (dimensions === 3) return ["X", "Y", "Z"];
	return Array.from({ length: dimensions }, (_, i) => `Col ${i + 1}`);
}

/**
 * Normalize value to 2D array format for editing.
 * Handles: undefined, number, number[], number[][], string, string[]
 * Returns empty array for undefined (required fields without defaults).
 */
export function normalizeToArray(
	value:
		| string
		| number
		| (string | number)[]
		| (string | number)[][]
		| undefined,
	schema: any,
): (string | number)[][] {
	const fieldInfo = getArrayFieldInfo(schema);
	if (!fieldInfo) {
		// Fallback: if we can't parse schema, return empty or wrap value
		if (value === undefined) return [];
		if (Array.isArray(value) && Array.isArray(value[0]))
			return value as (string | number)[][];
		if (Array.isArray(value)) return [value as (string | number)[]];
		return [[value as string | number]];
	}

	const { dimensions, isStringType } = fieldInfo;

	// Return empty array for undefined (required fields without schema defaults)
	if (value === undefined) {
		return [];
	}

	// For string types (like hex colors)
	if (isStringType) {
		// If it's already a 2D array of strings, return as-is
		if (Array.isArray(value) && Array.isArray(value[0])) {
			return value as string[][];
		}

		// If it's a 1D array of strings (e.g., ["#FF0000", "#00FF00"])
		if (Array.isArray(value)) {
			return (value as string[]).map((item) => [item]);
		}

		// If it's a single string (e.g., "#FF0000")
		if (typeof value === "string") {
			return [[value]];
		}

		// Malformed data - return empty array
		return [];
	}

	// For numeric types
	// If it's already a 2D array, return as-is
	if (Array.isArray(value) && Array.isArray(value[0])) {
		return value as number[][];
	}

	// If it's a 1D array
	if (Array.isArray(value)) {
		const arr = value as number[];

		// If dimensions match, wrap in array (single row)
		if (arr.length === dimensions) {
			return [arr];
		}
		// Otherwise, chunk by dimensions
		const result: number[][] = [];
		for (let i = 0; i < arr.length; i += dimensions) {
			result.push(arr.slice(i, i + dimensions));
		}
		return result;
	}

	// If it's a single number, create single row
	if (typeof value === "number") {
		return [Array(dimensions).fill(value)];
	}

	// Malformed data - return empty array
	return [];
}

/**
 * Convert 2D array back to appropriate format based on schema.
 * Always returns list format to match Pydantic types.
 */
export function denormalizeFromArray(
	arrayValue: (string | number)[][],
	schema: any,
): (string | number)[] | (string | number)[][] {
	const fieldInfo = getArrayFieldInfo(schema);

	// Empty array - return empty list
	if (arrayValue.length === 0) {
		return [];
	}

	// For string types (like hex colors) - return 1D array
	if (fieldInfo?.isStringType) {
		return arrayValue.map((row) => row[0] as string);
	}

	// For 1D numeric fields (like radius) - return 1D array
	if (fieldInfo?.dimensions === 1) {
		return arrayValue.map((row) => row[0] as number);
	}

	// For multi-dimensional fields (Vec2, Vec3) - return 2D array
	return arrayValue as number[][];
}

/**
 * Validate array data against schema.
 */
export function validateArrayData(
	arrayValue: (string | number)[][],
	schema: any,
): { valid: boolean; errors: string[] } {
	const errors: string[] = [];
	const fieldInfo = getArrayFieldInfo(schema);

	if (!fieldInfo) {
		return { valid: true, errors: [] };
	}

	const { dimensions, isStringType, itemMinimum, itemMaximum } = fieldInfo;
	const columnLabels = getColumnLabels(dimensions);

	// Check all rows have correct dimensions
	arrayValue.forEach((row, idx) => {
		if (row.length !== dimensions) {
			errors.push(
				`Row ${idx + 1} has ${row.length} values, expected ${dimensions}`,
			);
		}

		// For string types (like hex colors), validate format
		if (isStringType) {
			row.forEach((val, colIdx) => {
				if (typeof val !== "string") {
					errors.push(
						`Row ${idx + 1}, ${columnLabels[colIdx]}: expected string, got ${typeof val}`,
					);
				} else if (!isValidHexColor(val)) {
					errors.push(
						`Row ${idx + 1}, ${columnLabels[colIdx]}: invalid hex color "${val}"`,
					);
				}
			});
		} else {
			// Check value ranges for numeric types
			if (itemMinimum !== undefined || itemMaximum !== undefined) {
				row.forEach((val, colIdx) => {
					const numVal = val as number;
					if (itemMinimum !== undefined && numVal < itemMinimum) {
						errors.push(
							`Row ${idx + 1}, ${columnLabels[colIdx]}: value ${numVal} below minimum ${itemMinimum}`,
						);
					}
					if (itemMaximum !== undefined && numVal > itemMaximum) {
						errors.push(
							`Row ${idx + 1}, ${columnLabels[colIdx]}: value ${numVal} above maximum ${itemMaximum}`,
						);
					}
				});
			}
		}
	});

	return {
		valid: errors.length === 0,
		errors,
	};
}

/**
 * Validate hex color format
 */
function isValidHexColor(color: string): boolean {
	return /^#[0-9A-Fa-f]{6}$/.test(color);
}

/**
 * Create a new row by duplicating the last existing row,
 * or deriving structure from schema default.
 */
export function createDefaultRow(
	schema: any,
	existingRows?: (string | number)[][],
): (string | number)[] {
	// If there are existing rows, copy the last one
	if (existingRows && existingRows.length > 0) {
		return [...existingRows[existingRows.length - 1]];
	}

	// Try to derive from schema default
	if (schema?.default && Array.isArray(schema.default)) {
		const defaultValue = schema.default;
		// If default is 2D array, use first row
		if (Array.isArray(defaultValue[0])) {
			return [...defaultValue[0]];
		}
		// If default is 1D array (like radius [1.0]), wrap each item
		return [defaultValue[0]];
	}

	// Fallback: use schema info to create empty row
	const fieldInfo = getArrayFieldInfo(schema);
	if (fieldInfo) {
		const { dimensions, isStringType } = fieldInfo;
		if (isStringType) {
			return Array(dimensions).fill("#000000");
		}
		return Array(dimensions).fill(0);
	}

	// Ultimate fallback
	return [0];
}

/**
 * Parse clipboard data (supports CSV, TSV, JSON)
 */
export function parseClipboardData(text: string): number[][] | null {
	try {
		// Try JSON first
		const parsed = JSON.parse(text);
		if (Array.isArray(parsed)) {
			return parsed;
		}
	} catch {
		// Not JSON, try CSV/TSV
		const lines = text.trim().split("\n");
		const result: number[][] = [];

		for (const line of lines) {
			// Try tab-separated, then comma-separated
			const values = line.includes("\t") ? line.split("\t") : line.split(",");

			const numbers = values
				.map((v) => parseFloat(v.trim()))
				.filter((n) => !isNaN(n));

			if (numbers.length > 0) {
				result.push(numbers);
			}
		}

		if (result.length > 0) {
			return result;
		}
	}

	return null;
}

/**
 * Get human-readable label for array shape.
 * Simple version that just shows row count.
 */
export function getArrayShapeLabel(
	value: (string | number)[] | (string | number)[][],
	schema: any,
): string {
	const fieldInfo = getArrayFieldInfo(schema);

	// Convert to 2D array format for consistent handling
	let rows: (string | number)[][];
	if (Array.isArray(value[0])) {
		rows = value as (string | number)[][];
	} else {
		// Flat array - chunk by dimensions
		const flatArray = value as (string | number)[];
		const dimensions = fieldInfo?.dimensions ?? 1;
		rows = [];
		for (let i = 0; i < flatArray.length; i += dimensions) {
			rows.push(flatArray.slice(i, i + dimensions));
		}
	}

	const count = rows.length;
	const word = count === 1 ? "row" : "rows";

	return `${count} ${word}`;
}

/**
 * Get preview string for array tooltip.
 */
export function getArrayPreview(
	value: (string | number)[] | (string | number)[][],
	maxItems: number = 2,
): string {
	if (Array.isArray(value[0])) {
		const rows = value as (string | number)[][];
		const preview = rows
			.slice(0, maxItems)
			.map(
				(row) =>
					`[${row.map((v) => (typeof v === "number" ? v.toFixed(2) : v)).join(", ")}]`,
			)
			.join(", ");
		return rows.length > maxItems ? `${preview}, ...` : preview;
	} else {
		const flatArray = value as (string | number)[];
		const preview = flatArray
			.slice(0, maxItems * 3)
			.map((v) => (typeof v === "number" ? v.toFixed(2) : v))
			.join(", ");
		return flatArray.length > maxItems * 3 ? `${preview}, ...` : preview;
	}
}
