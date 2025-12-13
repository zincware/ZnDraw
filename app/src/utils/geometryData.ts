/**
 * Type definitions for geometry data props.
 * Backend uses strict list-only format:
 * - Position/Rotation/Scale: string (dynamic) | number[][] (Vec3 list)
 * - Color: string (dynamic) | string[]
 * - Size/Radius: string (dynamic) | number[]
 * Single-element lists broadcast to all instances.
 */

import { isTransform, type Transform } from "./transformProcessor";

// Fetched data types from server
// Uses any to handle TypedArray (Float32Array, etc.) since we call Array.from()
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export type FetchedNumeric = any;
// eslint-disable-next-line @typescript-eslint/no-explicit-any
export type FetchedStrings = any;

// Position is ALWAYS per-instance list of Vec3 or dynamic key or Transform
// Use number[][] instead of tuple for flexibility with JSON parsing
export type PositionProp = string | number[][] | Transform;

// Color is always hex strings - list of hex strings or dynamic reference or Transform
export type ColorProp = string | string[] | Transform;

// Size/radius is always list of floats or dynamic reference or Transform
export type SizeProp = string | number[] | Transform;

// Rotation is always list of Vec3 (Euler angles) or dynamic reference
export type RotationProp = string | number[][];

// Scale is always list of Vec3 or dynamic reference
export type ScaleProp = string | number[][];

// Size2D for planes: always list of [width, height] or dynamic reference
export type Size2DProp = string | number[][];

// Size3D for boxes: always list of [width, height, depth] or dynamic reference
export type Size3DProp = string | number[][];

// Generic data prop type (for processNumericAttribute)
export type DataProp =
	| string
	| number
	| number[]
	| number[][]
	| Transform
	| Record<string, unknown>;

/**
 * Constants for selection and hover visual effects.
 */
export const SELECTION_SCALE = 1.01;
export const HOVER_SCALE = 1.25;
export const SELECTION_COLOR = [1.0, 0.75, 0.8] as const; // Pink

/**
 * Process a Vec3 list attribute (position, rotation, scale, direction).
 * Backend always sends [[x,y,z], ...] format.
 * Single-element lists are broadcast to all instances.
 *
 * @param propValue - The prop value (string key or Vec3 list)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of values [x,y,z,x,y,z,...]
 */
export function processVec3Attribute(
	propValue: string | number[][],
	fetchedValue: FetchedNumeric,
	count: number,
): number[] {
	if (fetchedValue) {
		return Array.from(fetchedValue);
	}

	if (typeof propValue === "string") {
		// Dynamic reference not yet fetched
		return [];
	}

	// Backend sends [[x,y,z], ...] - broadcast if single element
	if (propValue.length === 1 && count > 1) {
		const [x, y, z] = propValue[0];
		const result: number[] = [];
		for (let i = 0; i < count; i++) {
			result.push(x, y, z);
		}
		return result;
	}

	return propValue.flat();
}

/**
 * Process position attribute.
 * Alias for processVec3Attribute for semantic clarity.
 */
export function processPositionAttribute(
	propValue: PositionProp,
	fetchedValue: FetchedNumeric,
	count = 0,
): number[] {
	// Transform should be evaluated before calling this
	if (isTransform(propValue)) {
		return [];
	}
	// For position, count is derived from the data itself
	const effectiveCount =
		count || (Array.isArray(propValue) ? propValue.length : 0);
	return processVec3Attribute(propValue, fetchedValue, effectiveCount);
}

/**
 * Process rotation attribute.
 * Backend sends [[rx,ry,rz], ...] Euler angles in radians.
 */
export function processRotationAttribute(
	propValue: RotationProp,
	fetchedValue: FetchedNumeric,
	count: number,
): number[] {
	return processVec3Attribute(propValue, fetchedValue, count);
}

/**
 * Process scale attribute.
 * Backend sends [[sx,sy,sz], ...] - always anisotropic.
 */
export function processScaleAttribute(
	propValue: ScaleProp,
	fetchedValue: FetchedNumeric,
	count: number,
): { values: number[]; isAnisotropic: boolean } {
	const values = processVec3Attribute(propValue, fetchedValue, count);
	// With strict list-only format, scale is always anisotropic (Vec3)
	return { values, isAnisotropic: true };
}

/**
 * Process a float list attribute (radius, etc.).
 * Backend sends [v, ...] format.
 * Single-element lists are broadcast to all instances.
 *
 * @param propValue - The prop value (string key or float list)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Array of values [v1, v2, ...]
 */
export function processFloatAttribute(
	propValue: SizeProp,
	fetchedValue: FetchedNumeric,
	count: number,
): number[] {
	if (fetchedValue) {
		return Array.from(fetchedValue);
	}

	// Transform should be evaluated before calling this
	if (isTransform(propValue)) {
		return [];
	}

	if (typeof propValue === "string") {
		// Dynamic reference not yet fetched
		return [];
	}

	// Backend sends [v, ...] - broadcast if single element
	if (propValue.length === 1 && count > 1) {
		return Array(count).fill(propValue[0]);
	}

	return propValue;
}

/**
 * Process a 2D size attribute (Plane geometry: width, height).
 * Backend sends [[w, h], ...] format.
 * Single-element lists are broadcast to all instances.
 *
 * @param propValue - The size prop value (string key or Vec2 list)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of 2D size values [w1,h1, w2,h2, ...]
 */
export function processSize2D(
	propValue: Size2DProp,
	fetchedValue: FetchedNumeric,
	count: number,
): number[] {
	if (fetchedValue) {
		return Array.from(fetchedValue);
	}

	if (typeof propValue === "string") {
		return [];
	}

	// Backend sends [[w,h], ...] - broadcast if single element
	if (propValue.length === 1 && count > 1) {
		const [w, h] = propValue[0];
		const result: number[] = [];
		for (let i = 0; i < count; i++) {
			result.push(w, h);
		}
		return result;
	}

	return propValue.flat();
}

/**
 * Process a 3D size attribute (Box geometry: width, height, depth).
 * Backend sends [[w, h, d], ...] format.
 * Single-element lists are broadcast to all instances.
 *
 * @param propValue - The size prop value (string key or Vec3 list)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of 3D size values [w1,h1,d1, w2,h2,d2, ...]
 */
export function processSize3D(
	propValue: Size3DProp,
	fetchedValue: FetchedNumeric,
	count: number,
): number[] {
	if (fetchedValue) {
		return Array.from(fetchedValue);
	}

	if (typeof propValue === "string") {
		return [];
	}

	// Backend sends [[w,h,d], ...] - broadcast if single element
	if (propValue.length === 1 && count > 1) {
		const [w, h, d] = propValue[0];
		const result: number[] = [];
		for (let i = 0; i < count; i++) {
			result.push(w, h, d);
		}
		return result;
	}

	return propValue.flat();
}

/**
 * Process color data - returns array of hex strings.
 * Backend sends ["#hex", ...] format.
 * Single-element lists are broadcast to all instances.
 *
 * @param propValue - The color prop value from backend
 * @param fetchedValue - Data fetched from server (if propValue is dynamic ref)
 * @param count - Expected number of instances
 * @returns Array of hex color strings (empty if dynamic ref not yet fetched)
 */
export function processColorData(
	propValue: ColorProp | Record<string, unknown>,
	fetchedValue: FetchedStrings,
	count: number,
): string[] {
	if (fetchedValue) {
		return Array.from(fetchedValue);
	}

	if (Array.isArray(propValue)) {
		// Broadcast if single element
		if (propValue.length === 1 && count > 1) {
			return Array(count).fill(propValue[0]);
		}
		return propValue;
	}

	// Dynamic reference or transform not yet fetched/evaluated
	return [];
}

/**
 * Legacy function for backwards compatibility.
 * Use processVec3Attribute or processFloatAttribute instead.
 */
export function processNumericAttribute(
	propValue: DataProp,
	fetchedValue: FetchedNumeric,
	count: number,
): number[] {
	if (fetchedValue) {
		return Array.from(fetchedValue);
	}

	if (typeof propValue === "number") {
		return Array(count).fill(propValue);
	}

	if (typeof propValue === "string") {
		return [];
	}

	if (Array.isArray(propValue) && propValue.length > 0) {
		const firstElem = propValue[0];
		if (Array.isArray(firstElem)) {
			// [[x,y,z], ...] -> broadcast if single element
			const tuples = propValue as number[][];
			if (tuples.length === 1 && count > 1) {
				const result: number[] = [];
				for (let i = 0; i < count; i++) {
					result.push(...tuples[0]);
				}
				return result;
			}
			return tuples.flat();
		}
		// [v, ...] - assume floats
		const values = propValue as number[];
		if (values.length === 1 && count > 1) {
			return Array(count).fill(values[0]);
		}
		return values;
	}

	return [];
}

/**
 * Validate that all attribute arrays have consistent lengths.
 *
 * @param arrays - Object mapping attribute names to their arrays
 * @param expectedCounts - Object mapping attribute names to their expected element count
 * @returns True if all arrays have correct lengths
 */
export function validateArrayLengths(
	arrays: Record<string, number[]>,
	expectedCounts: Record<string, number>,
): boolean {
	for (const [name, array] of Object.entries(arrays)) {
		const expectedCount = expectedCounts[name];
		if (array.length !== expectedCount) {
			console.error(
				`${name} has incorrect length: expected ${expectedCount}, got ${array.length}`,
			);
			return false;
		}
	}
	return true;
}

/**
 * Determine instance count from position attribute.
 * Backend always sends [[x,y,z], ...] format after normalization.
 *
 * @param positionProp - The position prop value
 * @param fetchedPosition - Fetched position data (if applicable)
 * @returns Number of instances (0 if invalid or not yet fetched)
 */
export function getInstanceCount(
	positionProp: PositionProp | Record<string, unknown>,
	fetchedPosition: FetchedNumeric,
): number {
	if (fetchedPosition) {
		if (fetchedPosition.length % 3 !== 0) {
			return 0;
		}
		return fetchedPosition.length / 3;
	}

	if (Array.isArray(positionProp)) {
		return positionProp.length;
	}

	return 0;
}

/**
 * Get color count to detect single-color mode for UI.
 *
 * @param colorProp - The color prop value
 * @returns Number of colors, or 0 if dynamic reference or Transform
 */
export function getColorCount(colorProp: ColorProp): number {
	if (typeof colorProp === "string" || isTransform(colorProp)) {
		return 0; // Dynamic reference or Transform
	}
	return colorProp.length;
}

/**
 * Expand a shared color (single color for all instances) to a full array.
 * If the color array has only one color and there are multiple instances,
 * replicate that color for all instances.
 *
 * @param colorHexArray - Array of hex color strings
 * @param count - Number of instances
 * @returns Array of hex colors (expanded if shared, unchanged otherwise)
 */
export function expandSharedColor(
	colorHexArray: string[],
	count: number,
): string[] {
	return colorHexArray.length === 1 && count > 1
		? Array(count).fill(colorHexArray[0])
		: colorHexArray;
}
