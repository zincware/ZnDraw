import { shouldFetchAsFrameData, hexToRgb } from "./colorUtils";

/**
 * Type definitions for geometry data props.
 * Props can be static values or dynamic keys that reference server data.
 */
export type StaticValue = number | number[] | number[][];
export type DataProp = string | StaticValue;

/**
 * Constants for selection and hover visual effects.
 */
export const SELECTION_SCALE = 1.01;
export const HOVER_SCALE = 1.25;
export const SELECTION_COLOR = [1.0, 0.75, 0.8] as const; // Pink

/**
 * Build list of keys that need to be fetched from the server.
 * Filters out static values and hex color strings.
 * Automatically deduplicates keys (e.g., if both position and direction are "arrays.positions").
 *
 * @param props - Object containing all component props
 * @param skipColorCheck - If true, includes all string color props (default: false)
 * @returns Deduplicated array of keys to fetch from server
 */
export function buildFetchKeys(
  props: Record<string, DataProp>,
  skipColorCheck: boolean = false
): string[] {
  const keys: string[] = [];

  for (const [key, value] of Object.entries(props)) {
    if (typeof value !== "string") continue;

    // Special handling for color props - skip hex colors unless forced
    if (key === "color" && !skipColorCheck && !shouldFetchAsFrameData(value)) {
      continue;
    }

    keys.push(value);
  }

  // Deduplicate keys - multiple props might reference the same data
  // e.g., position="arrays.positions" and direction="arrays.positions"
  return Array.from(new Set(keys));
}

/**
 * Process a numeric attribute (position, radius, scale, direction, etc.).
 * Handles both fetched data and static values.
 *
 * @param propValue - The prop value (string key or static data)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of numeric values
 */
export function processNumericAttribute(
  propValue: DataProp,
  fetchedValue: number[] | Float32Array | Float64Array | undefined,
  count: number
): number[] {
  let finalArray: number[] = [];

  if (fetchedValue) {
    finalArray = Array.from(fetchedValue);
  } else if (typeof propValue !== "string") {
    // Handle static values
    if (Array.isArray(propValue[0])) {
      // Array of arrays/vectors (e.g., [[1,2,3], [4,5,6]])
      finalArray = (propValue as number[][]).flat();
    } else if (Array.isArray(propValue) && propValue.length > 0) {
      // Check if it's a single vector or array to replicate
      const firstElem = propValue[0];
      if (typeof firstElem === "number") {
        // Could be single value [1,2,3] or value to replicate [5]
        if (propValue.length === 1 || propValue.length === 3) {
          // Single value to replicate for all instances
          finalArray = Array(count).fill(propValue).flat();
        } else {
          // Already a flat array
          finalArray = propValue as number[];
        }
      }
    } else {
      // Single number to replicate
      finalArray = Array(count).fill(propValue).flat();
    }
  }

  return finalArray;
}

/**
 * Process a color attribute with special handling for hex colors.
 * Handles hex strings, RGB arrays, and fetched data.
 *
 * @param propValue - The color prop value (hex string, RGB array, or fetch key)
 * @param fetchedValue - Data fetched from server (if propValue is a fetch key)
 * @param count - Expected number of instances
 * @returns Flattened array of RGB values (0-1 range)
 * @throws Error if hex color is invalid
 */
export function processColorAttribute(
  propValue: DataProp,
  fetchedValue: number[] | Float32Array | Float64Array | undefined,
  count: number
): number[] {
  let finalColors: number[] = [];

  if (fetchedValue) {
    // Data was fetched from server
    finalColors = Array.from(fetchedValue);
  } else if (typeof propValue === "string") {
    // Hex color string - convert to RGB and replicate for all instances
    const rgb = hexToRgb(propValue);
    if (rgb) {
      finalColors = Array(count).fill(rgb).flat();
    } else {
      throw new Error(`Invalid hex color: ${propValue}`);
    }
  } else if (Array.isArray(propValue)) {
    // Static color data
    if (Array.isArray(propValue[0])) {
      // Array of RGB arrays
      finalColors = (propValue as number[][]).flat();
    } else {
      // Single RGB value to replicate
      finalColors = Array(count).fill(propValue).flat();
    }
  }

  return finalColors;
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
  expectedCounts: Record<string, number>
): boolean {
  for (const [name, array] of Object.entries(arrays)) {
    const expectedCount = expectedCounts[name];
    if (array.length !== expectedCount) {
      console.error(
        `${name} has incorrect length: expected ${expectedCount}, got ${array.length}`
      );
      return false;
    }
  }
  return true;
}

/**
 * Determine instance count from position attribute.
 *
 * @param positionProp - The position prop value
 * @param fetchedPosition - Fetched position data (if applicable)
 * @returns Number of instances
 */
export function getInstanceCount(
  positionProp: DataProp,
  fetchedPosition: number[] | Float32Array | Float64Array | undefined
): number {
  if (fetchedPosition) {
    return fetchedPosition.length / 3;
  } else if (typeof positionProp !== "string") {
    if (Array.isArray(positionProp[0])) {
      return (positionProp as number[][]).length;
    } else if (Array.isArray(positionProp)) {
      return (positionProp as number[]).length / 3;
    } else {
      return 1;
    }
  }
  return 0;
}
