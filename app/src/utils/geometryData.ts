import { hexToRgb } from "./colorUtils";

/**
 * Type definitions for geometry data props.
 * Updated to match backend type system.
 */
// Position is ALWAYS per-instance (list of tuples or dynamic key)
export type PositionProp = string | [number, number, number][];

// Color can be shared (single hex) or per-instance (list of hex) or dynamic
export type ColorProp = string | string[];

// Size/radius can be shared (single value) or per-instance (list) or dynamic
export type SizeProp = string | number | number[];

// Rotation can be shared (single tuple) or per-instance (list of tuples) or dynamic
export type RotationProp = string | [number, number, number] | [number, number, number][];

/**
 * Constants for selection and hover visual effects.
 */
export const SELECTION_SCALE = 1.01;
export const HOVER_SCALE = 1.25;
export const SELECTION_COLOR = [1.0, 0.75, 0.8] as const; // Pink

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
  fetchedValue: any,
  count: number
): number[] {
  let finalArray: number[] = [];

  if (fetchedValue) {
    finalArray = Array.from(fetchedValue);
  } else if (typeof propValue !== "string" && typeof propValue !== "number") {
    // Handle static values (arrays)
    if (Array.isArray(propValue) && propValue.length > 0) {
      const firstElem = propValue[0];
      if (Array.isArray(firstElem)) {
        // Array of arrays/vectors (e.g., [[1,2,3], [4,5,6]])
        finalArray = (propValue as number[][]).flat();
      } else if (typeof firstElem === "number") {
        // Check if it's a single vector or array to replicate
        if (propValue.length === 1 || propValue.length === 3) {
          // Single value to replicate for all instances
          finalArray = Array(count).fill([...propValue]).flat();
        } else {
          // Already a flat array
          finalArray = propValue as number[];
        }
      }
    }
  } else if (typeof propValue === "number") {
    // Single number to replicate
    finalArray = Array(count).fill(propValue);
  }

  return finalArray;
}

/**
 * Process a color attribute with special handling for hex colors.
 * Handles hex strings (shared), list of hex strings (per-instance), RGB arrays, and fetched data.
 *
 * @param propValue - The color prop value (hex string, list of hex strings, RGB array, or fetch key)
 * @param fetchedValue - Data fetched from server (if propValue is a fetch key)
 * @param count - Expected number of instances
 * @returns Flattened array of RGB values (0-1 range)
 * @throws Error if hex color is invalid
 */
export function processColorAttribute(
  propValue: ColorProp,
  fetchedValue: any,
  count: number
): number[] {
  let finalColors: number[] = [];

  if (fetchedValue) {
    // Data was fetched from server
    finalColors = Array.from(fetchedValue);
  } else if (typeof propValue === "string") {
    // Shared hex color string - convert to RGB and replicate for all instances
    const rgb = hexToRgb(propValue);
    if (rgb) {
      finalColors = Array(count).fill(rgb).flat();
    } else {
      throw new Error(`Invalid hex color: ${propValue}`);
    }
  } else if (Array.isArray(propValue)) {
    // Check if it's an array of hex strings or RGB arrays
    if (typeof propValue[0] === "string") {
      // Array of hex strings - convert each to RGB
      finalColors = (propValue as string[]).flatMap(hex => {
        const rgb = hexToRgb(hex);
        if (!rgb) throw new Error(`Invalid hex color: ${hex}`);
        return rgb;
      });
    } else if (Array.isArray(propValue[0])) {
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
 * Process a size/radius attribute (1D - single value per instance).
 * Handles shared (single value for all) or per-instance (list) or dynamic data.
 *
 * @param propValue - The size prop value (string key, number, or list)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of size values [r1, r2, r3, ...]
 */
export function processSizeAttribute(
  propValue: SizeProp,
  fetchedValue: any,
  count: number
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (typeof propValue === "number") {
    // Shared value - replicate for all instances
    return Array(count).fill(propValue);
  }

  if (Array.isArray(propValue)) {
    // Per-instance list
    return propValue as number[];
  }

  return [];
}

/**
 * Process a 2D size attribute (Plane geometry: width, height).
 * Handles shared (single tuple for all) or per-instance (list of tuples) or dynamic data.
 *
 * @param propValue - The size prop value (string key, tuple, or list of tuples)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of 2D size values [w1,h1, w2,h2, ...]
 */
export function processSize2D(
  propValue: string | number[] | number[][],
  fetchedValue: any,
  count: number
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (typeof propValue !== "string" && Array.isArray(propValue)) {
    // Check if it's a list of tuples or a single tuple
    const firstElem = propValue[0];
    if (Array.isArray(firstElem)) {
      // List of tuples [[w,h], [w,h], ...]
      return (propValue as number[][]).flat();
    } else if (propValue.length === 2 && typeof firstElem === "number") {
      // Single tuple [w, h] - replicate for all instances
      return Array(count).fill([...propValue]).flat();
    }
  }

  return [];
}

/**
 * Process a 3D size attribute (Box geometry: width, height, depth).
 * Handles shared (single tuple for all) or per-instance (list of tuples) or dynamic data.
 *
 * @param propValue - The size prop value (string key, tuple, or list of tuples)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of 3D size values [w1,h1,d1, w2,h2,d2, ...]
 */
export function processSize3D(
  propValue: string | number[] | number[][],
  fetchedValue: any,
  count: number
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (typeof propValue !== "string" && Array.isArray(propValue)) {
    // Check if it's a list of tuples or a single tuple
    const firstElem = propValue[0];
    if (Array.isArray(firstElem)) {
      // List of tuples [[w,h,d], [w,h,d], ...]
      return (propValue as number[][]).flat();
    } else if (propValue.length === 3 && typeof firstElem === "number") {
      // Single tuple [w, h, d] - replicate for all instances
      return Array(count).fill([...propValue]).flat();
    }
  }

  return [];
}

/**
 * Process a rotation attribute.
 * Handles shared (single tuple for all) or per-instance (list of tuples) or dynamic data.
 *
 * @param propValue - The rotation prop value (string key, tuple, or list of tuples)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @param count - Expected number of instances
 * @returns Flattened array of rotation values [x,y,z,x,y,z,...]
 */
export function processRotationAttribute(
  propValue: string | number[] | number[][],
  fetchedValue: any,
  count: number
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (Array.isArray(propValue)) {
    // Check if it's a single tuple or list of tuples
    if (typeof propValue[0] === "number") {
      // Single tuple [x, y, z] - replicate for all instances
      return Array(count).fill([...propValue]).flat();
    } else {
      // List of tuples [[x,y,z], [x,y,z], ...] - flatten
      return (propValue as number[][]).flat();
    }
  }

  return [];
}

/**
 * Process position attribute (always per-instance).
 *
 * @param propValue - The position prop value (string key, list of tuples, or flat array)
 * @param fetchedValue - Data fetched from server (if propValue is a string)
 * @returns Flattened array of position values [x,y,z,x,y,z,...]
 */
export function processPositionAttribute(
  propValue: string | number[][] | number[],
  fetchedValue: any
): number[] {
  if (fetchedValue) {
    return Array.from(fetchedValue);
  }

  if (typeof propValue !== "string" && Array.isArray(propValue)) {
    // Check if it's already flat or needs flattening
    if (Array.isArray(propValue[0])) {
      // List of tuples [[x,y,z], [x,y,z], ...]
      return propValue.flat();
    }
    // Already flat array [x,y,z,x,y,z,...]
    return propValue as number[];
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
 * Position is always per-instance (list of tuples or dynamic key).
 *
 * @param positionProp - The position prop value
 * @param fetchedPosition - Fetched position data (if applicable)
 * @returns Number of instances
 */
export function getInstanceCount(
  positionProp: string | number[] | number[][],
  fetchedPosition: any
): number {
  if (fetchedPosition) {
    return fetchedPosition.length / 3;
  }

  if (Array.isArray(positionProp)) {
    // Check if it's a flat array or array of tuples
    if (Array.isArray(positionProp[0])) {
      // List of tuples [[x,y,z], [x,y,z], ...]
      return positionProp.length;
    } else {
      // Flat array [x,y,z,x,y,z,...] - divide by 3
      return positionProp.length / 3;
    }
  }

  // String (dynamic key) - unknown until fetched
  return 0;
}
