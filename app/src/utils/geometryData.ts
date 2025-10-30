/**
 * Type definitions for geometry data props.
 * Updated to match backend normalization - all non-dynamic fields are lists.
 */
// Position is ALWAYS per-instance (list of tuples or dynamic key)
export type PositionProp = string | [number, number, number][];

// Color is always hex strings - list of hex strings or dynamic reference
export type ColorProp = string | string[];

// Size/radius can be shared (single value) or per-instance (list) or dynamic
export type SizeProp = string | number | number[];

// Rotation can be shared (single tuple) or per-instance (list of tuples) or dynamic
export type RotationProp = string | [number, number, number] | [number, number, number][];

// Generic data prop type (for processNumericAttribute)
export type DataProp = string | number | number[] | number[][] | Record<string, any>;

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
 * Process color data - returns array of hex strings
 * Backend sends: string (dynamic ref) | string[] (hex list)
 * Fetched data: string[] (hex list from server)
 *
 * @param propValue - The color prop value from backend
 * @param fetchedValue - Data fetched from server (if propValue is dynamic ref)
 * @param count - Expected number of instances
 * @returns Array of hex color strings
 */
export function processColorData(
  propValue: string | string[] | Record<string, any>,
  fetchedValue: any,
  count: number
): string[] {
  // Fetched data is already hex strings from backend
  if (fetchedValue) {
    return fetchedValue as string[];
  }

  // Backend always sends list of hex strings [\"#FF0000\", ...]
  if (Array.isArray(propValue)) {
    return propValue;
  }

  // Dynamic ref or Transform not fetched - throw error
  if (typeof propValue === 'string') {
    throw new Error(`Dynamic color reference not fetched: ${propValue}`);
  }

  // Transform case - should be evaluated before calling this function
  throw new Error(`Transform not evaluated before processing color data`);
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
      // Check if it's a single tuple that needs replication
      if ((propValue as number[][]).length === 1 && count > 1) {
        // Single tuple in a list [(w,h)] - replicate for all instances
        return Array(count).fill([...(propValue as number[][])[0]]).flat();
      }
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
      // Check if it's a single tuple that needs replication
      if ((propValue as number[][]).length === 1 && count > 1) {
        // Single tuple in a list [(w,h,d)] - replicate for all instances
        return Array(count).fill([...(propValue as number[][])[0]]).flat();
      }
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
      // Check if it's a single tuple in a list that needs replication
      if ((propValue as number[][]).length === 1 && count > 1) {
        // Single tuple in a list [(x,y,z)] - replicate for all instances
        return Array(count).fill([...(propValue as number[][])[0]]).flat();
      }
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
 * Backend always sends [[x,y,z], ...] format after normalization.
 *
 * @param positionProp - The position prop value
 * @param fetchedPosition - Fetched position data (if applicable)
 * @returns Number of instances
 */
export function getInstanceCount(
  positionProp: string | number[][] | Record<string, any>,
  fetchedPosition: any
): number {
  if (fetchedPosition) {
    return fetchedPosition.length / 3;
  }

  if (Array.isArray(positionProp)) {
    return positionProp.length;  // Backend always sends [[x,y,z], ...]
  }

  return 0;  // Dynamic reference not yet fetched (string key or Transform)
}

/**
 * Get color count to detect single-color mode for UI
 *
 * @param colorProp - The color prop value
 * @returns Number of colors, or 0 if dynamic reference
 */
export function getColorCount(
  colorProp: string | string[]
): number {
  if (typeof colorProp === 'string') {
    return 0;  // Dynamic reference
  }
  return colorProp.length;  // List of hex strings
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
  count: number
): string[] {
  return colorHexArray.length === 1 && count > 1
    ? Array(count).fill(colorHexArray[0])
    : colorHexArray;
}
