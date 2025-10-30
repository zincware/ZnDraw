/**
 * Transform processor for filtering geometry data based on frame data.
 *
 * Transforms allow dynamic filtering of geometry properties (positions, radii, colors)
 * based on indices extracted from frame data (e.g., constraint indices).
 */

/**
 * InArrayTransform filters array data to only include elements at specified indices.
 *
 * Example: Filter positions to only show atoms with FixAtoms constraints
 * {
 *   type: "in_array",
 *   source: "constraints",           // Frame data key containing indices
 *   path: "0.kwargs.indices",        // Dot-separated path to extract indices
 *   filter: "arrays.positions"       // Frame data key to filter
 * }
 */
export interface InArrayTransform {
  type: "in_array";
  source: string;
  path: string;
  filter: string;
}

export type Transform = InArrayTransform;

/**
 * Check if a value is a Transform object.
 */
export function isTransform(value: unknown): value is Transform {
  return (
    typeof value === "object" &&
    value !== null &&
    "type" in value &&
    typeof (value as any).type === "string" &&
    (value as any).type === "in_array"
  );
}

/**
 * Extract a value from nested object using dot-separated path.
 *
 * @param obj - Object to extract from
 * @param path - Dot-separated path (e.g., "0.kwargs.indices" or "FixAtoms.indices")
 * @returns Extracted value or undefined if path is invalid
 *
 * @example
 * extractPath({ constraints: [{ kwargs: { indices: [0, 2] } }] }, "0.kwargs.indices")
 * // Returns: [0, 2]
 */
export function extractPath(obj: any, path: string): any {
  const parts = path.split(".");
  let current = obj;

  for (const part of parts) {
    if (current === null || current === undefined) {
      return undefined;
    }

    // Handle array indices (numeric strings)
    if (/^\d+$/.test(part)) {
      const index = parseInt(part, 10);
      if (Array.isArray(current) && index >= 0 && index < current.length) {
        current = current[index];
      } else {
        return undefined;
      }
    }
    // Handle object keys
    else if (typeof current === "object" && part in current) {
      current = current[part];
    } else {
      return undefined;
    }
  }

  return current;
}

/**
 * Evaluate a transform to filter array data based on extracted indices.
 *
 * @param transform - Transform definition
 * @param frameData - Frame data containing source and filter arrays
 * @returns Filtered Float32Array or null if transform cannot be evaluated
 *
 * @example
 * const frameData = {
 *   constraints: [{ kwargs: { indices: [0, 2, 4] } }],
 *   "arrays.positions": new Float32Array([...])  // 15 values (5 atoms * 3)
 * };
 *
 * const transform = {
 *   type: "in_array",
 *   source: "constraints",
 *   path: "0.kwargs.indices",
 *   filter: "arrays.positions"
 * };
 *
 * const filtered = evaluateTransform(transform, frameData);
 * // Returns: Float32Array with 9 values (3 atoms * 3) for indices [0, 2, 4]
 */
export function evaluateTransform(
  transform: Transform,
  frameData: Record<string, any>
): Float32Array | null {
  if (transform.type !== "in_array") {
    console.error(`[transformProcessor] Unknown transform type: ${(transform as any).type}`);
    return null;
  }

  // 1. Extract indices from frameData using source + path
  const sourceData = frameData[transform.source];
  if (!sourceData) {
    if (process.env.NODE_ENV !== "production") {
      console.debug(
        `[transformProcessor] Source key "${transform.source}" not found in frame data. Transform will return empty array.`
      );
    }
    return new Float32Array(0);
  }

  const indices = extractPath(sourceData, transform.path);
  if (!Array.isArray(indices)) {
    if (process.env.NODE_ENV !== "production") {
      console.debug(
        `[transformProcessor] Path "${transform.path}" did not extract array from source "${transform.source}". Transform will return empty array.`
      );
    }
    return new Float32Array(0);
  }

  if (indices.length === 0) {
    return new Float32Array(0);
  }

  // 2. Get the data to filter
  const filterData = frameData[transform.filter];
  if (!filterData) {
    if (process.env.NODE_ENV !== "production") {
      console.debug(
        `[transformProcessor] Filter key "${transform.filter}" not found in frame data. Transform will return empty array.`
      );
    }
    return new Float32Array(0);
  }

  // Convert to Float32Array if needed
  let sourceArray: Float32Array;
  if (filterData instanceof Float32Array) {
    sourceArray = filterData;
  } else if (filterData instanceof Float64Array || filterData instanceof Int32Array ||
             filterData instanceof Uint32Array || filterData instanceof Int16Array ||
             filterData instanceof Uint16Array || filterData instanceof Int8Array ||
             filterData instanceof Uint8Array) {
    // Handle other typed arrays - convert to Float32Array
    sourceArray = new Float32Array(filterData);
  } else if (Array.isArray(filterData)) {
    // Flatten nested arrays (e.g., [[x,y,z], [x,y,z]] -> [x,y,z,x,y,z])
    const flatArray = filterData.flat(Infinity) as number[];
    sourceArray = new Float32Array(flatArray);
  } else {
    console.error(
      `[transformProcessor] Filter data for "${transform.filter}" is not an array or Float32Array`
    );
    return null;
  }

  // 3. Determine stride (values per element)
  // For positions: 3 values per atom (x, y, z)
  // For radii: 1 value per atom
  // For colors (when stored as floats): 3-4 values per atom (r, g, b, a?)
  let stride: number;
  if (transform.filter.includes("position")) {
    stride = 3;
  } else if (transform.filter.includes("radii") || transform.filter.includes("radius")) {
    stride = 1;
  } else if (transform.filter.includes("color")) {
    // Colors in zndraw are hex strings, but if stored as floats, assume RGB
    stride = 3;
  } else {
    // Default to stride of 1 for unknown types
    stride = 1;
  }

  const elementCount = sourceArray.length / stride;

  // 4. Filter the array by indices
  const filteredValues: number[] = [];
  for (const index of indices) {
    if (typeof index !== "number" || index < 0 || index >= elementCount) {
      console.warn(
        `[transformProcessor] Invalid index ${index} (valid range: 0-${elementCount - 1}). Skipping.`
      );
      continue;
    }

    // Extract the values for this index
    const startIdx = index * stride;
    for (let i = 0; i < stride; i++) {
      filteredValues.push(sourceArray[startIdx + i]);
    }
  }

  return new Float32Array(filteredValues);
}

/**
 * Get all frame data keys required to evaluate a transform.
 * Used to fetch necessary data from the server.
 *
 * @param transform - Transform definition
 * @returns Array of frame data keys needed
 *
 * @example
 * getTransformSources({
 *   type: "in_array",
 *   source: "constraints",
 *   path: "0.kwargs.indices",
 *   filter: "arrays.positions"
 * })
 * // Returns: ["constraints", "arrays.positions"]
 */
export function getTransformSources(transform: Transform): string[] {
  if (transform.type === "in_array") {
    return [transform.source, transform.filter];
  }
  return [];
}
