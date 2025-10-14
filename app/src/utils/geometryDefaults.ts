import { merge } from "lodash";

/**
 * Utility functions for handling geometry defaults from Pydantic models.
 * Follows Single Responsibility Principle - Pydantic is the single source of truth for defaults.
 */

/**
 * Merges geometry data with defaults from Pydantic models.
 *
 * This ensures that all required fields have values, even if they weren't sent by the API.
 * Pydantic models are the single source of truth for defaults.
 *
 * @param data - The geometry data from the server (may be partial)
 * @param type - The geometry type (e.g., "Arrow", "Box", "Sphere")
 * @param defaults - The geometry defaults object from Zustand store
 * @returns Complete geometry data with all defaults applied
 *
 * @example
 * ```ts
 * const arrowData = getGeometryWithDefaults(
 *   { position: "arrays.positions", direction: "calc.forces" },
 *   "Arrow",
 *   geometryDefaults
 * );
 * // Result will include opacity, selecting, hovering from Arrow defaults
 * ```
 */
export function getGeometryWithDefaults<T extends Record<string, any>>(
  data: Partial<T>,
  type: string,
  defaults: Record<string, any>
): T {
  // If defaults aren't loaded yet (during app initialization), use data as-is
  if (!defaults || Object.keys(defaults).length === 0) {
    return data as T;
  }

  const defaultForType = defaults[type];

  // If this specific geometry type has no defaults, use data as-is
  if (!defaultForType) {
    return data as T;
  }

  // Deep merge: defaults first, then override with actual data
  // lodash merge mutates the first argument, so pass empty object
  return merge({}, defaultForType, data);
}
