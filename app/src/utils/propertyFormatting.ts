/**
 * Utility functions for formatting property values consistently across components.
 * Centralizes formatting logic to maintain DRY principle.
 */

import type { PropertyInfo } from "../types/property-inspector";

/**
 * Check if a value is array-like (includes regular arrays and TypedArrays).
 */
export const isArrayLike = (value: unknown): value is ArrayLike<unknown> => {
  return (
    (Array.isArray(value) || ArrayBuffer.isView(value)) &&
    (value as ArrayLike<unknown>).length > 0
  );
};

/**
 * Format a single value for display (numbers, bigints, etc.).
 */
export const formatValue = (value: unknown): string => {
  if (typeof value === "number") {
    return value.toFixed(3);
  }
  if (typeof value === "bigint") {
    return value.toString();
  }
  return String(value);
};

/**
 * Format an array value with optional truncation.
 */
export const formatArray = (
  values: ArrayLike<unknown>,
  maxItems: number = 5
): string => {
  const arrayValues = Array.from(values);

  if (arrayValues.length <= maxItems) {
    return `[${arrayValues.map(formatValue).join(", ")}]`;
  }

  return `[${arrayValues.slice(0, 3).map(formatValue).join(", ")}, ... +${arrayValues.length - 3}]`;
};

/**
 * Extract per-particle value from array-like data.
 * Handles both 1D arrays and multi-dimensional arrays based on property metadata.
 */
export const extractParticleValue = (
  value: ArrayLike<unknown>,
  particleId: number,
  propertyInfo?: PropertyInfo
): string => {
  const shape = propertyInfo?.metadata?.shape;

  // Check if multi-dimensional based on shape metadata
  if (shape && Array.isArray(shape) && shape.length > 1) {
    // Multi-dimensional: extract slice for this particle
    const numParticles = shape[0];
    const elementsPerParticle = shape[1];

    // Validate the particle ID is within the expected range
    if (particleId < 0 || particleId >= numParticles) {
      return "—";
    }

    const startIdx = particleId * elementsPerParticle;
    const endIdx = startIdx + elementsPerParticle;

    // Validate indices are within array bounds
    if (endIdx > value.length) {
      return "—";
    }

    const particleSlice = Array.from(
      (value as any).slice(startIdx, endIdx)
    );
    return `[${particleSlice.map(formatValue).join(", ")}]`;
  }

  // 1D array: direct indexing
  // Validate particleId is within bounds for 1D array
  if (particleId < 0 || particleId >= value.length) {
    return "—";
  }

  const particleValue = value[particleId];
  if (particleValue === undefined || particleValue === null) {
    return "—";
  }
  return formatValue(particleValue);
};

/**
 * Format property value for display.
 * Handles scalars, arrays, and per-particle extraction.
 */
export const formatPropertyValue = (
  value: unknown,
  particleId?: number,
  propertyInfo?: PropertyInfo
): string => {
  if (value === undefined || value === null) {
    return "—";
  }

  if (isArrayLike(value)) {
    if (particleId !== undefined && particleId !== null) {
      return extractParticleValue(value, particleId, propertyInfo);
    }
    return formatArray(value);
  }

  return formatValue(value);
};
