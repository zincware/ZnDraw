import { useEffect, useRef, useState, useMemo } from "react";
import * as THREE from "three";
import { useAppStore } from "../store";
import {
  getRelativePositions,
  applyTransformToPositions,
  isPositionStatic,
} from "../utils/geometryEditing";

/**
 * Hook for geometry components to handle editing mode with transform controls
 *
 * Responsibilities:
 * - Subscribe to transform matrix changes from MultiGeometryTransformControls
 * - Calculate relative positions from initial centroid
 * - Apply transforms to selected instances and update local Zustand state
 *
 * NOTE: This hook does NOT persist to server - that's the component's responsibility!
 * Components should watch their position data and debounce persistence themselves.
 *
 * @param geometryKey - Unique key for this geometry
 * @param positionData - Position data (must be static number[][])
 * @param selectedIndices - Array of selected instance indices
 * @param geometryType - Type of geometry (e.g., "Arrow", "Box")
 * @param fullGeometryData - Complete geometry data object
 */
export function useGeometryEditing(
  geometryKey: string,
  positionData: any, // Can be number[][], string, TypedArray, or undefined
  selectedIndices: number[],
  geometryType: string,
  fullGeometryData: any
) {
  const {
    mode,
    subscribeToEditing,
    updateGeometry,
    editingCombinedCentroid,
    geometries,
  } = useAppStore();

  const isEditing = mode === 'editing';

  // Early return if position is not static - this geometry can't be edited
  // This prevents any side effects for dynamic geometries like particles
  const isStatic = positionData && isPositionStatic(positionData);

  // Memoize selectedIndices to prevent infinite loops from array reference changes
  const stableSelectedIndices = useMemo(
    () => selectedIndices,
    [JSON.stringify(selectedIndices)]
  );

  // Use refs instead of state to avoid triggering re-renders and infinite loops
  const initialCentroidRef = useRef<[number, number, number] | null>(null);
  const relativePositionsRef = useRef<Map<number, THREE.Vector3>>(new Map());
  const hasSelection = stableSelectedIndices.length > 0;

  // Store the last combined centroid we used to detect when it changes
  const lastCombinedCentroidRef = useRef<string | null>(null);

  // Initialize centroid and relative positions when entering edit mode OR when combined centroid changes
  // IMPORTANT: Recalculates when combined centroid changes (selection changes) but NOT when positions change
  useEffect(() => {
    // Skip entirely if not static
    if (!isStatic) {
      return;
    }

    // Clear state when exiting edit mode
    if (!isEditing || !hasSelection || !editingCombinedCentroid) {
      initialCentroidRef.current = null;
      relativePositionsRef.current = new Map();
      lastCombinedCentroidRef.current = null;
      return;
    }

    // Check if combined centroid has changed
    const centroidKey = editingCombinedCentroid.join(',');
    if (lastCombinedCentroidRef.current === centroidKey) {
      // Same centroid, don't recalculate (prevents drift from position updates)
      return;
    }

    // Combined centroid changed - recalculate using CURRENT positions from store
    // This handles selection changes while in edit mode
    // Access geometries directly from store to avoid dependency issues
    const currentGeometry = geometries[geometryKey];
    const currentPositions = currentGeometry?.data?.position;

    // Use current positions from store if available (may have been transformed)
    // Otherwise fall back to original positionData
    const positionsToUse = (currentPositions && isPositionStatic(currentPositions))
      ? currentPositions
      : positionData;

    if (!positionsToUse) {
      return;
    }

    const centroid = editingCombinedCentroid;

    initialCentroidRef.current = centroid;

    // Calculate relative positions using the combined centroid from current positions
    const relative = getRelativePositions(positionsToUse, stableSelectedIndices, centroid);
    relativePositionsRef.current = relative;

    // Store this centroid to detect future changes
    lastCombinedCentroidRef.current = centroidKey;
  }, [isEditing, hasSelection, stableSelectedIndices, isStatic, editingCombinedCentroid, positionData, geometryKey]);

  // Subscribe to transform changes
  useEffect(() => {
    // Skip entirely if not static
    if (!isStatic) {
      return;
    }

    if (!isEditing || !hasSelection || !initialCentroidRef.current || relativePositionsRef.current.size === 0) {
      return;
    }

    if (!positionData) {
      return;
    }

    // Subscribe to transform changes from MultiGeometryTransformControls
    const unsubscribe = subscribeToEditing(geometryKey, (matrix: THREE.Matrix4) => {
      // Get current refs (they may have been updated by selection changes)
      const currentCentroid = initialCentroidRef.current;
      const currentRelativePositions = relativePositionsRef.current;

      if (!currentCentroid || currentRelativePositions.size === 0) {
        return;
      }

      // Apply transform to this geometry's selected instances
      const newPositions = applyTransformToPositions(
        positionData,
        stableSelectedIndices,
        matrix,
        currentCentroid,
        currentRelativePositions
      );

      // Update local Zustand state immediately for responsive feedback
      // Mark as 'local' so component knows to persist to server
      updateGeometry(geometryKey, {
        type: geometryType,
        data: {
          ...fullGeometryData,
          position: newPositions,
        },
      }, 'local');
    });

    return () => {
      unsubscribe();
    };
  }, [
    isEditing,
    hasSelection,
    positionData,
    stableSelectedIndices,
    subscribeToEditing,
    updateGeometry,
    geometryKey,
    geometryType,
    fullGeometryData,
    isStatic,
  ]);
}
