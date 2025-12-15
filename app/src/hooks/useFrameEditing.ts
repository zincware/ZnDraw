import { useEffect, useRef } from "react";
import { useQueryClient } from "@tanstack/react-query";
import * as THREE from "three";
import { useAppStore } from "../store";
import {
	getRelativePositions,
	applyTransformToPositions,
	convertFloat32ToPositionArray,
	convertPositionArrayToFloat32,
} from "../utils/geometryEditing";

/**
 * Compare two number arrays for equality
 */
function arraysEqual(a: number[], b: number[]): boolean {
	if (a.length !== b.length) return false;
	for (let i = 0; i < a.length; i++) {
		if (a[i] !== b[i]) return false;
	}
	return true;
}

/**
 * Hook for geometry components to handle editing mode with dynamic (frame) positions.
 *
 * Similar to useGeometryEditing but works with frame data fetched from the server.
 * Updates are applied optimistically to the TanStack Query cache for instant feedback,
 * then persisted to the server via the store's saveFrameEdits action.
 *
 * Key differences from useGeometryEditing:
 * - Works with Float32Array positions from frame data (not static number[][] from geometry)
 * - Updates the TanStack Query cache for instant dependent geometry updates (bonds, etc.)
 * - Registers pending edits with store for batch persistence
 *
 * @param geometryKey - Unique key for this geometry
 * @param positionKey - Frame data key for positions (e.g., "arrays.positions")
 * @param positionData - Current position data (Float32Array or null if not loaded)
 * @param selectedIndices - Array of selected instance indices
 * @param instanceCount - Total number of instances
 */
export function useFrameEditing(
	geometryKey: string,
	positionKey: string,
	positionData: Float32Array | null,
	selectedIndices: number[],
	instanceCount: number,
) {
	const queryClient = useQueryClient();
	const {
		mode,
		transformMode,
		currentFrame,
		roomId,
		subscribeToEditing,
		editingCombinedCentroid,
		setPendingFrameEdit,
		incrementEditingFrameDataCount,
		decrementEditingFrameDataCount,
	} = useAppStore();

	const isEditing = mode === "editing";

	// Track previous selectedIndices to detect changes without JSON.stringify
	const prevSelectedIndicesRef = useRef<number[]>(selectedIndices);
	const stableSelectedIndicesRef = useRef<number[]>(selectedIndices);

	// Update stable reference only when content actually changes
	if (!arraysEqual(prevSelectedIndicesRef.current, selectedIndices)) {
		prevSelectedIndicesRef.current = selectedIndices;
		stableSelectedIndicesRef.current = selectedIndices;
	}

	const stableSelectedIndices = stableSelectedIndicesRef.current;

	// Use refs instead of state to avoid triggering re-renders
	const initialCentroidRef = useRef<[number, number, number] | null>(null);
	const relativePositionsRef = useRef<Map<number, THREE.Vector3>>(new Map());
	const initialPositionsRef = useRef<Float32Array | null>(null);

	// Track the last combined centroid we used to detect when it changes
	const lastCombinedCentroidRef = useRef<string | null>(null);

	const hasSelection = stableSelectedIndices.length > 0;

	// Combined effect: Initialize AND subscribe
	useEffect(() => {
		// Skip if not editing or no position data
		if (!positionData) {
			return;
		}

		// Clear state when not in valid editing state
		if (!isEditing || !hasSelection || !editingCombinedCentroid) {
			initialCentroidRef.current = null;
			relativePositionsRef.current = new Map();
			initialPositionsRef.current = null;
			lastCombinedCentroidRef.current = null;
			return;
		}

		// Increment frame editing count (for reference counting with multiple hooks)
		incrementEditingFrameDataCount();

		// Check if combined centroid has changed
		const centroidKey = editingCombinedCentroid.join(",");
		const centroidChanged = lastCombinedCentroidRef.current !== centroidKey;

		// Initialize or re-initialize if centroid changed
		if (centroidChanged) {
			const centroid = editingCombinedCentroid;
			initialCentroidRef.current = centroid;

			// Store initial positions for transform calculations
			initialPositionsRef.current = new Float32Array(positionData);

			// Convert to number[][] for utility functions
			const positionsArray = convertFloat32ToPositionArray(positionData);

			// Calculate relative positions using the combined centroid
			const relative = getRelativePositions(
				positionsArray,
				stableSelectedIndices,
				centroid,
			);
			relativePositionsRef.current = relative;

			// Store this centroid to detect future changes
			lastCombinedCentroidRef.current = centroidKey;
		}

		// Now subscribe (refs are guaranteed to be set at this point)
		if (
			!initialCentroidRef.current ||
			relativePositionsRef.current.size === 0 ||
			!initialPositionsRef.current
		) {
			return;
		}

		// Subscribe to transform changes from MultiGeometryTransformControls
		const unsubscribe = subscribeToEditing(
			geometryKey,
			(matrix: THREE.Matrix4) => {
				// Handle translate and rotate modes (positions only)
				// Rotation is orbital rotation around centroid - positions change, not object rotation
				// (particles are spheres, so they don't have individual rotation)
				// Scale mode is not supported for dynamic positions
				const currentTransformMode = useAppStore.getState().transformMode;
				if (
					currentTransformMode !== "translate" &&
					currentTransformMode !== "rotate"
				) {
					return;
				}

				// Get current refs
				const currentCentroid = initialCentroidRef.current;
				const currentRelativePositions = relativePositionsRef.current;
				const currentInitialPositions = initialPositionsRef.current;

				if (
					!currentCentroid ||
					currentRelativePositions.size === 0 ||
					!currentInitialPositions
				) {
					return;
				}

				// Convert initial positions to number[][] format
				const initialPositionsArray = convertFloat32ToPositionArray(
					currentInitialPositions,
				);

				// Apply transform to get new positions
				const newPositionsArray = applyTransformToPositions(
					initialPositionsArray,
					stableSelectedIndices,
					matrix,
					currentCentroid,
					currentRelativePositions,
				);

				// Convert back to Float32Array
				const newPositions = convertPositionArrayToFloat32(newPositionsArray);

				// Update TanStack Query cache optimistically
				// This immediately updates all components using this data (including bonds)
				const currentRoomId = useAppStore.getState().roomId;
				const frame = useAppStore.getState().currentFrame;

				if (currentRoomId) {
					// Update the specific key in the frame response cache
					const queryKey = ["frame", currentRoomId, frame, positionKey];
					queryClient.setQueryData(queryKey, newPositions);

					// Also update any cache that might have the full frame response
					// This ensures consistency across different query patterns
					queryClient.setQueriesData(
						{ queryKey: ["frame", currentRoomId, frame] },
						(oldData: any) => {
							if (oldData && typeof oldData === "object") {
								return {
									...oldData,
									[positionKey]: newPositions,
								};
							}
							return oldData;
						},
					);

					// Register this as a pending edit for batch persistence
					setPendingFrameEdit(frame, positionKey, newPositions);
				}
			},
		);

		return () => {
			unsubscribe();
			// Decrement frame editing count when this hook stops editing
			decrementEditingFrameDataCount();
		};
	}, [
		isEditing,
		hasSelection,
		positionData,
		stableSelectedIndices,
		subscribeToEditing,
		geometryKey,
		positionKey,
		editingCombinedCentroid,
		queryClient,
		roomId,
		currentFrame,
		setPendingFrameEdit,
		incrementEditingFrameDataCount,
		decrementEditingFrameDataCount,
		transformMode,
	]);

	// Return whether this hook is actively handling frame editing
	return {
		isFrameEditing: isEditing && hasSelection && positionData !== null,
	};
}
