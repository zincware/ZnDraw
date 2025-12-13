import { useEffect, useRef, useMemo } from "react";
import * as THREE from "three";
import { useAppStore } from "../store";
import {
	getRelativePositions,
	applyTransformToPositions,
	isPositionStatic,
	isRotationStatic,
	isScaleStatic,
	getInitialRotations,
	getInitialScales,
	applyTransformToRotations,
	applyTransformToScales,
} from "../utils/geometryEditing";

import type { TypedArray } from "../myapi/client";
import type { Transform } from "../utils/transformProcessor";

// Type definitions for geometry data
// Accept TypedArray from server, number[][] from local state, or Transform
type PositionData = TypedArray | number[][] | Transform | null | undefined;
type RotationData = TypedArray | number[][] | null | undefined;
type ScaleData =
	| TypedArray
	| string
	| number
	| number[]
	| number[][]
	| [number, number, number]
	| [number, number, number][]
	| null
	| undefined;

// Type for the geometry data object passed to updateGeometry
interface GeometryUpdateData {
	position?: number[][];
	rotation?: number[][];
	scale?: number | [number, number, number][];
	[key: string]: unknown; // Allow other properties from fullGeometryData
}

/**
 * Hook for geometry components to handle editing mode with transform controls
 *
 * Responsibilities:
 * - Subscribe to transform matrix changes from MultiGeometryTransformControls
 * - Calculate relative positions/rotations/scales from initial centroid
 * - Apply transforms to selected instances and update local Zustand state
 *
 * NOTE: This hook does NOT persist to server - that's the component's responsibility!
 * Components should watch their data and debounce persistence themselves.
 *
 * @param geometryKey - Unique key for this geometry
 * @param positionData - Position data (must be static number[][] for editing)
 * @param rotationData - Rotation data (must be static number[][] for editing)
 * @param scaleData - Scale data (can be number, number[], or [number, number, number][] for editing)
 * @param selectedIndices - Array of selected instance indices
 * @param geometryType - Type of geometry (e.g., "Arrow", "Box")
 * @param fullGeometryData - Complete geometry data object (any type to allow flexibility)
 * @param instanceCount - Number of instances (for scale normalization)
 */
export function useGeometryEditing(
	geometryKey: string,
	positionData: PositionData,
	rotationData: RotationData,
	scaleData: ScaleData,
	selectedIndices: number[],
	geometryType: string,
	// eslint-disable-next-line @typescript-eslint/no-explicit-any
	fullGeometryData: any,
	instanceCount: number,
) {
	const {
		mode,
		transformMode,
		subscribeToEditing,
		updateGeometry,
		editingCombinedCentroid,
		geometries,
	} = useAppStore();

	const isEditing = mode === "editing";

	// Check which properties are static and editable
	const isPositionStaticVal = positionData && isPositionStatic(positionData);
	const isRotationStaticVal = rotationData && isRotationStatic(rotationData);
	const isScaleStaticVal = scaleData !== undefined && isScaleStatic(scaleData);

	// For editing mode, position must be static
	const canEdit = isPositionStaticVal;

	// Memoize selectedIndices to prevent infinite loops from array reference changes
	const stableSelectedIndices = useMemo(
		() => selectedIndices,
		[JSON.stringify(selectedIndices)],
	);

	// Use refs instead of state to avoid triggering re-renders and infinite loops
	const initialCentroidRef = useRef<[number, number, number] | null>(null);
	const relativePositionsRef = useRef<Map<number, THREE.Vector3>>(new Map());
	const initialRotationsRef = useRef<Map<number, THREE.Vector3>>(new Map());
	const initialScalesRef = useRef<Map<number, THREE.Vector3>>(new Map());
	const hasSelection = stableSelectedIndices.length > 0;

	// Store the last combined centroid we used to detect when it changes
	const lastCombinedCentroidRef = useRef<string | null>(null);

	// COMBINED effect: Initialize AND subscribe in one effect
	// This ensures subscription happens immediately after initialization
	useEffect(() => {
		// Skip entirely if not editable
		if (!canEdit) {
			return;
		}

		// Clear state when not in valid editing state
		if (!isEditing || !hasSelection || !editingCombinedCentroid) {
			initialCentroidRef.current = null;
			relativePositionsRef.current = new Map();
			initialRotationsRef.current = new Map();
			initialScalesRef.current = new Map();
			lastCombinedCentroidRef.current = null;
			return;
		}

		// Check if combined centroid has changed
		const centroidKey = editingCombinedCentroid.join(",");
		const centroidChanged = lastCombinedCentroidRef.current !== centroidKey;

		// Initialize or re-initialize if centroid changed
		if (centroidChanged) {
			// Access geometries directly from store to get current data
			const currentGeometry = geometries[geometryKey];
			const currentPositions = currentGeometry?.data?.position;
			const currentRotations = currentGeometry?.data?.rotation;
			const currentScales = currentGeometry?.data?.scale;

			// Use current data from store if available (may have been transformed)
			const positionsToUse =
				currentPositions && isPositionStatic(currentPositions)
					? currentPositions
					: positionData;

			if (!positionsToUse) {
				return;
			}

			const centroid = editingCombinedCentroid;
			initialCentroidRef.current = centroid;

			// Calculate relative positions using the combined centroid
			const relative = getRelativePositions(
				positionsToUse,
				stableSelectedIndices,
				centroid,
			);
			relativePositionsRef.current = relative;

			// Initialize rotations if available and static
			const rotationsToUse =
				currentRotations && isRotationStatic(currentRotations)
					? currentRotations
					: rotationData;
			if (rotationsToUse && isRotationStatic(rotationsToUse)) {
				initialRotationsRef.current = getInitialRotations(
					rotationsToUse,
					stableSelectedIndices,
				);
			}

			// Initialize scales if available and static
			const scalesToUse = currentScales ?? scaleData;
			if (scalesToUse !== undefined && isScaleStatic(scalesToUse)) {
				initialScalesRef.current = getInitialScales(
					scalesToUse,
					stableSelectedIndices,
					instanceCount,
				);
			}

			// Store this centroid to detect future changes
			lastCombinedCentroidRef.current = centroidKey;
		}

		// Now subscribe (refs are guaranteed to be set at this point)
		if (
			!initialCentroidRef.current ||
			relativePositionsRef.current.size === 0
		) {
			return;
		}

		if (!positionData) {
			return;
		}

		// Subscribe to transform changes from MultiGeometryTransformControls
		const unsubscribe = subscribeToEditing(
			geometryKey,
			(matrix: THREE.Matrix4) => {
				// Get current refs
				const currentCentroid = initialCentroidRef.current;
				const currentRelativePositions = relativePositionsRef.current;
				const currentInitialRotations = initialRotationsRef.current;
				const currentInitialScales = initialScalesRef.current;

				if (!currentCentroid || currentRelativePositions.size === 0) {
					return;
				}

				// Decompose the matrix
				const position = new THREE.Vector3();
				const quaternion = new THREE.Quaternion();
				const scale = new THREE.Vector3();
				matrix.decompose(position, quaternion, scale);

				// Build updated data based on current transform mode
				const updatedData: GeometryUpdateData = {
					...fullGeometryData,
				};

				// Get current mode from store (need fresh value in callback)
				const currentTransformMode = useAppStore.getState().transformMode;

				switch (currentTransformMode) {
					case "translate": {
						// Apply position transform
						const newPositions = applyTransformToPositions(
							positionData,
							stableSelectedIndices,
							matrix,
							currentCentroid,
							currentRelativePositions,
						);
						updatedData.position = newPositions;
						break;
					}
					case "rotate": {
						// Rotation mode: orbital rotation around centroid
						// 1. Update positions (objects orbit around centroid)
						const newPositions = applyTransformToPositions(
							positionData,
							stableSelectedIndices,
							matrix,
							currentCentroid,
							currentRelativePositions,
						);
						updatedData.position = newPositions;

						// 2. Update individual rotations if rotation data is static
						if (
							isRotationStaticVal &&
							rotationData &&
							currentInitialRotations.size > 0
						) {
							const newRotations = applyTransformToRotations(
								rotationData,
								stableSelectedIndices,
								quaternion,
								currentInitialRotations,
							);
							updatedData.rotation = newRotations;
						}
						break;
					}
					case "scale": {
						if (
							isScaleStaticVal &&
							scaleData !== undefined &&
							currentInitialScales.size > 0
						) {
							const newScales = applyTransformToScales(
								scaleData,
								stableSelectedIndices,
								scale,
								currentInitialScales,
								instanceCount,
							);
							updatedData.scale = newScales;
						}
						break;
					}
				}

				// Update local Zustand state immediately for responsive feedback
				updateGeometry(
					geometryKey,
					{
						type: geometryType,
						data: updatedData,
					},
					"local",
				);
			},
		);

		return () => {
			unsubscribe();
		};
	}, [
		isEditing,
		hasSelection,
		positionData,
		rotationData,
		scaleData,
		stableSelectedIndices,
		subscribeToEditing,
		updateGeometry,
		geometryKey,
		geometryType,
		fullGeometryData,
		canEdit,
		isRotationStaticVal,
		isScaleStaticVal,
		instanceCount,
		editingCombinedCentroid,
		geometries,
	]);
}
