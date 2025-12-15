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
	type ScaleData as UtilScaleData,
} from "../utils/geometryEditing";

import type { TypedArray } from "../myapi/client";
import type { Transform } from "../utils/transformProcessor";
import type { ScaleProp as BackendScaleProp } from "../utils/geometryData";

// Type definitions for geometry data
// Accept TypedArray from server, number[][] from local state, or Transform
type PositionData = TypedArray | number[][] | Transform | null | undefined;
type RotationData = TypedArray | number[][] | null | undefined;
// ScaleData accepts backend prop (string | number[][]), TypedArray from server, or static number
type ScaleData = BackendScaleProp | TypedArray | number | null | undefined;

// Type for the geometry data object passed to updateGeometry
interface GeometryUpdateData {
	position?: number[][];
	rotation?: number[][];
	scale?: [number, number, number][]; // Always anisotropic after transform
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
		updateSelectionForGeometry,
		editingCombinedCentroid,
		geometries,
	} = useAppStore();

	const isEditing = mode === "editing";

	// Check which properties are static and editable
	const isPositionStaticVal = positionData && isPositionStatic(positionData);
	const staticPositions = isPositionStaticVal
		? (positionData as number[][])
		: undefined;

	const isRotationStaticVal = rotationData && isRotationStatic(rotationData);
	const staticRotations = isRotationStaticVal
		? (rotationData as number[][])
		: undefined;

	// ScaleData can be a number (fallback 1.0) or ScaleProp (string | number[][])
	const isScaleStaticVal = scaleData !== undefined && isScaleStatic(scaleData);
	// Cast for staticScales to what the utility functions expect
	const staticScales = isScaleStaticVal
		? (scaleData as UtilScaleData)
		: undefined;

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
		// Skip if position data was explicitly passed as null (dynamic positions handled elsewhere)
		// This allows useFrameEditing to handle dynamic positions without this hook deselecting
		if (positionData === null) {
			return;
		}

		// 1. Validation Logic: Deselect if dynamic attributes are targeted by transform mode
		// "As soon as a transform controls action is applied... If any of these are dynamic, they should be deselected"
		if (isEditing && hasSelection) {
			let shouldDeselect = false;

			// Position must always be static for any editing (when this hook handles positions)
			if (!staticPositions) {
				shouldDeselect = true;
			}
			// If rotating and geometry has rotation data, rotation must be static
			// (geometries without rotation data, like Curve, can still orbit positions)
			else if (
				transformMode === "rotate" &&
				rotationData !== null &&
				!staticRotations
			) {
				shouldDeselect = true;
			}
			// If scaling, scale must be static
			else if (transformMode === "scale" && staticScales === undefined) {
				shouldDeselect = true;
			}

			if (shouldDeselect) {
				// Clear selection for this geometry
				updateSelectionForGeometry(geometryKey, []);
				return;
			}
		}

		// Skip if not editable (this covers the general static position check)
		if (!staticPositions) {
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
					: staticPositions;

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
			let rotationsToUse =
				currentRotations && isRotationStatic(currentRotations)
					? currentRotations
					: staticRotations;

			if (rotationsToUse) {
				// Handle broadcasted rotation (single value for multiple instances)
				// We must expand it before getting initial rotations for specific indices
				if (rotationsToUse.length === 1 && instanceCount > 1) {
					const template = rotationsToUse[0];
					rotationsToUse = Array.from({ length: instanceCount }, () => [
						...template,
					]);
				}

				initialRotationsRef.current = getInitialRotations(
					rotationsToUse,
					stableSelectedIndices,
				);
			}

			// Initialize scales if available and static
			const scalesToUse = currentScales ?? staticScales;
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

		if (!staticPositions) {
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
							staticPositions,
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
							staticPositions,
							stableSelectedIndices,
							matrix,
							currentCentroid,
							currentRelativePositions,
						);
						updatedData.position = newPositions;

						// 2. Update individual rotations if rotation data is static
						if (staticRotations && currentInitialRotations.size > 0) {
							// Handle broadcasted rotation: expand if needed
							let rotationDataToUse = staticRotations;
							if (staticRotations.length === 1 && instanceCount > 1) {
								const template = staticRotations[0];
								rotationDataToUse = Array.from(
									{ length: instanceCount },
									() => [...template],
								);
							}

							const newRotations = applyTransformToRotations(
								rotationDataToUse,
								stableSelectedIndices,
								quaternion,
								currentInitialRotations,
							);
							updatedData.rotation = newRotations;
						}
						break;
					}
					case "scale": {
						if (staticScales !== undefined && currentInitialScales.size > 0) {
							const newScales = applyTransformToScales(
								staticScales,
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
		staticPositions,
		staticRotations,
		staticScales,
		stableSelectedIndices,
		subscribeToEditing,
		updateGeometry,
		updateSelectionForGeometry,
		geometryKey,
		geometryType,
		fullGeometryData,
		instanceCount,
		editingCombinedCentroid,
		geometries,
		transformMode,
	]);
}
