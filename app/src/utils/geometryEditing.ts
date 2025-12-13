import * as THREE from "three";

/**
 * Utility functions for geometry editing with transform controls
 */

// Type for static scale data from backend (always per-instance anisotropic)
export type ScaleData = [number, number, number][];

/**
 * Check if a position property is static (number[][]) rather than dynamic (string reference or TypedArray)
 * @param position - Position property from geometry data
 * @returns true if position is static and can be edited
 */
export function isPositionStatic(position: any): position is number[][] {
	return (
		Array.isArray(position) && position.length > 0 && Array.isArray(position[0])
	);
}

/**
 * Calculate the centroid (center point) of selected instances
 * @param positions - Array of all position tuples [[x, y, z], ...]
 * @param indices - Array of selected instance indices
 * @returns Centroid as [x, y, z] tuple
 */
export function calculateCentroid(
	positions: number[][],
	indices: number[],
): [number, number, number] {
	if (indices.length === 0 || positions.length === 0) {
		return [0, 0, 0];
	}

	let sumX = 0;
	let sumY = 0;
	let sumZ = 0;
	let count = 0;

	for (const idx of indices) {
		if (idx >= 0 && idx < positions.length) {
			const [x, y, z] = positions[idx];
			sumX += x;
			sumY += y;
			sumZ += z;
			count++;
		}
	}

	if (count === 0) {
		return [0, 0, 0];
	}

	return [sumX / count, sumY / count, sumZ / count];
}

/**
 * Get positions relative to a centroid
 * Used to maintain relative positions when transforming around centroid
 * @param positions - Array of all position tuples
 * @param indices - Array of selected instance indices
 * @param centroid - Centroid position [x, y, z]
 * @returns Map of index -> relative position vector
 */
export function getRelativePositions(
	positions: number[][],
	indices: number[],
	centroid: [number, number, number],
): Map<number, THREE.Vector3> {
	const relativePositions = new Map<number, THREE.Vector3>();

	for (const idx of indices) {
		if (idx >= 0 && idx < positions.length) {
			const [x, y, z] = positions[idx];
			const relative = new THREE.Vector3(
				x - centroid[0],
				y - centroid[1],
				z - centroid[2],
			);
			relativePositions.set(idx, relative);
		}
	}

	return relativePositions;
}

/**
 * Apply a transform matrix to selected positions
 * @param positions - Array of all position tuples (will be modified in place)
 * @param indices - Array of selected instance indices to transform
 * @param transformMatrix - THREE.Matrix4 representing the transformation
 * @param centroid - Original centroid position before transform
 * @param relativePositions - Map of relative positions from centroid
 * @returns New positions array with transformations applied
 */
export function applyTransformToPositions(
	positions: number[][],
	indices: number[],
	transformMatrix: THREE.Matrix4,
	initialCentroid: [number, number, number],
	relativePositions: Map<number, THREE.Vector3>,
): number[][] {
	// Create a copy to avoid mutating original
	const newPositions = positions.map((pos) => [...pos]) as number[][];

	// Decompose the transform matrix into position, rotation, and scale
	const position = new THREE.Vector3();
	const quaternion = new THREE.Quaternion();
	const scale = new THREE.Vector3();
	transformMatrix.decompose(position, quaternion, scale);

	// Calculate the initial centroid position
	const initialCentroidVec = new THREE.Vector3(
		initialCentroid[0],
		initialCentroid[1],
		initialCentroid[2],
	);

	// Calculate DELTA movement from initial position
	const deltaMovement = position.clone().sub(initialCentroidVec);

	// Create rotation/scale matrix from decomposed components
	const rotationScaleMatrix = new THREE.Matrix4();
	rotationScaleMatrix.compose(new THREE.Vector3(0, 0, 0), quaternion, scale);

	// For each selected instance
	for (const idx of indices) {
		if (idx >= 0 && idx < newPositions.length) {
			const relativePos = relativePositions.get(idx);
			if (relativePos) {
				// 1. Apply rotation/scale to relative position
				const transformedRelative = relativePos.clone();
				transformedRelative.applyMatrix4(rotationScaleMatrix);

				// 2. Add back to INITIAL centroid, then add movement delta
				const newPos = initialCentroidVec
					.clone()
					.add(transformedRelative)
					.add(deltaMovement);

				newPositions[idx] = [newPos.x, newPos.y, newPos.z];
			}
		}
	}

	return newPositions;
}

/**
 * Validate that all selected instances have static positions
 * Returns indices that need to be deselected (have dynamic positions)
 * @param geometries - Record of all geometries
 * @param selections - Record of selections per geometry
 * @returns Object with { validSelections, invalidGeometries }
 */
export function validateStaticPositions(
	geometries: Record<string, any>,
	selections: Record<string, number[]>,
): {
	validSelections: Record<string, number[]>;
	invalidGeometries: string[];
} {
	const validSelections: Record<string, number[]> = {};
	const invalidGeometries: string[] = [];

	for (const [geometryKey, selectedIndices] of Object.entries(selections)) {
		// Skip if no selection
		if (!selectedIndices || selectedIndices.length === 0) {
			continue;
		}

		const geometry = geometries[geometryKey];
		if (!geometry) {
			continue;
		}

		const position = geometry.data?.position;

		// Check if position is static
		if (isPositionStatic(position)) {
			validSelections[geometryKey] = selectedIndices;
		} else {
			invalidGeometries.push(geometryKey);
		}
	}

	return { validSelections, invalidGeometries };
}

/**
 * Calculate a combined centroid for all selected instances across multiple geometries
 * @param geometries - Record of all geometries
 * @param validSelections - Record of valid (static) selections per geometry
 * @returns Combined centroid position or null if no valid selections
 */
export function calculateCombinedCentroid(
	geometries: Record<string, any>,
	validSelections: Record<string, number[]>,
): [number, number, number] | null {
	let sumX = 0;
	let sumY = 0;
	let sumZ = 0;
	let totalCount = 0;

	for (const [geometryKey, selectedIndices] of Object.entries(
		validSelections,
	)) {
		if (selectedIndices.length === 0) {
			continue;
		}

		const geometry = geometries[geometryKey];
		if (!geometry) {
			continue;
		}

		const positions = geometry.data?.position;
		if (!isPositionStatic(positions)) {
			continue;
		}

		// Add up positions for this geometry
		for (const idx of selectedIndices) {
			if (idx >= 0 && idx < positions.length) {
				const [x, y, z] = positions[idx];
				sumX += x;
				sumY += y;
				sumZ += z;
				totalCount++;
			}
		}
	}

	if (totalCount === 0) {
		return null;
	}

	return [sumX / totalCount, sumY / totalCount, sumZ / totalCount];
}

/**
 * Get all relative positions for selected instances across multiple geometries
 * @param geometries - Record of all geometries
 * @param validSelections - Record of valid selections per geometry
 * @param centroid - Combined centroid position
 * @returns Map of geometryKey -> Map of index -> relative position
 */
export function getAllRelativePositions(
	geometries: Record<string, any>,
	validSelections: Record<string, number[]>,
	centroid: [number, number, number],
): Map<string, Map<number, THREE.Vector3>> {
	const allRelativePositions = new Map<string, Map<number, THREE.Vector3>>();

	for (const [geometryKey, selectedIndices] of Object.entries(
		validSelections,
	)) {
		const geometry = geometries[geometryKey];
		if (!geometry) {
			continue;
		}

		const positions = geometry.data?.position;
		if (!isPositionStatic(positions)) {
			continue;
		}

		const relativePositions = getRelativePositions(
			positions,
			selectedIndices,
			centroid,
		);

		allRelativePositions.set(geometryKey, relativePositions);
	}

	return allRelativePositions;
}

/**
 * Check if a rotation property is static (number[][]) rather than dynamic (string reference)
 * @param rotation - Rotation property from geometry data
 * @returns true if rotation is static and can be edited
 */
export function isRotationStatic(rotation: any): rotation is number[][] {
	return (
		Array.isArray(rotation) && rotation.length > 0 && Array.isArray(rotation[0])
	);
}

/**
 * Check if scale property is static (not a dynamic string reference)
 * @param scale - Scale property from geometry data
 * @returns true if scale is static [[x,y,z],...] and can be edited
 */
export function isScaleStatic(scale: any): scale is ScaleData {
	if (typeof scale === "string") return false; // Dynamic reference
	return Array.isArray(scale);
}

/**
 * Get the initial rotation values for selected indices
 * @param rotations - Array of all rotation tuples [[x, y, z], ...]
 * @param indices - Array of selected instance indices
 * @returns Map of index -> initial Euler rotation as Vector3
 */
export function getInitialRotations(
	rotations: number[][],
	indices: number[],
): Map<number, THREE.Vector3> {
	const initialRotations = new Map<number, THREE.Vector3>();

	for (const idx of indices) {
		if (idx >= 0 && idx < rotations.length) {
			const [rx, ry, rz] = rotations[idx];
			initialRotations.set(idx, new THREE.Vector3(rx, ry, rz));
		}
	}

	return initialRotations;
}

/**
 * Get the initial scale values for selected indices
 * Normalizes scale to per-instance anisotropic format for editing
 * @param scale - Scale data (various formats)
 * @param indices - Array of selected instance indices
 * @param instanceCount - Total number of instances
 * @returns Map of index -> initial scale as Vector3 [sx, sy, sz]
 */
export function getInitialScales(
	scale: ScaleData,
	indices: number[],
	instanceCount: number,
): Map<number, THREE.Vector3> {
	const initialScales = new Map<number, THREE.Vector3>();

	// Normalize scale to per-instance anisotropic format
	const scaleArray = normalizeScaleToArray(scale, instanceCount);

	for (const idx of indices) {
		if (idx >= 0 && idx < scaleArray.length) {
			const [sx, sy, sz] = scaleArray[idx];
			initialScales.set(idx, new THREE.Vector3(sx, sy, sz));
		}
	}

	return initialScales;
}

/**
 * Normalize scale data to per-instance anisotropic array format
 * @param scale - Scale data as [number, number, number][]
 * @param instanceCount - Number of instances
 * @returns Array of [sx, sy, sz] tuples with exactly instanceCount entries
 */
export function normalizeScaleToArray(
	scale: ScaleData,
	instanceCount: number,
): [number, number, number][] {
	if (scale.length === 0) {
		// Empty array: default to [1,1,1] for all instances
		return Array.from(
			{ length: instanceCount },
			() => [1, 1, 1] as [number, number, number],
		);
	}

	// Handle broadcasted scale: [[sx,sy,sz]] -> [[sx,sy,sz], [sx,sy,sz], ...]
	if (scale.length === 1) {
		const template = scale[0];
		return Array.from(
			{ length: instanceCount },
			() => [...template] as [number, number, number],
		);
	}

	// Clone each entry and ensure exactly instanceCount entries
	// Truncate if longer, pad with [1,1,1] if shorter
	const result: [number, number, number][] = scale
		.slice(0, instanceCount)
		.map((s) => [...s] as [number, number, number]);

	while (result.length < instanceCount) {
		result.push([1, 1, 1]);
	}

	return result;
}

/**
 * Apply a rotation delta (quaternion) to selected rotations
 * @param rotations - Array of all rotation tuples (Euler angles)
 * @param indices - Array of selected instance indices
 * @param deltaQuaternion - Rotation delta to apply
 * @param initialRotations - Map of initial rotation values
 * @returns New rotations array with transformations applied
 */
export function applyTransformToRotations(
	rotations: number[][],
	indices: number[],
	deltaQuaternion: THREE.Quaternion,
	initialRotations: Map<number, THREE.Vector3>,
): number[][] {
	// Create a copy to avoid mutating original
	const newRotations = rotations.map((rot) => [...rot]) as number[][];

	for (const idx of indices) {
		const initialRot = initialRotations.get(idx);
		if (initialRot && idx >= 0 && idx < newRotations.length) {
			// Convert initial Euler to quaternion
			const euler = new THREE.Euler(initialRot.x, initialRot.y, initialRot.z);
			const initialQuat = new THREE.Quaternion().setFromEuler(euler);

			// Apply delta rotation
			const newQuat = deltaQuaternion.clone().multiply(initialQuat);

			// Convert back to Euler
			const newEuler = new THREE.Euler().setFromQuaternion(newQuat);
			newRotations[idx] = [newEuler.x, newEuler.y, newEuler.z];
		}
	}

	return newRotations;
}

/**
 * Apply scale factors to selected instances' scales
 * @param scale - Current scale data
 * @param indices - Array of selected instance indices
 * @param scaleFactors - Scale factors to apply [sx, sy, sz]
 * @param initialScales - Map of initial scale values
 * @param instanceCount - Total number of instances
 * @returns New scale array with transformations applied (always anisotropic format)
 */
export function applyTransformToScales(
	scale: ScaleData,
	indices: number[],
	scaleFactors: THREE.Vector3,
	initialScales: Map<number, THREE.Vector3>,
	instanceCount: number,
): [number, number, number][] {
	// Normalize current scale to per-instance anisotropic format
	const newScales = normalizeScaleToArray(scale, instanceCount);

	for (const idx of indices) {
		const initialScale = initialScales.get(idx);
		if (initialScale && idx >= 0 && idx < newScales.length) {
			// Apply scale factors to initial scale
			newScales[idx] = [
				initialScale.x * scaleFactors.x,
				initialScale.y * scaleFactors.y,
				initialScale.z * scaleFactors.z,
			];
		}
	}

	return newScales;
}

/**
 * Apply uniform scale factor for geometries with a single scale value
 * Uses geometric mean of scale factors for uniform scaling
 * @param initialScale - Original uniform scale value
 * @param scaleFactors - Scale factors from transform [sx, sy, sz]
 * @returns New uniform scale value
 */
export function applyUniformScale(
	initialScale: number,
	scaleFactors: THREE.Vector3,
): number {
	// Use geometric mean for uniform scaling from potentially non-uniform transform
	const geometricMean = Math.cbrt(
		scaleFactors.x * scaleFactors.y * scaleFactors.z,
	);
	return initialScale * geometricMean;
}
