/**
 * Camera utilities for handling CurveAttachment and position resolution.
 *
 * Cameras can use either direct coordinates or CurveAttachment for position/target.
 * These utilities help resolve and type-check camera position types.
 */

/**
 * CurveAttachment references a Curve geometry with progress along the curve.
 * Similar to InArrayTransform, this allows referencing external geometry data.
 */
export interface CurveAttachment {
	type: "curve_attachment";
	geometry_key: string;
	progress: number;
}

/**
 * Position can be either direct coordinates or a CurveAttachment.
 */
export type PositionType = [number, number, number] | CurveAttachment;

/**
 * Check if a value is a CurveAttachment.
 *
 * @param value - The value to check
 * @returns True if value is a CurveAttachment
 */
export function isCurveAttachment(value: unknown): value is CurveAttachment {
	return (
		typeof value === "object" &&
		value !== null &&
		"type" in value &&
		(value as CurveAttachment).type === "curve_attachment"
	);
}

/**
 * Interpolate a position along a curve at a given progress.
 *
 * @param curvePositions - Array of [x, y, z] positions defining the curve
 * @param progress - Progress along the curve (0.0 to 1.0)
 * @returns Interpolated [x, y, z] position
 */
export function interpolateCurvePosition(
	curvePositions: [number, number, number][],
	progress: number,
): [number, number, number] {
	if (curvePositions.length === 0) {
		return [0, 0, 0];
	}
	if (curvePositions.length === 1) {
		return curvePositions[0];
	}

	// Clamp progress to [0, 1]
	const t = Math.max(0, Math.min(1, progress));

	// Map progress to segment index
	const segmentCount = curvePositions.length - 1;
	const segmentProgress = t * segmentCount;
	const segmentIndex = Math.min(Math.floor(segmentProgress), segmentCount - 1);
	const localT = segmentProgress - segmentIndex;

	// Linear interpolation between segment points
	const p0 = curvePositions[segmentIndex];
	const p1 = curvePositions[segmentIndex + 1];

	return [
		p0[0] + (p1[0] - p0[0]) * localT,
		p0[1] + (p1[1] - p0[1]) * localT,
		p0[2] + (p1[2] - p0[2]) * localT,
	];
}

/**
 * Resolve a PositionType to actual coordinates.
 *
 * If position is a CurveAttachment, looks up the referenced curve geometry
 * and interpolates the position at the specified progress.
 *
 * @param position - Direct coordinates or CurveAttachment
 * @param geometries - All geometries in the scene (for CurveAttachment lookup)
 * @returns Resolved [x, y, z] coordinates
 */
export function resolvePosition(
	position: PositionType,
	geometries: Record<string, any>,
): [number, number, number] {
	if (isCurveAttachment(position)) {
		const curve = geometries[position.geometry_key];
		if (!curve || curve.type !== "Curve") {
			console.warn(
				`CurveAttachment references missing/invalid curve: ${position.geometry_key}`,
			);
			return [0, 0, 0]; // Fallback to origin
		}

		const curvePositions = curve.data?.position as
			| [number, number, number][]
			| undefined;
		if (!curvePositions || curvePositions.length === 0) {
			return [0, 0, 0]; // Empty curve fallback
		}

		return interpolateCurvePosition(curvePositions, position.progress);
	}

	// Direct coordinates - ensure it's a proper tuple
	if (Array.isArray(position) && position.length >= 3) {
		return [position[0], position[1], position[2]];
	}

	console.warn("Invalid position type:", position);
	return [0, 0, 0];
}

/**
 * Session camera state as received from/sent to backend.
 */
export interface SessionCameraState {
	position: PositionType;
	target: PositionType;
	fov: number;
	near: number;
	far: number;
	zoom: number;
}
