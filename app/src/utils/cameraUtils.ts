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

// Note: interpolateCurvePosition and resolvePosition were removed as unused.
// Camera position resolution is handled directly in Camera.tsx using THREE.js curves.
