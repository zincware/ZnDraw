/**
 * Hook for syncing OrbitControls changes back to Camera geometry.
 *
 * When attached to a Camera geometry with XYZ (direct coordinate) position/target,
 * user interactions via OrbitControls update the geometry data.
 */

import { useCallback, useEffect, useRef, useMemo } from "react";
import type { Camera as ThreeCamera } from "three";
import type { OrbitControls as OrbitControlsImpl } from "three-stdlib";
import { debounce } from "lodash";
import { useAppStore } from "../store";
import { createGeometry } from "../myapi/client";
import { isCurveAttachment } from "../utils/cameraUtils";
import type { ControlsState } from "./useCameraControls";

const DEBOUNCE_MS = 100; // Debounce geometry updates
const POSITION_EPSILON = 0.0001; // Tolerance for position comparison

/**
 * Check if two position arrays are approximately equal.
 */
function positionsEqual(
	a: number[] | undefined,
	b: number[] | undefined,
): boolean {
	if (!a || !b || a.length !== b.length) return false;
	return a.every((v, i) => Math.abs(v - b[i]) < POSITION_EPSILON);
}

interface UseGeometryCameraSyncProps {
	camera: ThreeCamera | null;
	controlsRef: React.RefObject<OrbitControlsImpl | null>;
	controlsState: ControlsState;
}

/**
 * Syncs OrbitControls changes to the attached Camera geometry.
 */
export function useGeometryCameraSync({
	camera,
	controlsRef,
	controlsState,
}: UseGeometryCameraSyncProps) {
	const roomId = useAppStore((state) => state.roomId);
	const geometries = useAppStore((state) => state.geometries);
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);

	// Track last sent values for echo-back detection (value-based, not time-based)
	const lastSentRef = useRef<{
		position: number[] | null;
		target: number[] | null;
	}>({ position: null, target: null });

	// Ref for latest persist function (avoids stale closure)
	const persistRef = useRef<
		((position: number[], target: number[]) => Promise<void>) | null
	>(null);

	// Update persist function ref
	persistRef.current = async (position: number[], target: number[]) => {
		if (!roomId || !attachedCameraKey) return;

		const cameraGeometry = geometries[attachedCameraKey];
		if (!cameraGeometry || cameraGeometry.type !== "Camera") return;

		// Build updated data, preserving CurveAttachments
		const currentData = cameraGeometry.data || {};
		const updatedData = { ...currentData };

		// Only update position if it's editable (XYZ mode)
		if (
			controlsState.positionEditable &&
			!isCurveAttachment(currentData.position)
		) {
			updatedData.position = position;
		}

		// Only update target if it's editable (XYZ mode)
		if (
			controlsState.targetEditable &&
			!isCurveAttachment(currentData.target)
		) {
			updatedData.target = target;
		}

		// Store what we're sending for echo-back detection
		lastSentRef.current = {
			position: updatedData.position as number[] | null,
			target: updatedData.target as number[] | null,
		};

		try {
			// No lock token needed - server handles via @check_lock
			await createGeometry(roomId, attachedCameraKey, "Camera", updatedData);
		} catch (error) {
			console.error(
				"[GeometryCameraSync] Failed to update camera geometry:",
				error,
			);
			// Clear last sent on error so next update will proceed
			lastSentRef.current = { position: null, target: null };
		}
	};

	/**
	 * Check if incoming geometry data matches what we last sent (echo-back detection).
	 * Returns true if this appears to be our own update echoed back.
	 */
	const isEchoBack = useCallback(
		(geomData: { position?: number[]; target?: number[] }): boolean => {
			const { position, target } = lastSentRef.current;
			if (!position && !target) return false;

			// Check if position matches (for array positions only)
			const positionMatches =
				!position ||
				!Array.isArray(geomData.position) ||
				positionsEqual(position, geomData.position);

			// Check if target matches (for array targets only)
			const targetMatches =
				!target ||
				!Array.isArray(geomData.target) ||
				positionsEqual(target, geomData.target);

			// If both match, it's likely our echo-back
			if (positionMatches && targetMatches) {
				// Clear the tracking after detecting echo-back
				lastSentRef.current = { position: null, target: null };
				return true;
			}

			return false;
		},
		[],
	);

	// Debounced persist
	const debouncedPersist = useMemo(
		() =>
			debounce((position: number[], target: number[]) => {
				persistRef.current?.(position, target);
			}, DEBOUNCE_MS),
		[],
	);

	// Cleanup on unmount
	useEffect(() => {
		return () => {
			debouncedPersist.flush();
		};
	}, [debouncedPersist]);

	/**
	 * Sync function to call when OrbitControls changes.
	 * Updates the camera geometry (session camera or explicitly attached camera).
	 */
	const syncToGeometry = useCallback(() => {
		const controls = controlsRef.current;
		if (!camera || !controls) return;
		if (!attachedCameraKey) return;

		// Skip if nothing is editable
		if (!controlsState.positionEditable && !controlsState.targetEditable)
			return;

		const position = camera.position.toArray();
		const target = controls.target.toArray();

		debouncedPersist(position, target);
	}, [
		camera,
		controlsRef,
		attachedCameraKey,
		controlsState.positionEditable,
		controlsState.targetEditable,
		debouncedPersist,
	]);

	return { syncToGeometry, isEchoBack };
}
