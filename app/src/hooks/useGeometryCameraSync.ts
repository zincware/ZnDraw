/**
 * Hook for syncing OrbitControls changes back to Camera geometry.
 *
 * When attached to a Camera geometry with XYZ (direct coordinate) position/target,
 * user interactions via OrbitControls update the geometry data.
 *
 * With session cameras, there's always an attached camera - either the user's
 * session camera (cam:session:<id>) or a camera they explicitly switched to.
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
	const sessionId = useAppStore((state) => state.sessionId);
	const geometries = useAppStore((state) => state.geometries);
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	const lock = useAppStore((state) => state.lock);

	// Flag to prevent echo-back when geometry update triggers camera update
	const isUpdatingGeometryRef = useRef(false);

	// Derive the effective camera key - use attached or fall back to session camera
	const sessionCameraKey = sessionId ? `cam:session:${sessionId}` : null;
	const effectiveCameraKey = attachedCameraKey || sessionCameraKey;

	// Ref for latest persist function (avoids stale closure)
	const persistRef = useRef<
		((position: number[], target: number[]) => Promise<void>) | null
	>(null);

	// Update persist function ref
	persistRef.current = async (position: number[], target: number[]) => {
		if (!roomId || !effectiveCameraKey) return;

		const cameraGeometry = geometries[effectiveCameraKey];
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

		try {
			isUpdatingGeometryRef.current = true;
			await createGeometry(
				roomId,
				effectiveCameraKey,
				"Camera",
				updatedData,
				lock?.token,
			);
		} catch (error) {
			console.error(
				"[GeometryCameraSync] Failed to update camera geometry:",
				error,
			);
		} finally {
			// Reset flag after a short delay to allow the update to propagate
			setTimeout(() => {
				isUpdatingGeometryRef.current = false;
			}, 50);
		}
	};

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
		if (!effectiveCameraKey) return;

		// Skip if nothing is editable
		if (!controlsState.positionEditable && !controlsState.targetEditable)
			return;

		// Skip if we're in the middle of updating geometry (prevent loops)
		if (isUpdatingGeometryRef.current) return;

		const position = camera.position.toArray();
		const target = controls.target.toArray();

		debouncedPersist(position, target);
	}, [
		camera,
		controlsRef,
		effectiveCameraKey,
		controlsState.positionEditable,
		controlsState.targetEditable,
		debouncedPersist,
	]);

	return { syncToGeometry, isUpdatingGeometry: isUpdatingGeometryRef };
}
