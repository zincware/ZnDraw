/**
 * Hook for syncing OrbitControls changes back to Camera geometry via REST.
 *
 * When attached to a Camera geometry with XYZ (direct coordinate) position/target,
 * user interactions via OrbitControls update the geometry data.
 */

import { debounce } from "lodash";
import { useCallback, useEffect, useMemo, useRef } from "react";
import type { Camera as ThreeCamera } from "three";
import type { OrbitControls as OrbitControlsImpl } from "three-stdlib";
import { createGeometry } from "../myapi/client";
import { useAppStore } from "../store";
import { isCurveAttachment } from "../utils/cameraUtils";
import type { ControlsState } from "./useCameraControls";

const DEBOUNCE_MS = 100; // Debounce geometry updates

interface UseGeometryCameraSyncProps {
	camera: ThreeCamera | null;
	controlsRef: React.RefObject<OrbitControlsImpl | null>;
	controlsState: ControlsState;
}

/**
 * Syncs OrbitControls changes to the attached Camera geometry via REST.
 */
export function useGeometryCameraSync({
	camera,
	controlsRef,
	controlsState,
}: UseGeometryCameraSyncProps) {
	const roomId = useAppStore((state) => state.roomId);
	const userId = useAppStore((state) => state.user?.id ?? null);
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	// Select only the attached camera geometry to avoid recreating callbacks
	// when unrelated geometries (other users' cameras) change.
	const attachedCameraGeometry = useAppStore((state) =>
		state.attachedCameraKey ? state.geometries[state.attachedCameraKey] : null,
	);

	// Ref for latest persist function (avoids stale closure)
	const persistRef = useRef<
		((position: number[], target: number[]) => Promise<void>) | null
	>(null);

	// Update persist function ref
	persistRef.current = async (position: number[], target: number[]) => {
		if (!roomId || !attachedCameraKey) return;

		if (!attachedCameraGeometry || attachedCameraGeometry.type !== "Camera")
			return;

		// Build updated data, preserving CurveAttachments
		const currentData = attachedCameraGeometry.data || {};
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
			await createGeometry(roomId, attachedCameraKey, "Camera", updatedData);
		} catch (error) {
			console.error(
				"[GeometryCameraSync] Failed to update camera geometry:",
				error,
			);
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
	 */
	const syncToGeometry = useCallback(() => {
		const controls = controlsRef.current;
		if (!camera || !controls) return;
		if (!attachedCameraKey) return;

		// Don't persist when following another user's camera (read-only)
		if (
			attachedCameraGeometry?.data?.owner &&
			attachedCameraGeometry.data.owner !== userId
		) {
			return;
		}

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
		userId,
		attachedCameraGeometry,
		controlsState.positionEditable,
		controlsState.targetEditable,
		debouncedPersist,
	]);

	return { syncToGeometry };
}
