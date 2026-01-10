import { useEffect, useState } from "react";
import { isCurveAttachment } from "../utils/cameraUtils";

export interface ControlsState {
	enabled: boolean;
	enablePan: boolean;
	enableRotate: boolean;
	enableZoom: boolean;
	/** Whether position can be changed (XYZ mode) */
	positionEditable: boolean;
	/** Whether target can be changed (XYZ mode) */
	targetEditable: boolean;
}

/**
 * Hook to determine camera control states based on camera attachments.
 *
 * @param attachedCameraKey - The key of the currently attached camera geometry
 * @param geometries - All geometries in the scene
 * @returns Control state configuration for OrbitControls
 */
export function useCameraControls(
	attachedCameraKey: string | null,
	geometries: Record<string, any>,
): ControlsState {
	const [controlsState, setControlsState] = useState<ControlsState>({
		enabled: true,
		enablePan: true,
		enableRotate: true,
		enableZoom: true,
		positionEditable: true,
		targetEditable: true,
	});

	useEffect(() => {
		if (!attachedCameraKey) {
			// Not attached - full controls (default interactive camera)
			setControlsState({
				enabled: true,
				enablePan: true,
				enableRotate: true,
				enableZoom: true,
				positionEditable: true,
				targetEditable: true,
			});
			return;
		}

		const camera = geometries[attachedCameraKey];
		if (!camera || camera.type !== "Camera") {
			// Invalid attachment, fallback to full controls
			setControlsState({
				enabled: true,
				enablePan: true,
				enableRotate: true,
				enableZoom: true,
				positionEditable: true,
				targetEditable: true,
			});
			return;
		}

		// Check if position/target are CurveAttachment (locked) or XYZ (editable)
		const positionLocked = isCurveAttachment(camera.data.position);
		const targetLocked = isCurveAttachment(camera.data.target);

		if (positionLocked && targetLocked) {
			// Both locked - fully constrained cinematic mode
			setControlsState({
				enabled: false,
				enablePan: false,
				enableRotate: false,
				enableZoom: false,
				positionEditable: false,
				targetEditable: false,
			});
		} else if (positionLocked && !targetLocked) {
			// Position locked, target free
			// All OrbitControls ops change position, so disable all
			setControlsState({
				enabled: false,
				enablePan: false,
				enableRotate: false,
				enableZoom: false,
				positionEditable: false,
				targetEditable: true,
			});
		} else if (!positionLocked && targetLocked) {
			// Position free, target locked
			// Rotate and zoom only change position (OK), pan changes target (NOT OK)
			setControlsState({
				enabled: true,
				enablePan: false,
				enableRotate: true,
				enableZoom: true,
				positionEditable: true,
				targetEditable: false,
			});
		} else {
			// Both free - full controls, sync changes to geometry
			setControlsState({
				enabled: true,
				enablePan: true,
				enableRotate: true,
				enableZoom: true,
				positionEditable: true,
				targetEditable: true,
			});
		}
	}, [attachedCameraKey, geometries]);

	return controlsState;
}
