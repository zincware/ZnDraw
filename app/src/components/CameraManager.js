/**
 * CameraManager - Updates camera properties from session camera geometry.
 *
 * Uses session camera data to configure the Three.js camera's near/far planes,
 * position, target, and other properties.
 *
 * Position/target sync only happens when viewing through the session camera
 * (not when attached to another camera geometry).
 */

import { useEffect, useRef } from "react";
import { useThree } from "@react-three/fiber";
import { useAppStore } from "../store";

function CameraManager({ sessionCameraData }) {
	const { camera, controls } = useThree();
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	const sessionId = useAppStore((state) => state.sessionId);

	// Derive session camera key
	const sessionCameraKey = sessionId ? `cam:session:${sessionId}` : null;

	// Track if this is the initial mount to avoid position sync on first render
	// (let the Canvas's initial camera position take precedence)
	const isInitialMount = useRef(true);

	useEffect(() => {
		if (!sessionCameraData) return;

		camera.near = sessionCameraData.near;
		camera.far = sessionCameraData.far;

		if (camera.isPerspectiveCamera && sessionCameraData.fov) {
			camera.fov = sessionCameraData.fov;
		}

		// CRITICAL: Re-calculate projection matrix with new values
		camera.updateProjectionMatrix();

		// Only sync position/target when viewing through the session camera
		// (not when attached to another camera) and not on initial mount
		const isViewingSessionCamera = attachedCameraKey === sessionCameraKey;

		if (isViewingSessionCamera && !isInitialMount.current) {
			const position = sessionCameraData.position;
			const target = sessionCameraData.target;

			// Sync position if it's a direct coordinate array
			if (Array.isArray(position) && position.length === 3) {
				camera.position.set(position[0], position[1], position[2]);
			}

			// Sync target if controls exist and target is a direct coordinate array
			if (controls?.target && Array.isArray(target) && target.length === 3) {
				controls.target.set(target[0], target[1], target[2]);
				controls.update();
			}
		}

		isInitialMount.current = false;
	}, [
		sessionCameraData,
		camera,
		controls,
		attachedCameraKey,
		sessionCameraKey,
	]);

	return null;
}

export default CameraManager;
