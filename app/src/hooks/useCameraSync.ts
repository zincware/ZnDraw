/**
 * Camera sync hook for synchronizing OrbitControls state with backend.
 *
 * Sends camera state updates when OrbitControls changes (throttled, 10Hz max).
 * Receives programmatic camera control from Python clients.
 *
 * Key features:
 * - Lazy sync: Only emits when state actually changes
 * - Throttled: Max 10 updates per second to prevent spam
 * - No echo-back: Ignores updates that came from programmatic control
 */

import { useCallback, useEffect, useRef, type RefObject } from "react";
import type {
	Camera as ThreeCamera,
	PerspectiveCamera,
	OrthographicCamera,
} from "three";
import type { OrbitControls as OrbitControlsImpl } from "three-stdlib";
import { socket } from "../socket";
import { useAppStore } from "../store";

// Type for cameras with projection properties
type ProjectionCamera = PerspectiveCamera | OrthographicCamera;

const THROTTLE_MS = 100; // 10Hz max update rate

interface CameraControlPayload {
	sessionId: string;
	camera: {
		position: [number, number, number];
		target: [number, number, number];
		fov?: number;
		zoom?: number;
	};
}

/**
 * Hook to synchronize camera state between frontend OrbitControls and backend.
 *
 * @param camera - Three.js camera
 * @param controlsRef - Ref to OrbitControls instance (ref because it may be null initially)
 * @returns Object with syncCamera function to call on controls change
 */
export function useCameraSync(
	camera: ThreeCamera | null,
	controlsRef: RefObject<OrbitControlsImpl | null>,
) {
	const sessionId = useAppStore((state) => state.sessionId);
	const lastState = useRef<string>("");
	const lastSync = useRef<number>(0);
	const ignoreNextUpdate = useRef<boolean>(false);

	/**
	 * Send camera state to backend (called when OrbitControls updates).
	 * Throttled and deduplicated to prevent spam.
	 */
	const syncCamera = useCallback(() => {
		const controls = controlsRef.current;
		if (!camera || !controls || !sessionId) return;

		// Skip if this update came from programmatic control
		if (ignoreNextUpdate.current) {
			ignoreNextUpdate.current = false;
			return;
		}

		const now = Date.now();
		if (now - lastSync.current < THROTTLE_MS) return;

		// Build state key for deduplication
		const state = JSON.stringify({
			position: camera.position.toArray(),
			target: controls.target.toArray(),
		});

		// Only emit if state actually changed (lazy sync)
		if (state === lastState.current) return;

		lastState.current = state;
		lastSync.current = now;

		// Emit camera state update to backend
		// Backend stores in session_cameras hash, does NOT broadcast
		const projCamera = camera as ProjectionCamera;
		socket.emit("camera:state:update", {
			sessionId,
			position: camera.position.toArray(),
			target: controls.target.toArray(),
			fov: "fov" in camera ? (camera as PerspectiveCamera).fov : 75,
			near: projCamera.near,
			far: projCamera.far,
			zoom: projCamera.zoom,
		});
	}, [camera, controlsRef, sessionId]);

	/**
	 * Listen for programmatic camera control from Python clients.
	 * Always registers the listener - null checks are done inside the handler.
	 * This ensures the listener is registered even if camera/controls aren't ready yet.
	 */
	useEffect(() => {
		const handleCameraControl = (data: CameraControlPayload) => {
			// Get current values from refs/state - these may have changed since effect ran
			const controls = controlsRef.current;
			const currentSessionId = useAppStore.getState().sessionId;

			// Null checks inside handler, not in effect setup
			if (!camera || !controls) {
				console.warn(
					"[CameraSync] camera:control received but camera/controls not ready",
				);
				return;
			}

			// Only apply if this is for our session
			if (data.sessionId !== currentSessionId) {
				return;
			}

			// Set flag to prevent echo-back
			ignoreNextUpdate.current = true;

			// Apply camera state
			const cameraData = data.camera;
			if (cameraData.position) {
				camera.position.fromArray(cameraData.position);
			}
			if (cameraData.target) {
				controls.target.fromArray(cameraData.target);
			}

			const projCamera = camera as ProjectionCamera;
			if (cameraData.fov && "fov" in camera) {
				(camera as PerspectiveCamera).fov = cameraData.fov;
				projCamera.updateProjectionMatrix();
			}
			if (cameraData.zoom !== undefined) {
				projCamera.zoom = cameraData.zoom;
				projCamera.updateProjectionMatrix();
			}

			// Update controls to reflect new camera state
			controls.update();
		};

		socket.on("camera:control", handleCameraControl);
		return () => {
			socket.off("camera:control", handleCameraControl);
		};
		// Note: camera is still a dependency because we need to re-register
		// if the Three.js camera instance changes. controlsRef doesn't change.
	}, [camera, controlsRef]);

	return { syncCamera };
}

/**
 * Register this session with the backend (called once on mount).
 *
 * @param sessionId - Session identifier
 * @param alias - Optional alias from URL parameter
 */
export function registerSession(sessionId: string, alias?: string | null) {
	socket.emit("session:register", {
		sessionId,
		alias: alias || null,
	});
}
