/**
 * ScreenshotProvider - Provides screenshot capture functionality to the Three.js scene.
 *
 * This component:
 * 1. Registers a screenshot capture function in the Zustand store
 * 2. Listens for Socket.IO screenshot:request events for programmatic capture
 * 3. Handles both standard rendering and path-tracing modes
 *
 * Key design insight: By calling gl.render() immediately before canvas.toBlob(),
 * we ensure the WebGL buffer has content WITHOUT needing preserveDrawingBuffer.
 *
 * Path-tracing mode: When path-tracing is enabled, the canvas already shows
 * the path-traced result (progressive rendering to fullscreen quad).
 * We skip the gl.render() call to preserve the accumulated samples.
 */
import { useCallback, useEffect } from "react";
import { useThree } from "@react-three/fiber";
import { useAppStore } from "../../store";
import { socket } from "../../socket";

export function ScreenshotProvider() {
	const { gl, scene, camera } = useThree();
	const setScreenshotCapture = useAppStore(
		(state) => state.setScreenshotCapture,
	);
	const roomId = useAppStore((state) => state.roomId);
	const settings = useAppStore((state) => state.settings);

	// Check if path tracing is enabled
	const isPathTracingEnabled = settings?.pathtracing?.enabled === true;

	/**
	 * Capture screenshot from the WebGL canvas.
	 *
	 * Standard mode: Forces a render before capture to ensure the buffer
	 * has content, eliminating the need for preserveDrawingBuffer.
	 *
	 * Path-tracing mode: Captures the canvas as-is since it already contains
	 * the accumulated path-traced result.
	 */
	const captureScreenshot = useCallback(async (): Promise<Blob> => {
		// Only force render in standard mode
		// In path-tracing mode, the canvas already has the accumulated result
		if (!isPathTracingEnabled) {
			gl.render(scene, camera);
		}

		// Convert canvas to blob
		return new Promise((resolve, reject) => {
			gl.domElement.toBlob(
				(blob) => {
					if (blob) {
						resolve(blob);
					} else {
						reject(new Error("Failed to create blob from canvas"));
					}
				},
				"image/png",
				1.0,
			);
		});
	}, [gl, scene, camera, isPathTracingEnabled]);

	// Register the capture function in the store
	useEffect(() => {
		setScreenshotCapture(captureScreenshot);

		return () => {
			setScreenshotCapture(null);
		};
	}, [captureScreenshot, setScreenshotCapture]);

	// Listen for programmatic screenshot requests from Socket.IO
	useEffect(() => {
		const handleScreenshotRequest = async ({
			requestId,
			uploadUrl,
		}: {
			requestId: string;
			uploadUrl: string;
		}) => {
			if (!roomId) {
				console.error(
					"[ScreenshotProvider] Cannot handle screenshot request: no roomId",
				);
				return;
			}

			try {
				const blob = await captureScreenshot();

				// Upload via existing endpoint with request_id for correlation
				const formData = new FormData();
				formData.append(
					"file",
					blob,
					`screenshot_${requestId}.png`,
				);
				formData.append("format", "png");
				formData.append("width", gl.domElement.width.toString());
				formData.append("height", gl.domElement.height.toString());
				formData.append("request_id", requestId);

				const response = await fetch(uploadUrl, {
					method: "POST",
					body: formData,
				});

				if (!response.ok) {
					const error = await response
						.json()
						.catch(() => ({ error: "Upload failed" }));
					console.error(
						"[ScreenshotProvider] Screenshot upload failed:",
						error,
					);
				}
			} catch (error) {
				console.error(
					"[ScreenshotProvider] Error handling screenshot request:",
					error,
				);
			}
		};

		socket.on("screenshot:request", handleScreenshotRequest);

		return () => {
			socket.off("screenshot:request", handleScreenshotRequest);
		};
	}, [captureScreenshot, roomId, gl]);

	// This component doesn't render anything
	return null;
}
