/**
 * ScreenshotProvider - Central screenshot orchestrator.
 *
 * This component:
 * 1. Registers standard-mode capture function (when pathtracing is disabled)
 * 2. Listens for programmatic screenshot requests from Socket.IO
 * 3. Routes capture to either standard or pathtracer mode based on pathtracerCapture availability
 *
 * Key design insight: By calling gl.render() immediately before canvas.toBlob(),
 * we ensure the WebGL buffer has content WITHOUT needing preserveDrawingBuffer.
 *
 * DRY design: PathtracingCaptureProvider only registers its capture function to the store.
 * This component handles ALL common logic (Socket.IO, upload) and routes to the correct
 * capture function based on what's available.
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
	const pathtracerCapture = useAppStore((state) => state.pathtracerCapture);
	const roomId = useAppStore((state) => state.roomId);

	/**
	 * Capture screenshot from the WebGL canvas (standard mode only).
	 *
	 * Forces a render before capture to ensure the buffer has content,
	 * eliminating the need for preserveDrawingBuffer.
	 */
	const captureStandard = useCallback(async (): Promise<Blob> => {
		// Force render to ensure WebGL buffer has current frame
		gl.render(scene, camera);

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
	}, [gl, scene, camera]);

	/**
	 * Get the appropriate capture function based on what's available.
	 * Pathtracer takes priority when registered.
	 */
	const getCapture = useCallback((): (() => Promise<Blob>) => {
		return pathtracerCapture ?? captureStandard;
	}, [pathtracerCapture, captureStandard]);

	// Register the main capture function in the store
	// This always points to the best available capture method
	useEffect(() => {
		const capture = getCapture();
		setScreenshotCapture(capture);

		return () => {
			setScreenshotCapture(null);
		};
	}, [getCapture, setScreenshotCapture]);

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
				// Use the best available capture function
				const capture = getCapture();
				const blob = await capture();

				// Upload via existing endpoint with request_id for correlation
				const formData = new FormData();
				formData.append("file", blob, `screenshot_${requestId}.png`);
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
	}, [getCapture, roomId, gl]);

	// This component doesn't render anything
	return null;
}
