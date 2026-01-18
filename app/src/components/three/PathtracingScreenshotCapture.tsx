/**
 * PathtracingCaptureProvider - Registers pathtracer capture capability with the store.
 *
 * This component MUST be rendered inside <Pathtracer> to signal pathtracing mode.
 * It registers a capture function that reads from the canvas (NOT the render target).
 *
 * Why canvas instead of render target?
 * - The pathtracer's internal render target uses HDR float textures
 * - When renderToCanvas=true (default), the pathtracer copies tone-mapped result to canvas
 * - The canvas has the correctly processed, displayable image
 * - We enable preserveDrawingBuffer in Canvas.tsx when pathtracing is active
 *
 * All other screenshot logic (Socket.IO, upload) is handled by ScreenshotProvider (DRY).
 */
import { useCallback, useEffect } from "react";
import { useThree } from "@react-three/fiber";
import { useAppStore } from "../../store";

export function PathtracingCaptureProvider() {
	const { gl } = useThree();
	const setPathtracerCapture = useAppStore(
		(state) => state.setPathtracerCapture,
	);

	/**
	 * Capture from canvas.
	 *
	 * The pathtracer renders progressively to canvas with tone mapping applied.
	 * With preserveDrawingBuffer=true (set in Canvas.tsx for pathtracing mode),
	 * we can reliably capture the accumulated result.
	 */
	const captureFromCanvas = useCallback(async (): Promise<Blob> => {
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
	}, [gl]);

	// Register/unregister the capture function
	useEffect(() => {
		setPathtracerCapture(captureFromCanvas);
		return () => setPathtracerCapture(null);
	}, [captureFromCanvas, setPathtracerCapture]);

	return null;
}
