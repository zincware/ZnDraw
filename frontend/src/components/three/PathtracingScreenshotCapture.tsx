import { useThree } from "@react-three/fiber";
/**
 * Registers pathtracer capture function. Must be rendered inside <Pathtracer>.
 * Captures from canvas (requires preserveDrawingBuffer=true in Canvas.tsx).
 */
import { useCallback, useEffect } from "react";
import { useAppStore } from "../../store";

export function PathtracingCaptureProvider() {
	const { gl } = useThree();
	const setPathtracerCapture = useAppStore(
		(state) => state.setPathtracerCapture,
	);

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

	useEffect(() => {
		setPathtracerCapture(captureFromCanvas);
		return () => setPathtracerCapture(null);
	}, [captureFromCanvas, setPathtracerCapture]);

	return null;
}
