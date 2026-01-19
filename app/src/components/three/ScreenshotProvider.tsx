/**
 * Central screenshot orchestrator. Registers capture function and handles
 * programmatic screenshot requests via Socket.IO.
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

	/** Standard mode capture - forces render before capture (no preserveDrawingBuffer needed). */
	const captureStandard = useCallback(async (): Promise<Blob> => {
		gl.render(scene, camera);
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

	const getCapture = useCallback((): (() => Promise<Blob>) => {
		return pathtracerCapture ?? captureStandard;
	}, [pathtracerCapture, captureStandard]);

	useEffect(() => {
		const capture = getCapture();
		setScreenshotCapture(capture);
		return () => setScreenshotCapture(null);
	}, [getCapture, setScreenshotCapture]);

	useEffect(() => {
		const handleScreenshotRequest = async ({
			requestId,
			uploadUrl,
		}: {
			requestId: string;
			uploadUrl: string;
		}) => {
			if (!roomId) {
				console.error("[ScreenshotProvider] No roomId");
				return;
			}

			try {
				const capture = getCapture();
				const blob = await capture();

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
					console.error("[ScreenshotProvider] Upload failed:", error);
				}
			} catch (error) {
				console.error("[ScreenshotProvider] Error:", error);
			}
		};

		socket.on("screenshot:request", handleScreenshotRequest);
		return () => socket.off("screenshot:request", handleScreenshotRequest);
	}, [getCapture, roomId, gl]);

	return null;
}
