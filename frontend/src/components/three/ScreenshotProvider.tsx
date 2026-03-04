import { useThree } from "@react-three/fiber";
/**
 * Central screenshot orchestrator. Registers capture function and handles
 * programmatic screenshot requests via Socket.IO.
 */
import { useCallback, useEffect } from "react";
import { completeScreenshot } from "../../myapi/client";
import { socket } from "../../socket";
import { useAppStore } from "../../store";

/**
 * Wait until all active geometries finish fetching and one rAF completes.
 * Mirrors the idle-detection pattern in useFrameLoadTime.ts.
 */
function waitForSceneIdle(timeout = 10000): Promise<void> {
	return new Promise((resolve, reject) => {
		const timer = setTimeout(
			() => reject(new Error("Scene idle timeout")),
			timeout,
		);
		const check = () => {
			const { geometryFetchingStates, geometries } = useAppStore.getState();
			const anyFetching = Object.entries(geometryFetchingStates).some(
				([key, isFetching]) => {
					if (!isFetching) return false;
					return geometries[key]?.data?.active !== false;
				},
			);
			if (!anyFetching) {
				clearTimeout(timer);
				unsub();
				requestAnimationFrame(() => resolve());
			}
		};
		const unsub = useAppStore.subscribe(() => check());
		check();
	});
}

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
			screenshot_id,
			upload_url,
		}: {
			screenshot_id: number;
			upload_url: string;
		}) => {
			if (!roomId) {
				console.error("[ScreenshotProvider] No roomId");
				return;
			}

			try {
				await waitForSceneIdle();
				const capture = getCapture();
				const blob = await capture();

				await completeScreenshot(
					upload_url,
					blob,
					"png",
					gl.domElement.width,
					gl.domElement.height,
				);
			} catch (error) {
				console.error("[ScreenshotProvider] Error:", error);
			}
		};

		socket.on("screenshot_request", handleScreenshotRequest);
		return () => {
			socket.off("screenshot_request", handleScreenshotRequest);
		};
	}, [getCapture, roomId, gl]);

	return null;
}
