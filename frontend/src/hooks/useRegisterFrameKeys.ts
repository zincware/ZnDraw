import { useEffect } from "react";
import { useAppStore } from "../store";

/**
 * Registers frame data keys needed by a geometry for prefetching.
 *
 * Call from each geometry component with the list of frame-level keys
 * it will query (e.g. "arrays.positions", "arrays.colors"). The
 * prefetcher reads these registrations to batch-fetch upcoming frames.
 *
 * Automatically unregisters on unmount or when keys change.
 */
export function useRegisterFrameKeys(
	geometryKey: string,
	frameKeys: string[],
): void {
	const registerFrameKeys = useAppStore((s) => s.registerFrameKeys);
	const unregisterFrameKeys = useAppStore((s) => s.unregisterFrameKeys);

	useEffect(() => {
		registerFrameKeys(geometryKey, frameKeys);
		return () => unregisterFrameKeys(geometryKey);
	}, [geometryKey, frameKeys, registerFrameKeys, unregisterFrameKeys]);
}
