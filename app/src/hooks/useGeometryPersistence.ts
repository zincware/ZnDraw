import { useCallback, useEffect, useMemo, useRef } from "react";
import { debounce } from "lodash";
import { useAppStore } from "../store";
import { createGeometry } from "../myapi/client";

/**
 * Hook for managing geometry persistence to the server
 *
 * Handles debounced persistence with proper cleanup and race condition prevention.
 * Only persists when the update source is "local" (from user edits, not server updates).
 *
 * @param geometryKey - Unique key for this geometry
 * @param geometryType - Type of geometry (e.g., "Arrow", "Box", "Sphere")
 * @param debounceMs - Debounce delay in milliseconds (default: 500)
 */
export function useGeometryPersistence(
	geometryKey: string,
	geometryType: string,
	debounceMs: number = 500,
) {
	const roomId = useAppStore((state) => state.roomId);
	const lock = useAppStore((state) => state.lock);
	const geometries = useAppStore((state) => state.geometries);
	const geometryUpdateSources = useAppStore(
		(state) => state.geometryUpdateSources,
	);

	// Ref to track if persistence is needed (dirty flag)
	const isDirtyRef = useRef(false);

	// Persistence callback - persists geometry changes to server
	const persistGeometryData = useCallback(async () => {
		if (!roomId) return;

		const currentGeometry = geometries[geometryKey];
		if (!currentGeometry || !currentGeometry.data) return;

		// Only persist if data is dirty
		if (!isDirtyRef.current) return;

		try {
			await createGeometry(
				roomId,
				geometryKey,
				geometryType,
				currentGeometry.data,
				lock?.token,
			);
			// Clear dirty flag after successful persistence
			isDirtyRef.current = false;
		} catch (error: unknown) {
			console.error(
				`[${geometryType}] Failed to persist ${geometryKey}:`,
				error,
			);
			// Keep dirty flag set on error so we can retry
		}
	}, [roomId, geometryKey, geometryType, geometries, lock]);

	// Create stable debounced function
	const debouncedPersist = useMemo(
		() => debounce(persistGeometryData, debounceMs),
		[persistGeometryData, debounceMs],
	);

	// Cleanup debounce on unmount and flush any pending changes
	useEffect(() => {
		return () => {
			// Cancel any pending debounced calls
			debouncedPersist.cancel();
			// Flush immediately if there are unsaved changes
			if (isDirtyRef.current) {
				persistGeometryData();
			}
		};
	}, [debouncedPersist, persistGeometryData]);

	// Watch geometry data changes and persist - only if source is 'local'
	useEffect(() => {
		const currentGeometry = geometries[geometryKey];
		if (!currentGeometry) return;

		// Only persist if update source is 'local' (not from server)
		const updateSource = geometryUpdateSources[geometryKey];
		if (updateSource !== "local") {
			return;
		}

		// Mark as dirty and trigger debounced persist
		isDirtyRef.current = true;
		debouncedPersist();
	}, [
		// Watch the geometry itself and the update source
		geometries[geometryKey],
		geometryUpdateSources[geometryKey],
		debouncedPersist,
		geometryKey,
	]);
}
