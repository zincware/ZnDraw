import { useEffect, useMemo, useRef } from "react";
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
	const geometries = useAppStore((state) => state.geometries);
	const geometryUpdateSources = useAppStore(
		(state) => state.geometryUpdateSources,
	);

	// Ref to track if persistence is needed (dirty flag)
	const isDirtyRef = useRef(false);

	// Ref to hold the latest persist function (avoids stale closure in debounce)
	const persistRef = useRef<(() => Promise<void>) | undefined>(undefined);

	// Update persistRef with latest closure
	persistRef.current = async () => {
		if (!roomId) return;

		// Get fresh data from store
		const currentGeometries = useAppStore.getState().geometries;
		const currentGeometry = currentGeometries[geometryKey];
		if (!currentGeometry || !currentGeometry.data) return;

		if (!isDirtyRef.current) return;

		try {
			// No lock token needed - server handles via @check_lock
			await createGeometry(
				roomId,
				geometryKey,
				geometryType,
				currentGeometry.data,
			);
			isDirtyRef.current = false;
		} catch (error: unknown) {
			console.error(
				`[${geometryType}] Failed to persist ${geometryKey}:`,
				error,
			);
		}
	};

	// Stable debounced function that calls through ref
	const debouncedPersist = useMemo(
		() =>
			debounce(() => {
				persistRef.current?.();
			}, debounceMs),
		[debounceMs],
	);

	// Cleanup on unmount
	useEffect(() => {
		return () => {
			debouncedPersist.flush();
		};
	}, [debouncedPersist]);

	// Watch geometry data changes and persist - only if source is 'local'
	useEffect(() => {
		const currentGeometry = geometries[geometryKey];
		if (!currentGeometry) return;

		const updateSource = geometryUpdateSources[geometryKey];
		if (updateSource !== "local") {
			return;
		}

		isDirtyRef.current = true;
		debouncedPersist();
	}, [
		geometries[geometryKey],
		geometryUpdateSources[geometryKey],
		debouncedPersist,
		geometryKey,
	]);
}
