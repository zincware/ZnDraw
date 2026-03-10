import {
	type GeometryData,
	getAllBookmarks,
	getGeometry,
	listGeometries,
	listSelectionGroups,
} from "../../myapi/client";
import { useAppStore } from "../../store";
import type { HandlerContext } from "./types";
import { createInvalidateHandler } from "./utils";

// --- Typed event interfaces ---

export interface GeometryInvalidateEvent {
	operation?: "set" | "delete";
	key?: string;
}

export interface DefaultCameraInvalidateEvent {
	room_id: string;
	default_camera: string | null;
}

export interface ActiveCameraUpdateEvent {
	active_camera: string;
}

// --- Factory ---

export function createGeometryHandlers(ctx: HandlerContext) {
	async function onGeometriesInvalidate(data: GeometryInvalidateEvent) {
		if (!ctx.roomId) return;

		try {
			const operation = data?.operation || "set"; // default to 'set' for backward compatibility

			if (operation === "delete") {
				// Handle geometry deletion
				const { key } = data;
				if (!key) {
					console.warn("Delete operation received without key");
					return;
				}

				// Remove from the store
				ctx.removeGeometry(key);

				// Remove the specific geometry from the cache
				ctx.queryClient.removeQueries({
					queryKey: ["geometries", ctx.roomId, "detail", key],
				});

				// Invalidate the list to update the UI
				ctx.queryClient.invalidateQueries({
					queryKey: ["geometries", ctx.roomId, "list"],
				});
			} else if (operation === "set") {
				// Handle geometry creation/update
				if (data && data.key) {
					const { key } = data;

					// Invalidate the specific geometry detail query
					ctx.queryClient.invalidateQueries({
						queryKey: ["geometries", ctx.roomId, "detail", key],
					});

					// Invalidate the list to update the geometry grid
					ctx.queryClient.invalidateQueries({
						queryKey: ["geometries", ctx.roomId, "list"],
					});

					// Fetch only the updated geometry
					try {
						const response = await getGeometry(ctx.roomId, key);
						// Update only this specific geometry in the store
						ctx.updateGeometry(key, response.geometry);

						// Auto-select newly created curves ONLY if no curve is currently selected
						const currentActiveCurve =
							useAppStore.getState().activeCurveForDrawing;
						if (
							!currentActiveCurve &&
							response.geometry.type === "Curve" &&
							response.geometry.data?.active !== false
						) {
							ctx.setActiveCurveForDrawing(key);
						}
					} catch (error) {
						console.error(`Error fetching geometry ${key}:`, error);
					}
				} else {
					// No specific key - refetch all geometries (fallback for backward compatibility)

					// Invalidate React Query cache for geometries
					ctx.queryClient.invalidateQueries({
						queryKey: ["geometries", ctx.roomId, "list"],
					});

					// Fetch list of geometry keys first
					const listResponse = await listGeometries(ctx.roomId);
					const keys = Object.keys(listResponse.items || {});

					// Fetch all geometries in parallel
					const geometryPromises = keys.map(async (key: string) => {
						try {
							const response = await getGeometry(ctx.roomId!, key);
							return { key, geometry: response.geometry };
						} catch (error) {
							console.error(`Error fetching geometry ${key}:`, error);
							return null;
						}
					});

					const geometries = await Promise.all(geometryPromises);

					// Build geometries object from results
					const geometriesObj: Record<string, GeometryData> = {};
					for (const item of geometries) {
						if (item && item.geometry) {
							geometriesObj[item.key] = item.geometry;
						}
					}

					ctx.setGeometries(geometriesObj);
				}
			}
		} catch (error) {
			console.error("Error handling geometry invalidation:", error);
		}
	}

	const onSelectionsInvalidate = createInvalidateHandler(
		listGeometries,
		(response) => {
			const geos = response.items || {};
			const selectionsFromGeos: Record<string, number[]> = {};
			for (const [key, geo] of Object.entries(geos)) {
				if (geo.type === "Camera") continue;
				selectionsFromGeos[key] = geo.selection ?? [];
			}
			ctx.setSelections(selectionsFromGeos);
		},
		"selections",
		() => ctx.roomId,
	);

	const onSelectionGroupsInvalidate = createInvalidateHandler(
		listSelectionGroups,
		(response) => {
			ctx.setSelectionGroups(response.items || {});
		},
		"selection_groups",
		() => ctx.roomId,
	);

	const onBookmarksInvalidate = createInvalidateHandler(
		getAllBookmarks,
		(response) => {
			ctx.setBookmarks(response.items || {});
		},
		"bookmarks",
		() => ctx.roomId,
	);

	function onDefaultCameraInvalidate(data: DefaultCameraInvalidateEvent) {
		ctx.queryClient.setQueryData(["defaultCamera", data.room_id], {
			default_camera: data.default_camera,
		});
	}

	function onActiveCameraUpdate(data: ActiveCameraUpdateEvent) {
		const { setAttachedCameraKey } = useAppStore.getState();
		setAttachedCameraKey(data.active_camera);
	}

	return {
		onGeometriesInvalidate,
		onSelectionsInvalidate,
		onSelectionGroupsInvalidate,
		onBookmarksInvalidate,
		onDefaultCameraInvalidate,
		onActiveCameraUpdate,
	};
}
