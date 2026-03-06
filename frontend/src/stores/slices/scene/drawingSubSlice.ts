import type * as THREE from "three";
import type { StateCreator } from "zustand";
import { createGeometry, getGeometry } from "../../../myapi/client";
import type { AppState } from "../../../store";
import { getActiveCurves, selectPreferredCurve } from "./geometrySubSlice";

export interface DrawingSubSlice {
	drawingPointerPosition: THREE.Vector3 | null;
	drawingIsValid: boolean;
	activeCurveForDrawing: string | null;

	enterDrawingMode: (queryClient?: any) => Promise<void>;
	exitDrawingMode: () => Promise<void>;
	setDrawingPointerPosition: (position: THREE.Vector3 | null) => void;
	setDrawingIsValid: (isValid: boolean) => void;
	setActiveCurveForDrawing: (key: string | null) => void;
}

export const createDrawingSubSlice: StateCreator<
	AppState,
	[],
	[],
	DrawingSubSlice
> = (set, get) => ({
	drawingPointerPosition: null,
	drawingIsValid: false,
	activeCurveForDrawing: null,

	enterDrawingMode: async (queryClient?: any) => {
		const state = get();
		const { geometries, mode, roomId } = state;

		if (mode !== "view") {
			console.warn(`Cannot enter drawing mode from ${mode} mode`);
			return;
		}

		const acquired = await get().acquireLock("drawing on curve");
		if (!acquired) return;

		const activeCurves = getActiveCurves(geometries);

		if (activeCurves.length === 0) {
			if (!roomId) {
				console.warn("Cannot create curve: no room ID");
				return;
			}
			set((state) => ({
				geometries: {
					...state.geometries,
					curve: { type: "Curve", data: { position: [] } },
				},
				activeCurveForDrawing: "curve",
				mode: "drawing",
			}));
			try {
				await createGeometry(roomId, "curve", "Curve", { position: [] });
				const response = await getGeometry(roomId, "curve");
				set((state) => ({
					geometries: { ...state.geometries, curve: response.geometry },
				}));
				if (queryClient && roomId) {
					queryClient.invalidateQueries({
						queryKey: ["geometries", roomId, "list"],
					});
				}
			} catch (error) {
				console.error("Error creating default curve:", error);
				set((state) => {
					const { curve: removed, ...rest } = state.geometries;
					return {
						geometries: rest,
						activeCurveForDrawing: null,
						mode: "view",
					};
				});
			}
			return;
		}

		if (activeCurves.length === 1) {
			set({ activeCurveForDrawing: activeCurves[0], mode: "drawing" });
			return;
		}

		const lastSelectedCurve = state.activeCurveForDrawing;
		const curveToActivate =
			lastSelectedCurve && activeCurves.includes(lastSelectedCurve)
				? lastSelectedCurve
				: selectPreferredCurve(activeCurves);
		set({ activeCurveForDrawing: curveToActivate, mode: "drawing" });
	},

	exitDrawingMode: async () => {
		const { mode } = get();
		if (mode !== "drawing") {
			console.warn(`Not in drawing mode, current mode: ${mode}`);
			return;
		}
		const released = await get().releaseLock();
		set({ mode: "view" });
		if (!released) {
			console.warn("Failed to release lock when exiting drawing mode");
		}
	},

	setDrawingPointerPosition: (position) =>
		set({ drawingPointerPosition: position }),

	setDrawingIsValid: (isValid) => set({ drawingIsValid: isValid }),

	setActiveCurveForDrawing: (key) => set({ activeCurveForDrawing: key }),
});
