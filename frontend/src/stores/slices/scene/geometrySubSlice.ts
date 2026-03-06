import type * as THREE from "three";
import type { StateCreator } from "zustand";
import { createGeometry, getGeometry, updateActiveCamera } from "../../../myapi/client";
import type { AppState } from "../../../store";

/** Returns keys of active curves. */
export const getActiveCurves = (geometries: Record<string, any>): string[] => {
	return Object.entries(geometries)
		.filter(([_, g]) => g.type === "Curve" && g.data?.active !== false)
		.map(([key, _]) => key);
};

/** Selects the preferred curve from a list of active curves. */
export const selectPreferredCurve = (activeCurves: string[]): string | null => {
	if (activeCurves.length === 0) return null;
	return activeCurves.includes("curve") ? "curve" : activeCurves[0];
};

export interface GeometrySubSlice {
	geometries: Record<string, any>;
	geometrySchemas: Record<string, any>;
	geometryDefaults: Record<string, any>;
	geometryUpdateSources: Record<string, "local" | "remote">;
	geometryFetchingStates: Record<string, boolean>;
	neededFrameKeys: Record<string, string[]>;
	mode: "view" | "drawing" | "editing";
	attachedCameraKey: string | null;
	curveRefs: Record<string, THREE.Curve<THREE.Vector3>>;
	hoveredGeometryInstance: {
		geometryKey: string;
		instanceId: number;
	} | null;
	particleCount: number;
	curveLength: number;
	loadedDynamicPositions: Map<
		string,
		{ positionKey: string; positions: Float32Array }
	>;

	setGeometries: (geometries: Record<string, any>) => void;
	setGeometrySchemas: (schemas: Record<string, any>) => void;
	setGeometryDefaults: (defaults: Record<string, any>) => void;
	updateGeometry: (
		key: string,
		geometry: any,
		source?: "local" | "remote",
	) => void;
	removeGeometry: (key: string) => void;
	setGeometryFetching: (geometryKey: string, isFetching: boolean) => void;
	removeGeometryFetching: (geometryKey: string) => void;
	registerFrameKeys: (geometryKey: string, keys: string[]) => void;
	unregisterFrameKeys: (geometryKey: string) => void;
	getIsFetching: () => boolean;
	setMode: (mode: "view" | "drawing" | "editing") => void;
	attachToCamera: (cameraKey: string) => void;
	setAttachedCameraKey: (cameraKey: string) => void;
	registerCurveRef: (key: string, curve: THREE.Curve<THREE.Vector3>) => void;
	unregisterCurveRef: (key: string) => void;
	setHoveredGeometryInstance: (
		geometryKey: string | null,
		instanceId: number | null,
	) => void;
	setParticleCount: (count: number) => void;
	setCurveLength: (length: number) => void;
	registerLoadedDynamicPositions: (
		geometryKey: string,
		positionKey: string,
		positions: Float32Array,
	) => void;
	unregisterLoadedDynamicPositions: (geometryKey: string) => void;
}

export const createGeometrySubSlice: StateCreator<
	AppState,
	[],
	[],
	GeometrySubSlice
> = (set, get) => ({
	geometries: {},
	geometrySchemas: {},
	geometryDefaults: {},
	geometryUpdateSources: {},
	geometryFetchingStates: {},
	neededFrameKeys: {},
	mode: "view",
	attachedCameraKey: null,
	curveRefs: {},
	hoveredGeometryInstance: null,
	particleCount: 0,
	curveLength: 0,
	loadedDynamicPositions: new Map(),

	setGeometries: (geometries) => set({ geometries }),
	setGeometrySchemas: (schemas) => set({ geometrySchemas: schemas }),
	setGeometryDefaults: (defaults) => set({ geometryDefaults: defaults }),

	updateGeometry: (key, geometry, source = "remote") =>
		set((state) => ({
			geometries: { ...state.geometries, [key]: geometry },
			geometryUpdateSources: {
				...state.geometryUpdateSources,
				[key]: source,
			},
		})),

	removeGeometry: (key) =>
		set((state) => {
			const { [key]: removed, ...rest } = state.geometries;
			const { [key]: removedSource, ...restSources } =
				state.geometryUpdateSources;
			const { [key]: removedSelection, ...restSelections } = state.selections;

			let newActiveCurve = state.activeCurveForDrawing;
			if (state.activeCurveForDrawing === key) {
				const remainingCurves = getActiveCurves(rest);
				newActiveCurve = selectPreferredCurve(remainingCurves);
			}

			const newAttachedCameraKey =
				state.attachedCameraKey === key ? null : state.attachedCameraKey;

			const newEditingCallbacks = new Map(state.editingCallbacks);
			newEditingCallbacks.delete(key);

			return {
				geometries: rest,
				geometryUpdateSources: restSources,
				selections: restSelections,
				activeCurveForDrawing: newActiveCurve,
				attachedCameraKey: newAttachedCameraKey,
				editingCallbacks: newEditingCallbacks,
			};
		}),

	setGeometryFetching: (geometryKey, isFetching) =>
		set((state) => ({
			geometryFetchingStates: {
				...state.geometryFetchingStates,
				[geometryKey]: isFetching,
			},
		})),

	removeGeometryFetching: (geometryKey) =>
		set((state) => {
			const newStates = { ...state.geometryFetchingStates };
			delete newStates[geometryKey];
			return { geometryFetchingStates: newStates };
		}),

	registerFrameKeys: (geometryKey, keys) =>
		set((state) => ({
			neededFrameKeys: { ...state.neededFrameKeys, [geometryKey]: keys },
		})),

	unregisterFrameKeys: (geometryKey) =>
		set((state) => {
			const { [geometryKey]: _, ...rest } = state.neededFrameKeys;
			return { neededFrameKeys: rest };
		}),

	getIsFetching: () => {
		const { geometries, geometryFetchingStates } = get();
		return Object.entries(geometryFetchingStates).some(([key, isFetching]) => {
			const geometry = geometries[key];
			const isActive = geometry?.data?.active !== false;
			return isFetching && isActive;
		});
	},

	setMode: (mode) => set({ mode }),

	attachToCamera: (cameraKey) => {
		const { geometries, roomId, sessionId } = get();
		const camera = geometries[cameraKey];
		if (!camera || camera.type !== "Camera") {
			console.warn(
				`Cannot attach to camera '${cameraKey}': not found or not a Camera`,
			);
			return;
		}
		set({ attachedCameraKey: cameraKey });
		if (roomId && sessionId) {
			updateActiveCamera(roomId, sessionId, cameraKey).catch((error) => {
				console.error("Failed to sync active camera to server:", error);
			});
		}
	},

	setAttachedCameraKey: (cameraKey) => set({ attachedCameraKey: cameraKey }),

	registerCurveRef: (key, curve) =>
		set((state) => ({ curveRefs: { ...state.curveRefs, [key]: curve } })),

	unregisterCurveRef: (key) =>
		set((state) => {
			const { [key]: removed, ...rest } = state.curveRefs;
			return { curveRefs: rest };
		}),

	setHoveredGeometryInstance: (geometryKey, instanceId) => {
		if (geometryKey === null || instanceId === null) {
			set({ hoveredGeometryInstance: null });
		} else {
			set({ hoveredGeometryInstance: { geometryKey, instanceId } });
		}
	},

	setParticleCount: (count) => set({ particleCount: count }),
	setCurveLength: (length) => set({ curveLength: length }),

	registerLoadedDynamicPositions: (geometryKey, positionKey, positions) => {
		const newMap = new Map(get().loadedDynamicPositions);
		newMap.set(geometryKey, { positionKey, positions });
		set({ loadedDynamicPositions: newMap });
	},

	unregisterLoadedDynamicPositions: (geometryKey) => {
		const newMap = new Map(get().loadedDynamicPositions);
		newMap.delete(geometryKey);
		set({ loadedDynamicPositions: newMap });
	},
});
