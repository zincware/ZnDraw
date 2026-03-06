import type * as THREE from "three";
import type { StateCreator } from "zustand";
import {
	acquireEditLock,
	createGeometry,
	getGeometry,
	partialUpdateFrame,
	releaseEditLock,
	updateActiveCamera,
	updateSelection as updateSelectionAPI,
} from "../../myapi/client";
import type { AppState } from "../../store";

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

export interface SceneSlice {
	geometries: Record<string, any>;
	geometrySchemas: Record<string, any>;
	geometryDefaults: Record<string, any>;
	geometryUpdateSources: Record<string, "local" | "remote">;
	geometryFetchingStates: Record<string, boolean>;
	neededFrameKeys: Record<string, string[]>;
	selections: Record<string, number[]>;
	selectionGroups: Record<string, Record<string, number[]>>;
	mode: "view" | "drawing" | "editing";
	transformMode: "translate" | "rotate" | "scale";
	editingSelectedAxis: "x" | "y" | "z" | null;
	drawingPointerPosition: THREE.Vector3 | null;
	drawingIsValid: boolean;
	editingCombinedCentroid: [number, number, number] | null;
	editingCallbacks: Map<string, Set<(matrix: THREE.Matrix4) => void>>;
	activeCurveForDrawing: string | null;
	attachedCameraKey: string | null;
	curveRefs: Record<string, THREE.CatmullRomCurve3>;
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
	pendingFrameEdits: {
		frameId: number;
		keys: Record<string, any>;
	} | null;
	editingFrameDataCount: number;

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
	setSelections: (selections: Record<string, number[]>) => void;
	updateSelectionForGeometry: (geometry: string, indices: number[]) => void;
	setSelectionGroups: (
		groups: Record<string, Record<string, number[]>>,
	) => void;
	loadSelectionGroup: (groupName: string) => void;
	updateSelections: (
		geometryKey: string,
		id: number,
		isShiftPressed: boolean,
	) => void;
	setMode: (mode: "view" | "drawing" | "editing") => void;
	enterDrawingMode: (queryClient?: any) => Promise<void>;
	exitDrawingMode: () => Promise<void>;
	enterEditingMode: () => Promise<void>;
	exitEditingMode: () => Promise<void>;
	setTransformMode: (mode: "translate" | "rotate" | "scale") => void;
	cycleTransformMode: () => Promise<void>;
	setEditingSelectedAxis: (axis: "x" | "y" | "z" | null) => void;
	setDrawingPointerPosition: (position: THREE.Vector3 | null) => void;
	setDrawingIsValid: (isValid: boolean) => void;
	setEditingCombinedCentroid: (
		centroid: [number, number, number] | null,
	) => void;
	subscribeToEditing: (
		geometryKey: string,
		callback: (matrix: THREE.Matrix4) => void,
	) => () => void;
	notifyEditingChange: (matrix: THREE.Matrix4) => void;
	setActiveCurveForDrawing: (key: string | null) => void;
	attachToCamera: (cameraKey: string) => void;
	setAttachedCameraKey: (cameraKey: string) => void;
	registerCurveRef: (key: string, curve: THREE.CatmullRomCurve3) => void;
	unregisterCurveRef: (key: string) => void;
	setHoveredGeometryInstance: (
		geometryKey: string | null,
		instanceId: number | null,
	) => void;
	setParticleCount: (count: number) => void;
	setCurveLength: (length: number) => void;
	setPendingFrameEdit: (frameId: number, key: string, data: any) => void;
	clearPendingFrameEdits: () => void;
	saveFrameEdits: () => Promise<void>;
	incrementEditingFrameDataCount: () => void;
	decrementEditingFrameDataCount: () => void;
	registerLoadedDynamicPositions: (
		geometryKey: string,
		positionKey: string,
		positions: Float32Array,
	) => void;
	unregisterLoadedDynamicPositions: (geometryKey: string) => void;
}

export const createSceneSlice: StateCreator<AppState, [], [], SceneSlice> = (
	set,
	get,
) => ({
	geometries: {},
	geometrySchemas: {},
	geometryDefaults: {},
	geometryUpdateSources: {},
	geometryFetchingStates: {},
	neededFrameKeys: {},
	selections: {},
	selectionGroups: {},
	mode: "view",
	transformMode: "translate",
	editingSelectedAxis: null,
	drawingPointerPosition: null,
	drawingIsValid: false,
	editingCombinedCentroid: null,
	editingCallbacks: new Map(),
	activeCurveForDrawing: null,
	attachedCameraKey: null,
	curveRefs: {},
	hoveredGeometryInstance: null,
	particleCount: 0,
	curveLength: 0,
	loadedDynamicPositions: new Map(),
	pendingFrameEdits: null,
	editingFrameDataCount: 0,

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

	setSelections: (selections) => set({ selections }),

	updateSelectionForGeometry: (geometry, indices) => {
		const roomId = get().roomId;
		if (!roomId) return;
		set((state) => ({
			selections: { ...state.selections, [geometry]: indices },
		}));
		updateSelectionAPI(roomId, geometry, indices).catch((error) => {
			console.error(`Failed to update selection for ${geometry}:`, error);
		});
	},

	setSelectionGroups: (groups) => set({ selectionGroups: groups }),

	loadSelectionGroup: (groupName) => {
		const { selectionGroups } = get();
		const group = selectionGroups[groupName];
		if (!group) return;
		for (const [geometryKey, indices] of Object.entries(group)) {
			get().updateSelectionForGeometry(geometryKey, indices);
		}
	},

	updateSelections: (geometryKey, id, isShiftPressed) => {
		const state = get();
		const roomId = state.roomId;
		if (!roomId) return;

		const currentSelection = state.selections[geometryKey] || [];
		const isCurrentlySelected = new Set(currentSelection).has(id);

		let newSelection: number[];

		if (isShiftPressed) {
			const selectionSet = new Set(currentSelection);
			if (isCurrentlySelected) {
				selectionSet.delete(id);
			} else {
				selectionSet.add(id);
			}
			newSelection = Array.from(selectionSet);
		} else {
			if (isCurrentlySelected) {
				newSelection = [];
			} else {
				newSelection = [id];
			}
		}

		set((state) => ({
			selections: { ...state.selections, [geometryKey]: newSelection },
		}));

		updateSelectionAPI(roomId, geometryKey, newSelection).catch((error) => {
			console.error(`Failed to update selection for ${geometryKey}:`, error);
		});
	},

	setMode: (mode) => set({ mode }),

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

	enterEditingMode: async () => {
		const state = get();
		if (state.mode !== "view") {
			console.warn(`Cannot enter editing mode from ${state.mode} mode`);
			return;
		}
		const acquired = await get().acquireLock("editing geometries");
		if (acquired) {
			set({ mode: "editing" });
			state.showSnackbar("Entered editing mode", "success");
		}
	},

	exitEditingMode: async () => {
		const state = get();
		if (state.mode !== "editing") {
			console.warn(`Not in editing mode, current mode: ${state.mode}`);
			return;
		}
		await get().saveFrameEdits();
		for (const geometryKey of Object.keys(state.selections)) {
			get().updateSelectionForGeometry(geometryKey, []);
		}
		const released = await get().releaseLock();
		set({
			mode: "view",
			transformMode: "translate",
			editingSelectedAxis: null,
			editingFrameDataCount: 0,
		});
		if (released) {
			state.showSnackbar("Exited editing mode", "info");
		} else {
			console.warn("Failed to release lock when exiting editing mode");
		}
	},

	setTransformMode: (mode) => set({ transformMode: mode }),

	cycleTransformMode: async () => {
		await get().saveFrameEdits();
		const { transformMode, showSnackbar } = get();
		const modes: Array<"translate" | "rotate" | "scale"> = [
			"translate",
			"rotate",
			"scale",
		];
		const currentIndex = modes.indexOf(transformMode);
		const nextMode = modes[(currentIndex + 1) % modes.length];
		set({ transformMode: nextMode });
		showSnackbar(`Transform mode: ${nextMode}`, "info");
	},

	setEditingSelectedAxis: (axis) => set({ editingSelectedAxis: axis }),

	setDrawingPointerPosition: (position) =>
		set({ drawingPointerPosition: position }),

	setDrawingIsValid: (isValid) => set({ drawingIsValid: isValid }),

	setEditingCombinedCentroid: (centroid) =>
		set({ editingCombinedCentroid: centroid }),

	subscribeToEditing: (geometryKey, callback) => {
		const callbacks = get().editingCallbacks;
		if (!callbacks.has(geometryKey)) {
			callbacks.set(geometryKey, new Set());
		}
		callbacks.get(geometryKey)!.add(callback);
		return () => {
			const currentCallbacks = get().editingCallbacks.get(geometryKey);
			if (currentCallbacks) {
				currentCallbacks.delete(callback);
				if (currentCallbacks.size === 0) {
					get().editingCallbacks.delete(geometryKey);
				}
			}
		};
	},

	notifyEditingChange: (matrix) => {
		const callbacks = get().editingCallbacks;
		callbacks.forEach((callbackSet) => {
			callbackSet.forEach((callback) => callback(matrix));
		});
	},

	setActiveCurveForDrawing: (key) => set({ activeCurveForDrawing: key }),

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

	setPendingFrameEdit: (frameId, key, data) => {
		const state = get();

		if (
			state.pendingFrameEdits &&
			state.pendingFrameEdits.frameId !== frameId
		) {
			const oldEdits = state.pendingFrameEdits;
			void (async () => {
				const { roomId, showSnackbar } = get();
				if (!roomId || Object.keys(oldEdits.keys).length === 0) return;
				try {
					const shapes = new Map<string, number[]>();
					for (const [k, d] of Object.entries(oldEdits.keys)) {
						if (k === "arrays.positions" && d instanceof Float32Array) {
							shapes.set(k, [d.length / 3, 3]);
						}
					}
					await partialUpdateFrame(
						roomId,
						oldEdits.frameId,
						oldEdits.keys,
						shapes,
					);
					showSnackbar(`Auto-saved frame ${oldEdits.frameId}`, "info");
				} catch (error) {
					console.error("[setPendingFrameEdit] Auto-save failed:", error);
					showSnackbar(
						`Failed to auto-save frame ${oldEdits.frameId}`,
						"warning",
					);
				}
			})();
		}

		set({
			pendingFrameEdits: {
				frameId,
				keys: {
					...(state.pendingFrameEdits?.frameId === frameId
						? state.pendingFrameEdits.keys
						: {}),
					[key]: data,
				},
			},
		});
	},

	clearPendingFrameEdits: () => set({ pendingFrameEdits: null }),

	saveFrameEdits: async () => {
		const { roomId, pendingFrameEdits, showSnackbar } = get();
		if (!roomId || !pendingFrameEdits) return;

		const { frameId, keys } = pendingFrameEdits;
		if (Object.keys(keys).length === 0) return;

		try {
			const shapes = new Map<string, number[]>();
			for (const [key, data] of Object.entries(keys)) {
				if (key === "arrays.positions" && data instanceof Float32Array) {
					shapes.set(key, [data.length / 3, 3]);
				}
			}
			await partialUpdateFrame(roomId, frameId, keys, shapes);
			set({ pendingFrameEdits: null });
			showSnackbar(`Frame ${frameId} saved`, "success");
		} catch (error) {
			console.error("[saveFrameEdits] Error saving frame edits:", error);
			showSnackbar("Failed to save frame edits", "error");
		}
	},

	incrementEditingFrameDataCount: () =>
		set((state) => ({
			editingFrameDataCount: state.editingFrameDataCount + 1,
		})),

	decrementEditingFrameDataCount: () =>
		set((state) => ({
			editingFrameDataCount: Math.max(0, state.editingFrameDataCount - 1),
		})),

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
