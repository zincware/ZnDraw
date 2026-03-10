import type * as THREE from "three";
import type { StateCreator } from "zustand";
import { partialUpdateFrame } from "../../../myapi/client";
import type { AppState } from "../../../store";

export interface EditingSubSlice {
	transformMode: "translate" | "rotate" | "scale";
	editingSelectedAxis: "x" | "y" | "z" | null;
	editingCombinedCentroid: [number, number, number] | null;
	editingCallbacks: Map<string, Set<(matrix: THREE.Matrix4) => void>>;
	pendingFrameEdits: {
		frameId: number;
		keys: Record<string, any>;
	} | null;
	editingFrameDataCount: number;

	enterEditingMode: () => Promise<void>;
	exitEditingMode: () => Promise<void>;
	setTransformMode: (mode: "translate" | "rotate" | "scale") => void;
	cycleTransformMode: () => Promise<void>;
	setEditingSelectedAxis: (axis: "x" | "y" | "z" | null) => void;
	setEditingCombinedCentroid: (
		centroid: [number, number, number] | null,
	) => void;
	subscribeToEditing: (
		geometryKey: string,
		callback: (matrix: THREE.Matrix4) => void,
	) => () => void;
	notifyEditingChange: (matrix: THREE.Matrix4) => void;
	setPendingFrameEdit: (frameId: number, key: string, data: any) => void;
	clearPendingFrameEdits: () => void;
	saveFrameEdits: () => Promise<void>;
	incrementEditingFrameDataCount: () => void;
	decrementEditingFrameDataCount: () => void;
}

export const createEditingSubSlice: StateCreator<
	AppState,
	[],
	[],
	EditingSubSlice
> = (set, get) => ({
	transformMode: "translate",
	editingSelectedAxis: null,
	editingCombinedCentroid: null,
	editingCallbacks: new Map(),
	pendingFrameEdits: null,
	editingFrameDataCount: 0,

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
});
