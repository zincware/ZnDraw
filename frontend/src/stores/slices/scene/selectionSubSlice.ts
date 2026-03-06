import type { StateCreator } from "zustand";
import { updateSelection as updateSelectionAPI } from "../../../myapi/client";
import type { AppState } from "../../../store";

export interface SelectionSubSlice {
	selections: Record<string, number[]>;
	selectionGroups: Record<string, Record<string, number[]>>;

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
}

export const createSelectionSubSlice: StateCreator<
	AppState,
	[],
	[],
	SelectionSubSlice
> = (set, get) => ({
	selections: {},
	selectionGroups: {},

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
});
