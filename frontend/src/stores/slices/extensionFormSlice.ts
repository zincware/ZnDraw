import type { StateCreator } from "zustand";
import type { AppState } from "../../store";

export interface ExtensionFormSlice {
	selectedExtensions: Record<string, string | null>;
	setSelectedExtension: (category: string, extension: string | null) => void;
}

export const createExtensionFormSlice: StateCreator<
	AppState,
	[],
	[],
	ExtensionFormSlice
> = (set) => ({
	selectedExtensions: {},
	setSelectedExtension: (category, extension) =>
		set((state) => ({
			selectedExtensions: {
				...state.selectedExtensions,
				[category]: extension,
			},
		})),
});
