// formStore.ts
import { create } from "zustand";
import { immer } from "zustand/middleware/immer";

type UiState = {
  selectedCategory: string | null;
  selectedExtensions: Record<string, string | null>;
};

type UiActions = {
  setSelectedCategory: (category: string | null) => void;
  setSelectedExtension: (category: string, extension: string | null) => void;
};

export const useFormStore = create<UiState & UiActions>()(
  immer((set) => ({
    // Initialize the state
    selectedCategory: null,
    selectedExtensions: {},

    // Define the actions
    setSelectedCategory: (category) => {
      set((state) => {
        state.selectedCategory = category;
      });
    },

    setSelectedExtension: (category, extension) => {
      set((state) => {
        state.selectedExtensions[category] = extension;
      });
    },
  })),
);
