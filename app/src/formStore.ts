// formStore.ts
import { create } from 'zustand';
import { immer } from 'zustand/middleware/immer';

type UiState = {
  selectedForm: string | null; // Add this back
  selectedMethods: Record<string, string | null>;
};

type UiActions = {
  setSelectedForm: (formId: string | null) => void; // Add this back
  setSelectedMethod: (formId: string, methodId: string | null) => void;
};

export const useFormStore = create<UiState & UiActions>()(
  immer((set) => ({
    // Initialize the state
    selectedForm: null,
    selectedMethods: {},

    // Define the actions
    setSelectedForm: (formId) => {
      set((state) => {
        state.selectedForm = formId;
      });
    },

    setSelectedMethod: (formId, methodId) => {
      set((state) => {
        state.selectedMethods[formId] = methodId;
      });
    },
    
  }))
);