// formStore.ts
import { create } from 'zustand';
import { immer } from 'zustand/middleware/immer';

type FormState = {
  formData: Record<string, any>;
  uiState: {
    selectedForm: string | null;
    selectedMethods: Record<string, string | null>;
  };
};

type FormActions = {
  setSelectedForm: (formId: string | null) => void;
  setSelectedMethod: (formId: string, methodId: string | null) => void;
  // This action now updates the local, temporary state of a form
  updateFormData: (compositeKey: string, data: any) => void;
  // New action to set the initial form data when it loads from the server
  initializeFormData: (compositeKey: string, data: any) => void;
};

export const useFormStore = create<FormState & FormActions>()(
  immer((set) => ({
    formConfigs: {},
    formData: {},
    uiState: {
      selectedForm: null,
      selectedMethods: {},
    },

    setSelectedForm: (formId) => {
      set((state) => {
        state.uiState.selectedForm = formId;
      });
    },

    setSelectedMethod: (formId, methodId) => {
      set((state) => {
        state.uiState.selectedMethods[formId] = methodId;
      });
    },

    updateFormData: (compositeKey, data) => {
      set((state) => {
        state.formData[compositeKey] = data;
      });
    },

    initializeFormData: (compositeKey, data) => {
        set((state) => {
            // Only initialize if it doesn't already have user edits
            if (state.formData[compositeKey] === undefined) {
                state.formData[compositeKey] = data;
            }
        });
    },
  }))
);