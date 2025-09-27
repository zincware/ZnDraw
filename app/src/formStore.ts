// formStore.ts
import { create } from 'zustand';
import { immer } from 'zustand/middleware/immer';

// The structure received from the server
// e.g., { selection: { all: { schema: ..., uischema: ... } } }
type FormConfigs = Record<string, Record<string, { schema: any; uischema: any }>>;

type FormState = {
  // 1. Holds the definitions from the server
  formConfigs: FormConfigs;

  // 2. Holds user input, keyed by a composite ID like "form_a.method1"
  formData: Record<string, any>;

  // 3. Holds the state of the select dropdowns
  uiState: {
    selectedForm: string | null;
    selectedMethods: Record<string, string | null>;
  };
};

type FormActions = {
  // Action to load the definitions from the server
  setFormConfigs: (configs: FormConfigs) => void;

  // Actions to manage UI selections
  setSelectedForm: (formId: string | null) => void;
  setSelectedMethod: (formId: string, methodId: string | null) => void;

  // Action to update the form data itself
  updateFormData: (formId: string, methodId: string, data: any) => void;
};

export const useFormStore = create<FormState & FormActions>()(
  immer((set) => ({
    formConfigs: {},
    formData: {},
    uiState: {
      selectedForm: null,
      selectedMethods: {},
    },

    setFormConfigs: (configs) => {
      set((state) => {
        state.formConfigs = configs;
      });
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

    updateFormData: (formId, methodId, data) => {
      const compositeKey = `${formId}.${methodId}`;
      set((state) => {
        state.formData[compositeKey] = data;
      });
    },
  }))
);