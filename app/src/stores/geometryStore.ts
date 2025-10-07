import { create } from "zustand";
import { immer } from "zustand/middleware/immer";

type GeometryMode = "list" | "create" | "edit";

type GeometryState = {
  mode: GeometryMode;
  selectedKey: string | null;
  formData: Record<string, any>;
  selectedType: string | null;
  searchFilter: string;
};

type GeometryActions = {
  setMode: (mode: GeometryMode) => void;
  setSelectedKey: (key: string | null) => void;
  setFormData: (data: Record<string, any>) => void;
  setSelectedType: (type: string | null) => void;
  setSearchFilter: (filter: string) => void;
  resetForm: () => void;
};

export const useGeometryStore = create<GeometryState & GeometryActions>()(
  immer((set) => ({
    // Initialize the state
    mode: "list",
    selectedKey: null,
    formData: {},
    selectedType: null,
    searchFilter: "",

    // Define the actions
    setMode: (mode) => {
      set((state) => {
        state.mode = mode;
      });
    },

    setSelectedKey: (key) => {
      set((state) => {
        state.selectedKey = key;
      });
    },

    setFormData: (data) => {
      set((state) => {
        state.formData = data;
      });
    },

    setSelectedType: (type) => {
      set((state) => {
        state.selectedType = type;
      });
    },

    setSearchFilter: (filter) => {
      set((state) => {
        state.searchFilter = filter;
      });
    },

    resetForm: () => {
      set((state) => {
        state.formData = {};
        state.selectedKey = null;
        state.selectedType = null;
        state.mode = "list";
      });
    },
  }))
);
