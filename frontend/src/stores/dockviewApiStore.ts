import type { DockviewApi } from "dockview-react";
import { create } from "zustand";

interface DockviewApiStore {
	api: DockviewApi | null;
	setApi: (api: DockviewApi | null) => void;
}

export const useDockviewApi = create<DockviewApiStore>((set) => ({
	api: null,
	setApi: (api) => set({ api }),
}));
