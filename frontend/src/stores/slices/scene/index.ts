import type { StateCreator } from "zustand";
import type { AppState } from "../../../store";

import type { DrawingSubSlice } from "./drawingSubSlice";
import { createDrawingSubSlice } from "./drawingSubSlice";
import type { EditingSubSlice } from "./editingSubSlice";
import { createEditingSubSlice } from "./editingSubSlice";
import type { GeometrySubSlice } from "./geometrySubSlice";
import { createGeometrySubSlice } from "./geometrySubSlice";
import type { SelectionSubSlice } from "./selectionSubSlice";
import { createSelectionSubSlice } from "./selectionSubSlice";

export type SceneSlice = GeometrySubSlice &
	SelectionSubSlice &
	EditingSubSlice &
	DrawingSubSlice;

export const createSceneSlice: StateCreator<AppState, [], [], SceneSlice> = (
	set,
	get,
	store,
) => ({
	...createGeometrySubSlice(set, get, store),
	...createSelectionSubSlice(set, get, store),
	...createEditingSubSlice(set, get, store),
	...createDrawingSubSlice(set, get, store),
});

export { getActiveCurves, selectPreferredCurve } from "./geometrySubSlice";
