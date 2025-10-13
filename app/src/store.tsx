import { create } from "zustand";
import { socket } from "./socket";
import * as THREE from "three";
import { updateSelection as updateSelectionAPI, loadSelectionGroup as loadSelectionGroupAPI } from "./myapi/client";

interface AppState {
  // Connection & Room
  roomId: string | null;
  userId: string | null;
  clientId: string | null;
  isConnected: boolean;
  isLoading: boolean;
  currentFrame: number;
  frameCount: number;
  skipFrames: number;
  selections: Record<string, number[]>; // per-geometry selections
  selectionGroups: Record<string, Record<string, number[]>>; // named selection groups
  activeSelectionGroup: string | null; // currently active group
  frame_selection: number[] | null;
  frameSelectionEnabled: boolean;
  bookmarks: Record<number, string> | null;
  playing: boolean;
  chatOpen: boolean;
  joinToken: string | null;
  geometries: Record<string, any>; // Store geometries by their IDs
  isDrawing: boolean;
  drawingPointerPosition: THREE.Vector3 | null; // 3D position of mouse cursor for drawing
  drawingIsValid: boolean; // Whether the drawing position is valid (over geometry)
  metadataLock: boolean; // Whether the room has a metadata lock (vis.lock - yellow)
  geometryFetchingStates: Record<string, boolean>; // Tracks fetching state per geometry key
  synchronizedMode: boolean; // Whether playback should wait for all active geometries to finish fetching
  showInfoBoxes: boolean; // Whether to show info boxes (toggled with 'i' key)
  hoveredParticleId: number | null; // ID of currently hovered particle
  particleCount: number; // Number of particles in current frame
  curveLength: number; // Length of the active curve in drawing mode

  // Actions (functions to modify the state)
  setRoomId: (roomId: string) => void;
  setUserId: (userId: string) => void;
  setConnected: (status: boolean) => void;
  setClientId: (clientId: string) => void;
  setJoinToken: (joinToken: string) => void;
  setCurrentFrame: (frame: number) => void;
  setFrameCount: (count: number) => void;
  setLoading: (loading: boolean) => void;
  setSkipFrames: (skip: number) => void;
  setSelections: (selections: Record<string, number[]>) => void; // set all selections
  updateSelectionForGeometry: (geometry: string, indices: number[]) => void; // update specific geometry
  setSelectionGroups: (groups: Record<string, Record<string, number[]>>) => void; // set all groups
  setActiveSelectionGroup: (groupName: string | null) => void; // set active group
  loadSelectionGroup: (groupName: string) => void; // load a group
  setFrameSelection: (selection: number[] | null) => void;
  setFrameSelectionEnabled: (enabled: boolean) => void;
  setBookmarks: (bookmark: Record<number, string> | null) => void;
  setPlaying: (playing: boolean) => void;
  setChatOpen: (open: boolean) => void;
  setGeometries: (geometries: Record<string, any>) => void;
  setIsDrawing: (isDrawing: boolean) => void;
  setDrawingPointerPosition: (position: THREE.Vector3 | null) => void;
  updateSelections: (geometryKey: string, id: number, isShiftPressed: boolean) => void;
  setDrawingIsValid: (isValid: boolean) => void;
  setMetadataLock: (locked: boolean) => void;
  setGeometryFetching: (geometryKey: string, isFetching: boolean) => void;
  removeGeometryFetching: (geometryKey: string) => void;
  getIsFetching: () => boolean; // Computed: returns true if any active geometry is fetching
  setSynchronizedMode: (enabled: boolean) => void;
  addBookmark: (frame: number, label?: string) => void;
  deleteBookmark: (frame: number) => void;
  toggleInfoBoxes: () => void;
  setHoveredParticleId: (id: number | null) => void;
  setParticleCount: (count: number) => void;
  setCurveLength: (length: number) => void;
}

export const useAppStore = create<AppState>((set, get) => ({
  // Initial State
  roomId: null,
  userId: null,
  clientId: null,
  isConnected: false,
  currentFrame: 0,
  frameCount: 0,
  isLoading: false,
  skipFrames: 1,
  selection: null,
  selections: {}, // New: per-geometry selections
  selectionGroups: {}, // New: named selection groups
  activeSelectionGroup: null, // New: currently active group
  frame_selection: null,
  frameSelectionEnabled: false,
  bookmarks: null,
  playing: false,
  chatOpen: false,
  joinToken: null,
  geometries: {},
  isDrawing: false,
  drawingPointerPosition: null,
  drawingIsValid: false,
  metadataLock: false,
  geometryFetchingStates: {},
  synchronizedMode: false,
  showInfoBoxes: false,
  hoveredParticleId: null,
  particleCount: 0,
  curveLength: 0,
  // Actions
  setConnected: (status) => set({ isConnected: status }),
  setRoomId: (roomId) => set({ roomId }),
  setUserId: (userId) => set({ userId }),
  setClientId: (clientId) => set({ clientId }),
  setCurrentFrame: (frame) => set({ currentFrame: frame }),
  setFrameCount: (count) =>
    set((state) => {
      // If currentFrame is beyond the new frameCount, reset to 0
      const newCurrentFrame =
        state.currentFrame >= count ? 0 : state.currentFrame;
      return { frameCount: count, currentFrame: newCurrentFrame };
    }),
  setLoading: (loading) => set({ isLoading: loading }),
  setSkipFrames: (skip) => set({ skipFrames: skip }),

  // New: set all selections (used when refetching from server)
  // TODO: socket emit / post request?
  setSelections: (selections) => {
    set({
      selections,
    });
  },

  // New: update selection for a specific geometry
  updateSelectionForGeometry: (geometry, indices) => {
    const roomId = get().roomId;
    if (!roomId) return;

    // Optimistic update
    set((state) => {
      const newSelections = {
        ...state.selections,
        [geometry]: indices,
      };
      return {
        selections: newSelections,
        selection: geometry === "particles" ? indices : state.selection,
      };
    });

    // Update via REST API
    updateSelectionAPI(roomId, geometry, indices).catch((error) => {
      console.error(`Failed to update selection for ${geometry}:`, error);
    });
  },

  // New: set all selection groups (used when refetching from server)
  setSelectionGroups: (groups) => set({ selectionGroups: groups }),

  // New: set the active selection group
  setActiveSelectionGroup: (groupName) => set({ activeSelectionGroup: groupName }),

  // New: load a selection group
  loadSelectionGroup: (groupName) => {
    const roomId = get().roomId;
    if (!roomId) return;

    // Call API to load the group - server will update selections and send invalidation
    loadSelectionGroupAPI(roomId, groupName).catch((error) => {
      console.error(`Failed to load selection group ${groupName}:`, error);
    });
  },

  setFrameSelection: (frame_selection) =>
    set({ frame_selection: frame_selection }),
  setFrameSelectionEnabled: (enabled) =>
    set({ frameSelectionEnabled: enabled }),
  setBookmarks: (bookmarks) => set({ bookmarks: bookmarks }),
  setPlaying: (playing) => set({ playing: playing }),
  setChatOpen: (open) => set({ chatOpen: open }),
  setJoinToken: (joinToken) => set({ joinToken }),
  setGeometries: (geometries) => set({ geometries: geometries }),
  updateGeometry: (key: string, geometry: any) =>
    set((state) => ({
      geometries: {
        ...state.geometries,
        [key]: geometry,
      },
    })),
  setIsDrawing: (isDrawing) => set({ isDrawing: isDrawing }),
  setDrawingPointerPosition: (position) => set({ drawingPointerPosition: position }),
  setDrawingIsValid: (isValid) => set({ drawingIsValid: isValid }),
  setMetadataLock: (locked) => set({ metadataLock: locked }),

  setGeometryFetching: (geometryKey, isFetching) =>
    set((state) => ({
      geometryFetchingStates: {
        ...state.geometryFetchingStates,
        [geometryKey]: isFetching,
      },
    })),

  removeGeometryFetching: (geometryKey) =>
    set((state) => {
      const newStates = { ...state.geometryFetchingStates };
      delete newStates[geometryKey];
      return { geometryFetchingStates: newStates };
    }),

  getIsFetching: () => {
    const state = get();
    const { geometries, geometryFetchingStates } = state;

    // Check if any active geometry is fetching
    return Object.entries(geometryFetchingStates).some(([key, isFetching]) => {
      const geometry = geometries[key];
      // Only count fetching if geometry exists and is active (or active is undefined for backwards compat)
      const isActive = geometry?.data?.active !== false;
      return isFetching && isActive;
    });
  },

  setSynchronizedMode: (enabled) => set({ synchronizedMode: enabled }),

  updateSelections: (geometryKey, id, isShiftPressed) => {
    const state = get();
    const roomId = state.roomId;
    if (!roomId) return;

    const currentSelection = state.selections[geometryKey] || [];
    // Use a Set for an efficient O(1) check to see if the item is already selected
    const isCurrentlySelected = new Set(currentSelection).has(id);

    let newSelection: number[];

    // --- LOGIC FOR SHIFT + CLICK ---
    if (isShiftPressed) {
      // Create a new Set from the current selection for easy manipulation
      const selectionSet = new Set(currentSelection);

      if (isCurrentlySelected) {
        // Rule 4: Shift-click on selected -> remove from selection set, keep the others
        selectionSet.delete(id);
      } else {
        // Rule 3: Shift-click on unselected -> add the element to selection
        selectionSet.add(id);
      }
      newSelection = Array.from(selectionSet);
    } else {
      // --- LOGIC FOR SIMPLE CLICK (no shift) ---
      if (isCurrentlySelected) {
        // Rule 2: Click on selected -> clear the selection
        newSelection = [];
      } else {
        // Rule 1: Click on unselected -> the element clicked on should be the ONLY element in selection
        newSelection = [id];
      }
    }

    // Optimistic update
    set((state) => ({
      selections: {
        ...state.selections,
        [geometryKey]: newSelection,
      },
    }));

    // Update via REST API
    updateSelectionAPI(roomId, geometryKey, newSelection).catch((error) => {
      console.error(`Failed to update selection for ${geometryKey}:`, error);
    });
  },

  addBookmark: (frame, label) => {
    const { bookmarks } = get();
    const currentBookmarks = bookmarks || {};

    // Add new bookmark
    const updatedBookmarks = {
      ...currentBookmarks,
      [frame]: label || `Bookmark ${frame}`
    };

    // Update local state immediately
    set({ bookmarks: updatedBookmarks });

    // Emit socket event to sync with other clients
    socket.emit("bookmarks:set", { bookmarks: updatedBookmarks });
  },

  deleteBookmark: (frame) => {
    const { bookmarks } = get();
    if (!bookmarks) return;

    // Remove bookmark
    const updatedBookmarks = { ...bookmarks };
    delete updatedBookmarks[frame];

    // Update local state immediately
    set({ bookmarks: updatedBookmarks });

    // Emit socket event to sync with other clients
    socket.emit("bookmarks:set", { bookmarks: updatedBookmarks });
  },

  toggleInfoBoxes: () => set((state) => ({ showInfoBoxes: !state.showInfoBoxes })),
  setHoveredParticleId: (id) => set({ hoveredParticleId: id }),
  setParticleCount: (count) => set({ particleCount: count }),
  setCurveLength: (length) => set({ curveLength: length }),
}));
