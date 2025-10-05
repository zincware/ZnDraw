import { create } from "zustand";
import { socket } from "./socket";

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
  selection: number[] | null;
  frame_selection: number[] | null;
  frameSelectionEnabled: boolean;
  bookmarks: Record<number, string> | null;
  playing: boolean;
  chatOpen: boolean;
  joinToken: string | null;
  geometries: Record<string, any>; // Store geometries by their IDs

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
  setSelection: (selection: number[] | null) => void;
  setFrameSelection: (selection: number[] | null) => void;
  setFrameSelectionEnabled: (enabled: boolean) => void;
  setBookmarks: (bookmark: Record<number, string> | null) => void;
  setPlaying: (playing: boolean) => void;
  setChatOpen: (open: boolean) => void;
  setGeometries: (geometries: Record<string, any>) => void;
  updateSelection: (id: number, isShiftPressed: boolean) => void;
}

export const useAppStore = create<AppState>((set) => ({
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
  frame_selection: null,
  frameSelectionEnabled: false,
  bookmarks: null,
  playing: false,
  chatOpen: false,
  joinToken: null,
  geometries: {},
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
  setSelection: (selection) => set({ selection: selection }),
  setFrameSelection: (frame_selection) =>
    set({ frame_selection: frame_selection }),
  setFrameSelectionEnabled: (enabled) =>
    set({ frameSelectionEnabled: enabled }),
  setBookmarks: (bookmarks) => set({ bookmarks: bookmarks }),
  setPlaying: (playing) => set({ playing: playing }),
  setChatOpen: (open) => set({ chatOpen: open }),
  setJoinToken: (joinToken) => set({ joinToken }),
  setGeometries: (geometries) => set({ geometries: geometries }),

  updateSelection: (id, isShiftPressed) =>
    set((state) => {
      const currentSelection = state.selection || [];
      // Use a Set for an efficient O(1) check to see if the item is already selected
      const isCurrentlySelected = new Set(currentSelection).has(id);

      // --- LOGIC FOR SHIFT + CLICK ---
      if (isShiftPressed) {
        // Create a new Set from the current selection for easy manipulation
        const selectionSet = new Set(currentSelection);
        
        if (isCurrentlySelected) {
          // Rule 4: Shift-click on selected -> remove from selection set, keep the others
          selectionSet.delete(id);
        } else {
          // Rule 3: Shift-click on unselected -> add the particle to selection
          selectionSet.add(id);
        }
        socket.emit("selection:set", {"indices": Array.from(selectionSet)});
        // Convert back to an array to store in state
        return { selection: Array.from(selectionSet) };
      }

      // --- LOGIC FOR SIMPLE CLICK (no shift) ---
      if (isCurrentlySelected) {
        // Rule 2: Click on selected -> no particles in the selection
        socket.emit("selection:set", {"indices": []});
        return { selection: [] };
      } else {
        // Rule 1: Click on unselected -> the particle clicked on should be the ONLY particle in selection
        socket.emit("selection:set", {"indices": [id]});
        return { selection: [id] };
      }
    }),
}));
