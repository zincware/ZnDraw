import { create } from "zustand";
import { socket } from "./socket";
import * as THREE from "three";
import { updateSelection as updateSelectionAPI, loadSelectionGroup as loadSelectionGroupAPI, setBookmark as setBookmarkAPI, deleteBookmark as deleteBookmarkAPI, createGeometry, getGeometry, type GlobalSettings, acquireLock, refreshLock, releaseLock } from "./myapi/client";
import type { UserRole } from "./utils/auth";

interface AppState {
  // Connection & Room
  roomId: string | null;
  userName: string | null;
  userRole: UserRole | null;
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
  geometries: Record<string, any>; // Store geometries by their IDs
  geometryDefaults: Record<string, any>; // Default values for each geometry type from Pydantic
  geometryUpdateSources: Record<string, 'local' | 'remote'>; // Track update source per geometry
  isDrawing: boolean;
  drawingPointerPosition: THREE.Vector3 | null; // 3D position of mouse cursor for drawing
  drawingIsValid: boolean; // Whether the drawing position is valid (over geometry)
  isEditing: boolean; // Whether the user is in geometry editing mode (transform controls)
  editingCallbacks: Map<string, Set<(matrix: THREE.Matrix4) => void>>; // Callbacks for geometry components to subscribe to transform changes, keyed by geometryKey
  lockMetadata: {
    locked: boolean;
    holder?: string;
    userName?: string;
    msg?: string;
    timestamp?: number;
    ttl?: number;
  } | null; // Lock metadata for trajectory:meta lock (vis.lock)
  lockRenewalIntervalId: number | null; // Interval ID for lock renewal
  lockRefreshInterval: number | null; // Server-provided refresh interval in milliseconds
  geometryFetchingStates: Record<string, boolean>; // Tracks fetching state per geometry key
  synchronizedMode: boolean; // Whether playback should wait for all active geometries to finish fetching
  showInfoBoxes: boolean; // Whether to show info boxes (toggled with 'i' key)
  hoveredGeometryInstance: { geometryKey: string; instanceId: number } | null; // Currently hovered geometry instance
  hoveredFrame: number | null; // Currently hovered frame (from plot hover)
  particleCount: number; // Number of particles in current frame
  curveLength: number; // Length of the active curve in drawing mode
  fps: number | null; // Current FPS during playback
  frameLoadTime: number | null; // Time to load current frame (ms) when not playing
  lastFrameChangeTime: number | null; // Timestamp of last frame change (for FPS calculation)
  activeCurveForDrawing: string | null; // The geometry key of the curve currently targeted for drawing
  attachedCameraKey: string | null; // The geometry key of the camera currently attached to
  curveRefs: Record<string, THREE.CatmullRomCurve3>; // Non-serializable refs to THREE.js curve objects
  pathtracingNeedsUpdate: boolean; // Flag to signal pathtracer that scene has changed
  chatUnreadCount: number; // Number of unread chat messages
  serverVersion: string | null; // Server version for display and compatibility checking
  globalSettings: GlobalSettings | null; // Global settings (simgen, etc.)
  snackbar: { open: boolean; message: string; severity: "success" | "info" | "warning" | "error" } | null; // Snackbar notification state

  // Actions (functions to modify the state)
  setRoomId: (roomId: string) => void;
  setUserName: (userName: string) => void;
  setUserRole: (role: UserRole) => void;
  setConnected: (status: boolean) => void;
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
  setGeometryDefaults: (defaults: Record<string, any>) => void;
  updateGeometry: (key: string, geometry: any, source?: 'local' | 'remote') => void; // update specific geometry
  geometryUpdateSources: Record<string, 'local' | 'remote'>; // track update source per geometry
  removeGeometry: (key: string) => void; // remove specific geometry
  setIsDrawing: (isDrawing: boolean) => void;
  setDrawingPointerPosition: (position: THREE.Vector3 | null) => void;
  updateSelections: (geometryKey: string, id: number, isShiftPressed: boolean) => void;
  setDrawingIsValid: (isValid: boolean) => void;
  setIsEditing: (isEditing: boolean) => void;
  editingCombinedCentroid: [number, number, number] | null;
  setEditingCombinedCentroid: (centroid: [number, number, number] | null) => void;
  subscribeToEditing: (geometryKey: string, callback: (matrix: THREE.Matrix4) => void) => () => void;
  notifyEditingChange: (matrix: THREE.Matrix4) => void;
  toggleEditingMode: () => Promise<void>;
  setLockMetadata: (metadata: {
    locked: boolean;
    holder?: string;
    userName?: string;
    msg?: string;
    timestamp?: number;
    ttl?: number;
  } | null) => void;
  startLockRenewal: () => void;
  stopLockRenewal: () => void;
  setGeometryFetching: (geometryKey: string, isFetching: boolean) => void;
  removeGeometryFetching: (geometryKey: string) => void;
  getIsFetching: () => boolean; // Computed: returns true if any active geometry is fetching
  setSynchronizedMode: (enabled: boolean) => void;
  addBookmark: (frame: number, label?: string) => void;
  deleteBookmark: (frame: number) => void;
  toggleInfoBoxes: () => void;
  setHoveredGeometryInstance: (geometryKey: string | null, instanceId: number | null) => void;
  setHoveredFrame: (frame: number | null) => void;
  setParticleCount: (count: number) => void;
  setCurveLength: (length: number) => void;
  setFps: (fps: number | null) => void;
  setFrameLoadTime: (time: number | null) => void;
  setLastFrameChangeTime: (time: number | null) => void;
  setActiveCurveForDrawing: (key: string | null) => void;
  toggleDrawingMode: (queryClient?: any) => Promise<void>;
  attachToCamera: (cameraKey: string) => void;
  detachFromCamera: () => void;
  registerCurveRef: (key: string, curve: THREE.CatmullRomCurve3) => void;
  unregisterCurveRef: (key: string) => void;
  requestPathtracingUpdate: () => void; // Signal that pathtracer needs to update its scene
  clearPathtracingUpdate: () => void; // Clear the update flag after update is processed
  incrementChatUnread: () => void;
  resetChatUnread: () => void;
  setServerVersion: (version: string | null) => void;
  setGlobalSettings: (settings: GlobalSettings | null) => void;
  showSnackbar: (message: string, severity?: "success" | "info" | "warning" | "error") => void;
  hideSnackbar: () => void;
}

// Helper functions (pure, exported for reuse across components)
export const getActiveCurves = (geometries: Record<string, any>): string[] => {
  return Object.entries(geometries)
    .filter(([_, g]) => g.type === "Curve" && g.data?.active !== false)
    .map(([key, _]) => key);
};

export const selectPreferredCurve = (activeCurves: string[]): string | null => {
  if (activeCurves.length === 0) return null;
  // Prefer "curve" if it exists, otherwise use first available
  return activeCurves.includes("curve") ? "curve" : activeCurves[0];
};

export const useAppStore = create<AppState>((set, get) => ({
  // Initial State
  roomId: null,
  userName: null,
  userRole: null,
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
  geometries: {},
  geometryDefaults: {},
  geometryUpdateSources: {},
  isDrawing: false,
  drawingPointerPosition: null,
  drawingIsValid: false,
  isEditing: false,
  editingCombinedCentroid: null,
  editingCallbacks: new Map(),
  lockMetadata: null,
  lockRenewalIntervalId: null,
  lockRefreshInterval: null,
  geometryFetchingStates: {},
  synchronizedMode: true,
  showInfoBoxes: false,
  hoveredGeometryInstance: null,
  hoveredFrame: null,
  particleCount: 0,
  curveLength: 0,
  fps: null,
  frameLoadTime: null,
  lastFrameChangeTime: null,
  activeCurveForDrawing: null,
  attachedCameraKey: null,
  pathtracingNeedsUpdate: false,
  chatUnreadCount: 0,
  serverVersion: null,
  globalSettings: null,
  snackbar: null,

  /**
   * Non-serializable THREE.js curve objects shared between Curve and Camera components.
   * ⚠️ Non-serializable - excluded from persistence/devtools.
   *
   * Managed via registerCurveRef/unregisterCurveRef lifecycle hooks in Curve.tsx.
   * Used by Camera.tsx to compute positions along curves without prop drilling.
   *
   * This is the correct pattern for Zustand - see https://docs.pmnd.rs/zustand/guides/flux-inspired-practice
   * If persistence is added later, use partialize() to exclude this field.
   */
  curveRefs: {},

  // Actions
  setConnected: (status) => set({ isConnected: status }),
  setRoomId: (roomId) => set({ roomId }),
  setUserName: (userName) => set({ userName }),
  setUserRole: (role) => set({ userRole: role }),
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
  setGeometries: (geometries) => set({ geometries: geometries }),
  setGeometryDefaults: (defaults) => set({ geometryDefaults: defaults }),
  updateGeometry: (key: string, geometry: any, source: 'local' | 'remote' = 'remote') =>
    set((state) => ({
      geometries: {
        ...state.geometries,
        [key]: geometry,
      },
      geometryUpdateSources: {
        ...state.geometryUpdateSources,
        [key]: source,
      },
    })),
  removeGeometry: (key: string) =>
    set((state) => {
      const { [key]: removed, ...rest } = state.geometries;
      const { [key]: removedSource, ...restSources } = state.geometryUpdateSources;
      const { [key]: removedSelection, ...restSelections } = state.selections;

      // If the deleted geometry was the active curve, auto-select next available curve
      let newActiveCurve = state.activeCurveForDrawing;
      if (state.activeCurveForDrawing === key) {
        const remainingCurves = getActiveCurves(rest);
        newActiveCurve = selectPreferredCurve(remainingCurves);
      }

      // Clean up editing callbacks for this geometry to prevent memory leaks
      const newEditingCallbacks = new Map(state.editingCallbacks);
      newEditingCallbacks.delete(key);

      return {
        geometries: rest,
        geometryUpdateSources: restSources,
        selections: restSelections,
        activeCurveForDrawing: newActiveCurve,
        editingCallbacks: newEditingCallbacks,
      };
    }),
  setIsDrawing: (isDrawing) => set({ isDrawing: isDrawing }),
  setDrawingPointerPosition: (position) => set({ drawingPointerPosition: position }),
  setDrawingIsValid: (isValid) => set({ drawingIsValid: isValid }),
  setIsEditing: (isEditing) => set({ isEditing: isEditing }),
  setEditingCombinedCentroid: (centroid) => set({ editingCombinedCentroid: centroid }),
  subscribeToEditing: (geometryKey, callback) => {
    const callbacks = get().editingCallbacks;

    if (!callbacks.has(geometryKey)) {
      callbacks.set(geometryKey, new Set());
    }

    callbacks.get(geometryKey)!.add(callback);

    // Return unsubscribe function
    return () => {
      const currentCallbacks = get().editingCallbacks.get(geometryKey);
      if (currentCallbacks) {
        currentCallbacks.delete(callback);
        // Clean up empty sets
        if (currentCallbacks.size === 0) {
          get().editingCallbacks.delete(geometryKey);
        }
      }
    };
  },
  notifyEditingChange: (matrix) => {
    const callbacks = get().editingCallbacks;
    // Iterate over all geometry callback sets
    callbacks.forEach((callbackSet) => {
      callbackSet.forEach((callback) => callback(matrix));
    });
  },
  setLockMetadata: (metadata) => set({ lockMetadata: metadata }),

  startLockRenewal: () => {
    const { lockRenewalIntervalId, lockRefreshInterval, roomId, stopLockRenewal } = get();

    // Clear existing interval if any
    if (lockRenewalIntervalId !== null) {
      stopLockRenewal();
    }

    // Use server-provided refresh interval (in seconds), default to 30s
    const refreshIntervalMs = (lockRefreshInterval || 30) * 1000;

    // Renew periodically to keep lock alive
    const intervalId = window.setInterval(async () => {
      const state = get();
      if (!state.isEditing || !state.roomId) {
        // Not editing anymore, stop renewal
        stopLockRenewal();
        return;
      }

      // Refresh lock via REST API (server controls TTL)
      try {
        const response = await refreshLock(state.roomId, "trajectory:meta");
        if (!response.success) {
          console.warn("Failed to refresh edit lock");
          // Lock was stolen or expired, exit editing mode
          get().setIsEditing(false);
          get().showSnackbar("Edit lock lost - exiting editing mode", "warning");
          get().stopLockRenewal();
        }
      } catch (error) {
        console.error("Error refreshing lock:", error);
        // Network error or other issue
        get().setIsEditing(false);
        get().showSnackbar("Failed to refresh lock - exiting editing mode", "error");
        get().stopLockRenewal();
      }
    }, refreshIntervalMs);

    set({ lockRenewalIntervalId: intervalId });
  },

  stopLockRenewal: () => {
    const { lockRenewalIntervalId } = get();
    if (lockRenewalIntervalId !== null) {
      window.clearInterval(lockRenewalIntervalId);
      set({ lockRenewalIntervalId: null });
    }
  },

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
    const roomId = get().roomId;
    if (!roomId) return;

    const bookmarkLabel = label || `Bookmark ${frame}`;

    // Optimistic update - add bookmark locally
    set((state) => ({
      bookmarks: {
        ...(state.bookmarks || {}),
        [frame]: bookmarkLabel,
      },
    }));

    // Update via REST API - server will emit invalidate event to sync with other clients
    setBookmarkAPI(roomId, frame, bookmarkLabel).catch((error) => {
      console.error(`Failed to set bookmark at frame ${frame}:`, error);
      // Rollback on error
      set((state) => {
        if (!state.bookmarks) return state;
        const { [frame]: removed, ...rest } = state.bookmarks;
        return { bookmarks: rest };
      });
    });
  },

  deleteBookmark: (frame) => {
    const roomId = get().roomId;
    if (!roomId) return;

    const { bookmarks } = get();
    if (!bookmarks || !(frame in bookmarks)) return;

    // Store the old value for rollback
    const oldLabel = bookmarks[frame];

    // Optimistic update - remove bookmark locally
    set((state) => {
      if (!state.bookmarks) return state;
      const { [frame]: removed, ...rest } = state.bookmarks;
      return { bookmarks: rest };
    });

    // Update via REST API - server will emit invalidate event to sync with other clients
    deleteBookmarkAPI(roomId, frame).catch((error) => {
      console.error(`Failed to delete bookmark at frame ${frame}:`, error);
      // Rollback on error
      set((state) => ({
        bookmarks: {
          ...(state.bookmarks || {}),
          [frame]: oldLabel,
        },
      }));
    });
  },

  toggleInfoBoxes: () => set((state) => ({ showInfoBoxes: !state.showInfoBoxes })),
  setHoveredGeometryInstance: (geometryKey, instanceId) => {
    if (geometryKey === null || instanceId === null) {
      set({ hoveredGeometryInstance: null });
    } else {
      set({ hoveredGeometryInstance: { geometryKey, instanceId } });
    }
  },
  setHoveredFrame: (frame) => set({ hoveredFrame: frame }),
  setParticleCount: (count) => set({ particleCount: count }),
  setCurveLength: (length) => set({ curveLength: length }),
  setFps: (fps) => set({ fps }),
  setFrameLoadTime: (time) => set({ frameLoadTime: time }),
  setLastFrameChangeTime: (time) => set({ lastFrameChangeTime: time }),
  setActiveCurveForDrawing: (key) => set({ activeCurveForDrawing: key }),

  toggleEditingMode: async () => {
    const state = get();
    const { isEditing, lockMetadata, userName, roomId } = state;

    // Import dynamically to avoid circular dependency
    const { socket } = await import("./socket");

    // If turning ON editing mode
    if (!isEditing) {
      // Check if room is already locked
      if (lockMetadata && lockMetadata.locked) {
        // Check if WE hold the lock
        if (lockMetadata.holder === userName) {
          // We already hold the lock, just enter editing mode
          set({ isEditing: true });
          state.showSnackbar("Entered editing mode", "success");

          // Start lock renewal
          get().startLockRenewal();
          return;
        }

        // Someone else holds the lock
        state.showSnackbar(
          `Cannot enter edit mode - room is locked by ${lockMetadata.userName || "another user"}`,
          "warning"
        );
        return;
      }

      // Try to acquire lock via REST API
      return (async () => {
        if (!roomId) {
          state.showSnackbar("Cannot acquire lock - no room ID", "error");
          return;
        }

        try {
          const response = await acquireLock(roomId, "trajectory:meta", "editing geometries");
          if (response.success) {
            // Store server-provided refresh interval
            if (response.refreshInterval) {
              set({ lockRefreshInterval: response.refreshInterval });
            }

            // Enter editing mode
            set({ isEditing: true });
            state.showSnackbar("Entered editing mode", "success");

            // Start lock renewal to keep lock alive
            get().startLockRenewal();
          } else {
            // Failed to acquire lock
            state.showSnackbar(
              response.error || "Cannot enter edit mode - room is locked by another user",
              "warning"
            );
          }
        } catch (error: any) {
          console.error("Error acquiring lock:", error);
          // Check if 423 Locked
          if (error.response?.status === 423) {
            state.showSnackbar(
              error.response.data?.error || "Cannot enter edit mode - room is locked by another user",
              "warning"
            );
          } else {
            state.showSnackbar("Failed to acquire lock", "error");
          }
        }
      })();
    }
    // If turning OFF editing mode
    else {
      // Stop lock renewal
      get().stopLockRenewal();

      // Release the lock via REST API
      return (async () => {
        if (!roomId) {
          set({ isEditing: false });
          return;
        }

        try {
          const response = await releaseLock(roomId, "trajectory:meta");
          // Exit editing mode regardless of release success
          set({ isEditing: false });

          if (response.success) {
            state.showSnackbar("Exited editing mode", "info");
          } else {
            console.warn("Failed to release lock, but exiting editing mode anyway");
          }
        } catch (error) {
          console.error("Error releasing lock:", error);
          // Exit editing mode regardless of release success
          set({ isEditing: false });
        }
      })();
    }
  },

  toggleDrawingMode: async (queryClient?: any) => {
    const state = get();
    const { geometries, activeCurveForDrawing, isDrawing, roomId } = state;

    // Get all active curves using helper function
    const activeCurves = getActiveCurves(geometries);

    // Case 1: No curves exist - create one named "curve"
    if (activeCurves.length === 0) {
      if (!roomId) {
        console.warn("Cannot create curve: no room ID");
        return;
      }

      console.log("No curves found, creating default curve named 'curve'");

      // Optimistic update - add geometry to store immediately
      set((state) => ({
        geometries: {
          ...state.geometries,
          curve: {
            type: "Curve",
            data: { position: [] }
          }
        },
        activeCurveForDrawing: "curve",
        isDrawing: true
      }));

      try {
        // Create geometry on backend
        await createGeometry(roomId, "curve", "Curve", {
          position: [],
        });
        console.log("Default curve created successfully");

        // Fetch the validated geometry from backend to get defaults/transformations
        // (socket event won't reach us due to skip_sid)
        const response = await getGeometry(roomId, "curve");
        console.log("Fetched validated curve data:", response.geometry);

        // Update with the backend-validated data
        set((state) => ({
          geometries: {
            ...state.geometries,
            curve: response.geometry
          }
        }));

        // Invalidate geometry list to refetch from server (this will update the sidebar)
        if (queryClient && roomId) {
          queryClient.invalidateQueries({
            queryKey: ["geometries", roomId, "list"]
          });
        }
      } catch (error) {
        console.error("Error creating default curve:", error);
        // Rollback on error - remove the optimistically added geometry
        set((state) => {
          const { curve: removed, ...rest } = state.geometries;
          return {
            geometries: rest,
            activeCurveForDrawing: null,
            isDrawing: false
          };
        });
      }
      return;
    }

    // Case 2: One curve exists - auto-select it and toggle
    if (activeCurves.length === 1) {
      const singleCurveKey = activeCurves[0];
      if (!activeCurveForDrawing) {
        set({ activeCurveForDrawing: singleCurveKey, isDrawing: !isDrawing });
      } else {
        set({ isDrawing: !isDrawing });
      }
      return;
    }

    // Case 3: Multiple curves exist
    if (!activeCurveForDrawing || !activeCurves.includes(activeCurveForDrawing)) {
      // No curve selected OR selected curve no longer exists - select preferred curve
      const preferredCurve = selectPreferredCurve(activeCurves);
      set({ activeCurveForDrawing: preferredCurve, isDrawing: true });
    } else {
      // Already have a valid selected curve - just toggle drawing
      set({ isDrawing: !isDrawing });
    }
  },

  attachToCamera: (cameraKey) => {
    const { geometries } = get();
    const camera = geometries[cameraKey];

    // Validate camera exists and is a Camera type
    if (!camera || camera.type !== "Camera") {
      console.warn(`Cannot attach to camera '${cameraKey}': not found or not a Camera`);
      return;
    }

    console.log(`Attaching to camera: ${cameraKey}`);
    set({ attachedCameraKey: cameraKey });
  },

  detachFromCamera: () => {
    console.log("Detaching from camera");
    set({ attachedCameraKey: null });
  },

  registerCurveRef: (key, curve) => {
    set((state) => ({
      curveRefs: {
        ...state.curveRefs,
        [key]: curve,
      },
    }));
  },

  unregisterCurveRef: (key) => {
    set((state) => {
      const { [key]: removed, ...rest } = state.curveRefs;
      return { curveRefs: rest };
    });
  },

  requestPathtracingUpdate: () => {
    set({ pathtracingNeedsUpdate: true });
  },

  clearPathtracingUpdate: () => {
    set({ pathtracingNeedsUpdate: false });
  },

  incrementChatUnread: () => set((state) => ({ chatUnreadCount: state.chatUnreadCount + 1 })),
  resetChatUnread: () => set({ chatUnreadCount: 0 }),
  setServerVersion: (version) => set({ serverVersion: version }),
  setGlobalSettings: (settings) => set({ globalSettings: settings }),
  showSnackbar: (message, severity = "info") => set({ snackbar: { open: true, message, severity } }),
  hideSnackbar: () => set((state) => state.snackbar ? { snackbar: { ...state.snackbar, open: false } } : {}),
}));
