import { create } from "zustand";
import { socket } from "./socket";
import * as THREE from "three";
import {
	updateSelection as updateSelectionAPI,
	updateFrameSelection as updateFrameSelectionAPI,
	loadSelectionGroup as loadSelectionGroupAPI,
	setBookmark as setBookmarkAPI,
	deleteBookmark as deleteBookmarkAPI,
	createGeometry,
	getGeometry,
	type GlobalSettings,
	acquireLock,
	releaseLock,
	partialUpdateFrame,
} from "./myapi/client";
import { getToken, type UserRole } from "./utils/auth";

/**
 * Lock state for a specific target
 * Simplified design: if you have the lock token for a target, you can use it for anything
 */
interface LockState {
	target: string; // e.g., "trajectory:meta"
	token: string;
	refreshInterval: number | null;
	acquiredAt: number;
	message: string; // What we're doing with the lock (for display to other users)
}

/**
 * Progress tracking state
 */
export interface Progress {
	progressId: string;
	roomId: string; // Room this progress tracker belongs to
	description: string;
	progress: number | null; // 0-100 for progress bar, null for spinner
}

/**
 * Error state for room initialization failures
 */
export interface InitializationError {
	message: string;
	details?: string;
}

interface AppState {
	// Connection & Room
	roomId: string | null;
	userName: string | null;
	userRole: UserRole | null;
	sessionId: string | null; // Session ID from /join response (identifies this browser tab)
	isConnected: boolean;
	isLoading: boolean;
	initializationError: InitializationError | null; // Error during room initialization
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
	geometrySchemas: Record<string, any>; // Pydantic JSON schemas for form generation
	geometryDefaults: Record<string, any>; // Pydantic default values for each geometry type
	geometryUpdateSources: Record<string, "local" | "remote">; // Track update source per geometry
	mode: "view" | "drawing" | "editing"; // Current interaction mode
	transformMode: "translate" | "rotate" | "scale"; // Transform mode within editing mode
	editingSelectedAxis: "x" | "y" | "z" | null; // Axis selected for keyboard editing (hold X/Y/Z)
	drawingPointerPosition: THREE.Vector3 | null; // 3D position of mouse cursor for drawing
	drawingIsValid: boolean; // Whether the drawing position is valid (over geometry)
	editingCallbacks: Map<string, Set<(matrix: THREE.Matrix4) => void>>; // Callbacks for geometry components to subscribe to transform changes, keyed by geometryKey
	lockMetadata: {
		locked: boolean;
		holder?: string;
		userName?: string;
		msg?: string;
		timestamp?: number;
		ttl?: number;
	} | null; // Lock metadata for trajectory:meta lock (vis.lock)

	// Centralized lock state - single source of truth
	lock: LockState | null; // Current lock ownership
	lockRenewalIntervalId: number | null; // Interval ID for lock renewal
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
	snackbar: {
		open: boolean;
		message: string;
		severity: "success" | "info" | "warning" | "error";
	} | null; // Snackbar notification state
	progressTrackers: Record<string, Progress>; // Active progress trackers keyed by progressId

	// Frame editing state for dynamic positions
	pendingFrameEdits: {
		frameId: number;
		keys: Record<string, any>; // e.g., { "arrays.positions": Float32Array }
	} | null;
	// Reference count for frame editing consumers (multiple hooks can be editing simultaneously)
	editingFrameDataCount: number;
	// Loaded dynamic positions registered by geometry components
	// Used by MultiGeometryTransformControls to calculate centroid for dynamic positions
	loadedDynamicPositions: Map<
		string,
		{ positionKey: string; positions: Float32Array }
	>;

	// Actions (functions to modify the state)
	setRoomId: (roomId: string) => void;
	setUserName: (userName: string) => void;
	setUserRole: (role: UserRole) => void;
	setSessionId: (sessionId: string | null) => void;
	setConnected: (status: boolean) => void;
	setInitializationError: (error: InitializationError | null) => void;
	setCurrentFrame: (frame: number) => void;
	setFrameCount: (count: number) => void;
	setLoading: (loading: boolean) => void;
	setSkipFrames: (skip: number) => void;
	setSelections: (selections: Record<string, number[]>) => void; // set all selections
	updateSelectionForGeometry: (geometry: string, indices: number[]) => void; // update specific geometry
	setSelectionGroups: (
		groups: Record<string, Record<string, number[]>>,
	) => void; // set all groups
	setActiveSelectionGroup: (groupName: string | null) => void; // set active group
	loadSelectionGroup: (groupName: string) => void; // load a group
	setFrameSelection: (selection: number[] | null) => void; // local state only (for socket)
	updateFrameSelection: (selection: number[] | null) => void; // local + API (for user actions)
	setFrameSelectionEnabled: (enabled: boolean) => void;
	setBookmarks: (bookmark: Record<number, string> | null) => void;
	setPlaying: (playing: boolean) => void;
	setChatOpen: (open: boolean) => void;
	setGeometries: (geometries: Record<string, any>) => void;
	setGeometrySchemas: (schemas: Record<string, any>) => void;
	setGeometryDefaults: (defaults: Record<string, any>) => void;
	updateGeometry: (
		key: string,
		geometry: any,
		source?: "local" | "remote",
	) => void; // update specific geometry
	removeGeometry: (key: string) => void; // remove specific geometry
	setMode: (mode: "view" | "drawing" | "editing") => void;
	enterDrawingMode: (queryClient?: any) => Promise<void>;
	exitDrawingMode: () => Promise<void>;
	enterEditingMode: () => Promise<void>;
	exitEditingMode: () => Promise<void>;
	setTransformMode: (mode: "translate" | "rotate" | "scale") => void;
	cycleTransformMode: () => Promise<void>;
	setEditingSelectedAxis: (axis: "x" | "y" | "z" | null) => void;
	setDrawingPointerPosition: (position: THREE.Vector3 | null) => void;
	updateSelections: (
		geometryKey: string,
		id: number,
		isShiftPressed: boolean,
	) => void;
	setDrawingIsValid: (isValid: boolean) => void;
	editingCombinedCentroid: [number, number, number] | null;
	setEditingCombinedCentroid: (
		centroid: [number, number, number] | null,
	) => void;
	subscribeToEditing: (
		geometryKey: string,
		callback: (matrix: THREE.Matrix4) => void,
	) => () => void;
	notifyEditingChange: (matrix: THREE.Matrix4) => void;
	setLockMetadata: (
		metadata: {
			locked: boolean;
			holder?: string;
			userName?: string;
			msg?: string;
			timestamp?: number;
			ttl?: number;
		} | null,
	) => void;
	startLockRenewal: () => void;
	stopLockRenewal: () => void;

	// Simplified lock management - target-based (aligns with server API)
	acquireLock: (target: string, msg: string) => Promise<boolean>;
	releaseLock: (target: string) => Promise<boolean>;
	hasLock: (target: string) => boolean;
	updateLockMessage: (target: string, msg: string) => Promise<boolean>;
	setGeometryFetching: (geometryKey: string, isFetching: boolean) => void;
	removeGeometryFetching: (geometryKey: string) => void;
	getIsFetching: () => boolean; // Computed: returns true if any active geometry is fetching
	setSynchronizedMode: (enabled: boolean) => void;
	addBookmark: (frame: number, label?: string) => void;
	deleteBookmark: (frame: number) => void;
	toggleInfoBoxes: () => void;
	setHoveredGeometryInstance: (
		geometryKey: string | null,
		instanceId: number | null,
	) => void;
	setHoveredFrame: (frame: number | null) => void;
	setParticleCount: (count: number) => void;
	setCurveLength: (length: number) => void;
	setFps: (fps: number | null) => void;
	setFrameLoadTime: (time: number | null) => void;
	setLastFrameChangeTime: (time: number | null) => void;
	setActiveCurveForDrawing: (key: string | null) => void;
	attachToCamera: (cameraKey: string) => void;
	setAttachedCameraKey: (cameraKey: string) => void;
	registerCurveRef: (key: string, curve: THREE.CatmullRomCurve3) => void;
	unregisterCurveRef: (key: string) => void;
	requestPathtracingUpdate: () => void; // Signal that pathtracer needs to update its scene
	clearPathtracingUpdate: () => void; // Clear the update flag after update is processed
	incrementChatUnread: () => void;
	resetChatUnread: () => void;
	setServerVersion: (version: string | null) => void;
	setGlobalSettings: (settings: GlobalSettings | null) => void;
	showSnackbar: (
		message: string,
		severity?: "success" | "info" | "warning" | "error",
	) => void;
	hideSnackbar: () => void;
	setProgressTrackers: (trackers: Record<string, Progress>) => void;
	addProgressTracker: (
		progressId: string,
		description: string,
		progress: number | null,
		roomId: string,
	) => void;
	updateProgressTracker: (
		progressId: string,
		description?: string,
		progress?: number | null,
	) => void;
	removeProgressTracker: (progressId: string) => void;

	// Frame editing actions
	setPendingFrameEdit: (frameId: number, key: string, data: any) => void;
	clearPendingFrameEdits: () => void;
	saveFrameEdits: () => Promise<void>;
	incrementEditingFrameDataCount: () => void;
	decrementEditingFrameDataCount: () => void;
	registerLoadedDynamicPositions: (
		geometryKey: string,
		positionKey: string,
		positions: Float32Array,
	) => void;
	unregisterLoadedDynamicPositions: (geometryKey: string) => void;
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
	sessionId: null, // Will be set after /join
	isConnected: false,
	initializationError: null,
	currentFrame: 0,
	frameCount: 0,
	isLoading: false,
	skipFrames: 1,
	selections: {}, // per-geometry selections
	selectionGroups: {}, // New: named selection groups
	activeSelectionGroup: null, // New: currently active group
	frame_selection: null,
	frameSelectionEnabled: false,
	bookmarks: null,
	playing: false,
	chatOpen: false,
	geometries: {},
	geometrySchemas: {},
	geometryDefaults: {},
	geometryUpdateSources: {},
	mode: "view",
	transformMode: "translate",
	editingSelectedAxis: null,
	drawingPointerPosition: null,
	drawingIsValid: false,
	editingCombinedCentroid: null,
	editingCallbacks: new Map(),
	lockMetadata: null,
	lock: null, // Centralized lock state
	lockRenewalIntervalId: null,
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
	progressTrackers: {},

	// Frame editing state
	pendingFrameEdits: null,
	editingFrameDataCount: 0,
	loadedDynamicPositions: new Map(),

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
	setInitializationError: (error) => set({ initializationError: error }),
	setRoomId: (roomId) => set({ roomId }),
	setUserName: (userName) => set({ userName }),
	setUserRole: (role) => set({ userRole: role }),
	setSessionId: (sessionId) => set({ sessionId }),
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
	setActiveSelectionGroup: (groupName) =>
		set({ activeSelectionGroup: groupName }),

	// New: load a selection group
	loadSelectionGroup: (groupName) => {
		const roomId = get().roomId;
		if (!roomId) return;

		// Call API to load the group - server will update selections and send invalidation
		loadSelectionGroupAPI(roomId, groupName).catch((error) => {
			console.error(`Failed to load selection group ${groupName}:`, error);
		});
	},

	// For socket broadcasts - local state only, no API call
	setFrameSelection: (frame_selection) =>
		set({ frame_selection: frame_selection }),

	// For user-initiated changes - local state + API call
	updateFrameSelection: (frame_selection) => {
		const roomId = get().roomId;
		if (!roomId) return;

		// Optimistic update
		set({ frame_selection: frame_selection });

		// Update via REST API
		updateFrameSelectionAPI(roomId, frame_selection || []).catch((error) => {
			console.error("Failed to update frame selection:", error);
		});
	},
	setFrameSelectionEnabled: (enabled) =>
		set({ frameSelectionEnabled: enabled }),
	setBookmarks: (bookmarks) => set({ bookmarks: bookmarks }),
	setPlaying: (playing) => set({ playing: playing }),
	setChatOpen: (open) => set({ chatOpen: open }),
	setGeometries: (geometries) => set({ geometries: geometries }),
	setGeometrySchemas: (schemas) => set({ geometrySchemas: schemas }),
	setGeometryDefaults: (defaults) => set({ geometryDefaults: defaults }),
	updateGeometry: (
		key: string,
		geometry: any,
		source: "local" | "remote" = "remote",
	) =>
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
			const { [key]: removedSource, ...restSources } =
				state.geometryUpdateSources;
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
	setMode: (mode) => set({ mode }),
	setDrawingPointerPosition: (position) =>
		set({ drawingPointerPosition: position }),
	setDrawingIsValid: (isValid) => set({ drawingIsValid: isValid }),
	setEditingCombinedCentroid: (centroid) =>
		set({ editingCombinedCentroid: centroid }),
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
		const { lockRenewalIntervalId, lock, roomId, stopLockRenewal } = get();

		if (!lock) {
			console.warn("[startLockRenewal] No lock to renew");
			return;
		}

		// Clear existing interval if any
		if (lockRenewalIntervalId !== null) {
			stopLockRenewal();
		}

		// Use server-provided refresh interval (in seconds), default to 30s
		const refreshIntervalMs = (lock.refreshInterval || 30) * 1000;

		// Renew periodically to keep lock alive using idempotent acquire
		const intervalId = window.setInterval(async () => {
			const state = get();
			const currentLock = state.lock;

			if (!currentLock || !state.roomId) {
				// No lock anymore, stop renewal
				stopLockRenewal();
				return;
			}

			try {
				// Use idempotent acquireLock to refresh TTL
				const response = await acquireLock(
					state.roomId,
					currentLock.target,
					currentLock.message,
				);
				if (!response.success) {
					console.warn(
						`[startLockRenewal] Failed to refresh lock for ${currentLock.target}`,
					);
					// Lock was stolen or expired
					set({ lock: null });
					get().showSnackbar("Lock lost - returning to view mode", "warning");
					// Exit to view mode if active
					set({ mode: "view" });
					get().stopLockRenewal();
				}
			} catch (error) {
				console.error("[startLockRenewal] Error refreshing lock:", error);
				set({ lock: null });
				get().showSnackbar("Failed to refresh lock", "error");
				set({ mode: "view" });
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

	// Simplified Lock Manager - Target-based (aligns with server API)
	acquireLock: async (target, msg) => {
		const state = get();
		const { roomId, lock } = state;

		if (!roomId) {
			console.error(`[acquireLock] No roomId for target: ${target}`);
			return false;
		}

		// If we already have the lock for this target, just update the message
		if (lock && lock.target === target) {
			await get().updateLockMessage(target, msg);
			return true;
		}

		// Acquire new lock
		try {
			const response = await acquireLock(roomId, target, msg);

			if (!response.success || !response.lockToken) {
				return false;
			}

			// Store lock state
			set({
				lock: {
					target,
					token: response.lockToken,
					refreshInterval: response.refreshInterval || null,
					acquiredAt: Date.now(),
					message: msg,
				},
			});

			// Start lock renewal
			get().startLockRenewal();

			return true;
		} catch (error) {
			console.error(`[acquireLock] Error acquiring lock for ${target}:`, error);
			return false;
		}
	},

	releaseLock: async (target) => {
		const state = get();
		const { roomId, lock } = state;

		if (!lock) {
			return true;
		}

		// Check if we have the lock for this target
		if (lock.target !== target) {
			return false;
		}

		if (!roomId) {
			console.error(`[releaseLock] No roomId`);
			return false;
		}

		try {
			await releaseLock(roomId, target, lock.token);

			// Clear lock state
			set({ lock: null });
			get().stopLockRenewal();

			return true;
		} catch (error) {
			console.error(`[releaseLock] Error releasing lock for ${target}:`, error);
			return false;
		}
	},

	hasLock: (target) => {
		const { lock } = get();
		return lock !== null && lock.target === target;
	},

	updateLockMessage: async (target, msg) => {
		const state = get();
		const { roomId, lock } = state;

		if (!lock || lock.target !== target) {
			console.warn(`[updateLockMessage] Don't have lock for ${target}`);
			return false;
		}

		if (!roomId) {
			console.error(`[updateLockMessage] No roomId`);
			return false;
		}

		try {
			// Use idempotent acquireLock to update message and refresh TTL
			const response = await acquireLock(roomId, target, msg);

			if (response.success) {
				// Update local message
				set({
					lock: {
						...lock,
						message: msg,
					},
				});
				return true;
			}

			return false;
		} catch (error) {
			console.error(
				`[updateLockMessage] Error updating message for ${target}:`,
				error,
			);
			return false;
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

	toggleInfoBoxes: () =>
		set((state) => ({ showInfoBoxes: !state.showInfoBoxes })),
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

	enterEditingMode: async () => {
		const state = get();
		const { mode } = state;

		// Only enter from view mode
		if (mode !== "view") {
			console.warn(`Cannot enter editing mode from ${mode} mode`);
			return;
		}

		const acquired = await get().acquireLock(
			"trajectory:meta",
			"editing geometries",
		);
		if (acquired) {
			set({ mode: "editing" });
			state.showSnackbar("Entered editing mode", "success");
		}
	},

	exitEditingMode: async () => {
		const state = get();
		const { mode } = state;

		if (mode !== "editing") {
			console.warn(`Not in editing mode, current mode: ${mode}`);
			return;
		}

		// Save any pending frame edits before exiting
		await get().saveFrameEdits();

		// Clear selections when exiting editing mode
		const { selections } = state;
		for (const geometryKey of Object.keys(selections)) {
			get().updateSelectionForGeometry(geometryKey, []);
		}

		const released = await get().releaseLock("trajectory:meta");
		set({
			mode: "view",
			transformMode: "translate",
			editingSelectedAxis: null,
			editingFrameDataCount: 0, // Reset frame editing count
		}); // Reset when exiting
		if (released) {
			state.showSnackbar("Exited editing mode", "info");
		} else {
			// Even if release failed, exit editing mode
			console.warn("Failed to release lock when exiting editing mode");
		}
	},

	setTransformMode: (mode: "translate" | "rotate" | "scale") => {
		set({ transformMode: mode });
	},

	setEditingSelectedAxis: (axis: "x" | "y" | "z" | null) => {
		set({ editingSelectedAxis: axis });
	},

	cycleTransformMode: async () => {
		// Save any pending frame edits before switching modes
		await get().saveFrameEdits();

		const { transformMode, showSnackbar } = get();
		const modes: Array<"translate" | "rotate" | "scale"> = [
			"translate",
			"rotate",
			"scale",
		];
		const currentIndex = modes.indexOf(transformMode);
		const nextIndex = (currentIndex + 1) % modes.length;
		const nextMode = modes[nextIndex];
		set({ transformMode: nextMode });
		showSnackbar(`Transform mode: ${nextMode}`, "info");
	},

	enterDrawingMode: async (queryClient?: any) => {
		const state = get();
		const { geometries, mode, roomId } = state;

		// Only enter from view mode
		if (mode !== "view") {
			console.warn(`Cannot enter drawing mode from ${mode} mode`);
			return;
		}

		// Acquire lock for drawing mode
		const acquired = await get().acquireLock(
			"trajectory:meta",
			"drawing on curve",
		);
		if (!acquired) {
			return;
		}

		// Get all active curves
		const activeCurves = getActiveCurves(geometries);

		// Case 1: No curves exist - create one named "curve"
		if (activeCurves.length === 0) {
			if (!roomId) {
				console.warn("Cannot create curve: no room ID");
				return;
			}

			// Optimistic update - add geometry to store immediately
			set((state) => ({
				geometries: {
					...state.geometries,
					curve: {
						type: "Curve",
						data: { position: [] },
					},
				},
				activeCurveForDrawing: "curve",
				mode: "drawing",
			}));

			try {
				// Create geometry on backend - no lock needed, server handles via @check_lock
				await createGeometry(roomId, "curve", "Curve", {
					position: [],
				});

				// Fetch the validated geometry from backend to get defaults/transformations
				// (socket event won't reach us due to skip_sid)
				const response = await getGeometry(roomId, "curve");

				// Update with the backend-validated data
				set((state) => ({
					geometries: {
						...state.geometries,
						curve: response.geometry,
					},
				}));

				// Invalidate geometry list to refetch from server (this will update the sidebar)
				if (queryClient && roomId) {
					queryClient.invalidateQueries({
						queryKey: ["geometries", roomId, "list"],
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
						mode: "view",
					};
				});
			}
			return;
		}

		// Case 2: One curve exists - auto-select it and enter drawing mode
		if (activeCurves.length === 1) {
			const singleCurveKey = activeCurves[0];
			set({ activeCurveForDrawing: singleCurveKey, mode: "drawing" });
			return;
		}

		// Case 3: Multiple curves exist - preserve last selection if valid, otherwise select preferred
		const lastSelectedCurve = state.activeCurveForDrawing;
		const curveToActivate =
			lastSelectedCurve && activeCurves.includes(lastSelectedCurve)
				? lastSelectedCurve
				: selectPreferredCurve(activeCurves);
		set({ activeCurveForDrawing: curveToActivate, mode: "drawing" });
	},

	exitDrawingMode: async () => {
		const state = get();
		const { mode } = state;

		if (mode !== "drawing") {
			console.warn(`Not in drawing mode, current mode: ${mode}`);
			return;
		}

		const released = await get().releaseLock("trajectory:meta");
		// Keep activeCurveForDrawing to remember the last selected curve
		set({ mode: "view" });
		if (!released) {
			console.warn("Failed to release lock when exiting drawing mode");
		}
	},

	attachToCamera: (cameraKey) => {
		const { geometries, roomId, sessionId } = get();
		const camera = geometries[cameraKey];

		// Validate camera exists and is a Camera type
		if (!camera || camera.type !== "Camera") {
			console.warn(
				`Cannot attach to camera '${cameraKey}': not found or not a Camera`,
			);
			return;
		}

		set({ attachedCameraKey: cameraKey });

		// Sync to server so Python can read the current value
		if (roomId && sessionId) {
			const token = getToken();
			const headers: Record<string, string> = {
				"Content-Type": "application/json",
			};
			if (token) {
				headers["Authorization"] = `Bearer ${token}`;
			}
			fetch(`/api/rooms/${roomId}/sessions/${sessionId}/active-camera`, {
				method: "PUT",
				headers,
				body: JSON.stringify({ active_camera: cameraKey }),
			}).catch((error) => {
				console.error("Failed to sync active camera to server:", error);
			});
		}
	},

	setAttachedCameraKey: (cameraKey) => {
		// Direct setter for socket handler - no server sync to avoid loops
		set({ attachedCameraKey: cameraKey });
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

	incrementChatUnread: () =>
		set((state) => ({ chatUnreadCount: state.chatUnreadCount + 1 })),
	resetChatUnread: () => set({ chatUnreadCount: 0 }),
	setServerVersion: (version) => set({ serverVersion: version }),
	setGlobalSettings: (settings) => set({ globalSettings: settings }),
	showSnackbar: (message, severity = "info") =>
		set({ snackbar: { open: true, message, severity } }),
	hideSnackbar: () =>
		set((state) =>
			state.snackbar ? { snackbar: { ...state.snackbar, open: false } } : {},
		),
	setProgressTrackers: (progressTrackers) => set({ progressTrackers }),
	addProgressTracker: (progressId, description, progress, roomId) =>
		set((state) => ({
			progressTrackers: {
				...state.progressTrackers,
				[progressId]: { progressId, roomId, description, progress },
			},
		})),
	updateProgressTracker: (progressId, description, progress) =>
		set((state) => {
			const tracker = state.progressTrackers[progressId];
			if (!tracker) return {};
			return {
				progressTrackers: {
					...state.progressTrackers,
					[progressId]: {
						...tracker,
						...(description !== undefined && { description }),
						...(progress !== undefined && { progress }),
					},
				},
			};
		}),
	removeProgressTracker: (progressId) =>
		set((state) => {
			const { [progressId]: _, ...remainingTrackers } = state.progressTrackers;
			return { progressTrackers: remainingTrackers };
		}),

	// Frame editing actions
	setPendingFrameEdit: (frameId, key, data) => {
		const state = get();

		// If there are existing pending edits for a different frame, auto-save them first
		if (
			state.pendingFrameEdits &&
			state.pendingFrameEdits.frameId !== frameId
		) {
			// Schedule async save of old edits (fire and forget)
			const oldEdits = state.pendingFrameEdits;
			void (async () => {
				const { roomId, showSnackbar } = get();
				if (!roomId || Object.keys(oldEdits.keys).length === 0) return;

				try {
					const shapes = new Map<string, number[]>();
					for (const [k, d] of Object.entries(oldEdits.keys)) {
						if (k === "arrays.positions" && d instanceof Float32Array) {
							shapes.set(k, [d.length / 3, 3]);
						}
					}
					await partialUpdateFrame(
						roomId,
						oldEdits.frameId,
						oldEdits.keys,
						shapes,
					);
					showSnackbar(`Auto-saved frame ${oldEdits.frameId}`, "info");
				} catch (error) {
					console.error("[setPendingFrameEdit] Auto-save failed:", error);
					showSnackbar(
						`Failed to auto-save frame ${oldEdits.frameId}`,
						"warning",
					);
				}
			})();
		}

		set({
			pendingFrameEdits: {
				frameId,
				keys: {
					...(state.pendingFrameEdits?.frameId === frameId
						? state.pendingFrameEdits.keys
						: {}),
					[key]: data,
				},
			},
		});
	},

	clearPendingFrameEdits: () => set({ pendingFrameEdits: null }),

	saveFrameEdits: async () => {
		const state = get();
		const { roomId, pendingFrameEdits, showSnackbar } = state;

		if (!roomId) {
			console.warn("[saveFrameEdits] No roomId");
			return;
		}

		if (!pendingFrameEdits) {
			// Nothing to save
			return;
		}

		const { frameId, keys } = pendingFrameEdits;

		if (Object.keys(keys).length === 0) {
			// No keys to save
			return;
		}

		try {
			// Create shapes map for position data (N atoms x 3 dimensions)
			const shapes = new Map<string, number[]>();
			for (const [key, data] of Object.entries(keys)) {
				if (key === "arrays.positions" && data instanceof Float32Array) {
					// Positions are Nx3
					const numAtoms = data.length / 3;
					shapes.set(key, [numAtoms, 3]);
				}
			}

			// Call the API to save the partial update
			await partialUpdateFrame(roomId, frameId, keys, shapes);

			// Clear pending edits on success
			set({ pendingFrameEdits: null });

			showSnackbar(`Frame ${frameId} saved`, "success");
		} catch (error) {
			console.error("[saveFrameEdits] Error saving frame edits:", error);
			showSnackbar("Failed to save frame edits", "error");
		}
	},

	incrementEditingFrameDataCount: () =>
		set((state) => ({
			editingFrameDataCount: state.editingFrameDataCount + 1,
		})),

	decrementEditingFrameDataCount: () =>
		set((state) => ({
			editingFrameDataCount: Math.max(0, state.editingFrameDataCount - 1),
		})),

	registerLoadedDynamicPositions: (geometryKey, positionKey, positions) => {
		const newMap = new Map(get().loadedDynamicPositions);
		newMap.set(geometryKey, { positionKey, positions });
		set({ loadedDynamicPositions: newMap });
	},

	unregisterLoadedDynamicPositions: (geometryKey) => {
		const newMap = new Map(get().loadedDynamicPositions);
		newMap.delete(geometryKey);
		set({ loadedDynamicPositions: newMap });
	},
}));
