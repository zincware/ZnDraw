import { create } from "zustand";
import type { ConnectionSlice } from "./stores/slices/connectionSlice";
import { createConnectionSlice } from "./stores/slices/connectionSlice";
import type { LockSlice } from "./stores/slices/lockSlice";
import { createLockSlice } from "./stores/slices/lockSlice";
import type { PlaybackSlice } from "./stores/slices/playbackSlice";
import { createPlaybackSlice } from "./stores/slices/playbackSlice";
import type { SceneSlice } from "./stores/slices/sceneSlice";
import { createSceneSlice } from "./stores/slices/sceneSlice";
import type { UISlice } from "./stores/slices/uiSlice";
import { createUISlice } from "./stores/slices/uiSlice";

// Re-export helpers that live in slices
export {
	getActiveCurves,
	selectPreferredCurve,
} from "./stores/slices/sceneSlice";

/**
 * Progress tracking state — mirrors server ProgressResponse (snake_case).
 */
export interface Progress {
	progress_id: string;
	description: string;
	n: number;
	total: number | null;
	elapsed: number;
	unit: string;
}

/**
 * Error state for room initialization failures
 */
export interface InitializationError {
	message: string;
	details?: string;
}

export type AppState = ConnectionSlice &
	PlaybackSlice &
	SceneSlice &
	LockSlice &
	UISlice;

// --- Selectors ---

/**
 * Selector: true when the current user cannot write to the room.
 * Combines admin lock + edit lock, with superuser bypass.
 */
export const selectIsRoomReadOnly = (state: AppState): boolean => {
	const isSuperuser = state.user?.is_superuser ?? false;
	if (isSuperuser) return false;
	if (state.superuserLock) return true;
	if (state.userLock) {
		const userEmail = state.user?.email ?? null;
		return state.userLock !== userEmail;
	}
	return false;
};

/**
 * Derived selector: returns the name of the selection group that exactly matches
 * the current selections, or null if none match.
 */
export const selectActiveSelectionGroup = (state: AppState): string | null => {
	const { selections, selectionGroups } = state;
	for (const [name, group] of Object.entries(selectionGroups)) {
		if (selectionsEqual(selections, group)) return name;
	}
	return null;
};

/** Compare two selection dicts for equality (empty/missing keys are equivalent). */
function selectionsEqual(
	a: Record<string, number[]>,
	b: Record<string, number[]>,
): boolean {
	const allKeys = new Set([...Object.keys(a), ...Object.keys(b)]);
	for (const key of allKeys) {
		const aArr = a[key] ?? [];
		const bArr = b[key] ?? [];
		if (aArr.length !== bArr.length) return false;
		const aSorted = [...aArr].sort((x, y) => x - y);
		const bSorted = [...bArr].sort((x, y) => x - y);
		for (let i = 0; i < aSorted.length; i++) {
			if (aSorted[i] !== bSorted[i]) return false;
		}
	}
	return true;
}

// --- Store ---

export const useAppStore = create<AppState>((...a) => ({
	...createConnectionSlice(...a),
	...createPlaybackSlice(...a),
	...createSceneSlice(...a),
	...createLockSlice(...a),
	...createUISlice(...a),
}));
