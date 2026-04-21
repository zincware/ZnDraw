import type { QueryClient } from "@tanstack/react-query";
import type { GlobalSettings } from "../../myapi/client";
import type { InitializationError, Progress } from "../../store";
import type { UserInfo } from "../../utils/auth";

/**
 * Flat dependency context passed to all handler factories.
 *
 * Fields mirror the selectors extracted at the top of useSocketManager.
 * Values that change during the effect lifetime (e.g. `playing`, active panel)
 * are NOT included -- handlers read those via `useAppStore.getState()`.
 */
export interface HandlerContext {
	// Identity / routing
	roomId: string | undefined;
	appStoreRoomId: string | null;
	isOverview: boolean;
	isCancelled: () => boolean;

	// Query client
	queryClient: QueryClient;

	// Connection setters
	setConnected: (status: boolean) => void;
	setInitializationError: (error: InitializationError | null) => void;
	setUser: (user: UserInfo) => void;
	setSessionId: (sessionId: string | null) => void;
	setCameraKey: (cameraKey: string | null) => void;
	setServerVersion: (version: string | null) => void;
	setGlobalSettings: (settings: GlobalSettings | null) => void;

	// Playback setters
	setFrameCount: (count: number) => void;
	setCurrentFrame: (frame: number) => void;
	setFrameSelection: (selection: number[] | null) => void;
	setBookmarks: (bookmarks: Record<number, string> | null) => void;

	// Scene/geometry setters
	setSelections: (selections: Record<string, number[]>) => void;
	setSelectionGroups: (
		groups: Record<string, Record<string, number[]>>,
	) => void;
	setGeometries: (geometries: Record<string, any>) => void;
	setGeometrySchemas: (schemas: Record<string, any>) => void;
	setGeometryDefaults: (defaults: Record<string, any>) => void;
	updateGeometry: (
		key: string,
		geometry: any,
		source?: "local" | "remote",
	) => void;
	removeGeometry: (key: string) => void;
	setActiveCurveForDrawing: (key: string | null) => void;

	// Lock setters
	setSuperuserLock: (locked: boolean) => void;
	setUserLock: (email: string | null, message?: string | null) => void;

	// Progress setters
	setProgressTrackers: (trackers: Record<string, Progress>) => void;
	addProgressTracker: (tracker: Progress) => void;
	updateProgressTracker: (
		update: Partial<Progress> & { progress_id: string },
	) => void;
	removeProgressTracker: (progressId: string) => void;
}
