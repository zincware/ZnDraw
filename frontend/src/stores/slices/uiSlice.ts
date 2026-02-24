import type { StateCreator } from "zustand";
import type { AppState, Progress } from "../../store";

export interface UISlice {
	chatOpen: boolean;
	chatUnreadCount: number;
	typingUsers: Set<string>;
	snackbar: {
		open: boolean;
		message: string;
		severity: "success" | "info" | "warning" | "error";
	} | null;
	progressTrackers: Record<string, Progress>;
	showInfoBoxes: boolean;
	hoveredFrame: number | null;
	screenshotCapture: (() => Promise<Blob>) | null;
	pathtracerCapture: (() => Promise<Blob>) | null;
	pathtracingNeedsUpdate: boolean;

	setChatOpen: (open: boolean) => void;
	incrementChatUnread: () => void;
	resetChatUnread: () => void;
	addTypingUser: (email: string) => void;
	removeTypingUser: (email: string) => void;
	showSnackbar: (
		message: string,
		severity?: "success" | "info" | "warning" | "error",
	) => void;
	hideSnackbar: () => void;
	setProgressTrackers: (trackers: Record<string, Progress>) => void;
	addProgressTracker: (tracker: Progress) => void;
	updateProgressTracker: (update: Partial<Progress> & { progress_id: string }) => void;
	removeProgressTracker: (progressId: string) => void;
	toggleInfoBoxes: () => void;
	setHoveredFrame: (frame: number | null) => void;
	setScreenshotCapture: (fn: (() => Promise<Blob>) | null) => void;
	setPathtracerCapture: (fn: (() => Promise<Blob>) | null) => void;
	requestPathtracingUpdate: () => void;
	clearPathtracingUpdate: () => void;
}

export const createUISlice: StateCreator<AppState, [], [], UISlice> = (
	set,
) => ({
	chatOpen: false,
	chatUnreadCount: 0,
	typingUsers: new Set(),
	snackbar: null,
	progressTrackers: {},
	showInfoBoxes: false,
	hoveredFrame: null,
	screenshotCapture: null,
	pathtracerCapture: null,
	pathtracingNeedsUpdate: false,

	setChatOpen: (open) => set({ chatOpen: open }),

	incrementChatUnread: () =>
		set((state) => ({ chatUnreadCount: state.chatUnreadCount + 1 })),

	resetChatUnread: () => set({ chatUnreadCount: 0 }),

	addTypingUser: (email) =>
		set((state) => {
			const next = new Set(state.typingUsers);
			next.add(email);
			return { typingUsers: next };
		}),

	removeTypingUser: (email) =>
		set((state) => {
			const next = new Set(state.typingUsers);
			next.delete(email);
			return { typingUsers: next };
		}),

	showSnackbar: (message, severity = "info") =>
		set({ snackbar: { open: true, message, severity } }),

	hideSnackbar: () =>
		set((state) =>
			state.snackbar ? { snackbar: { ...state.snackbar, open: false } } : {},
		),

	setProgressTrackers: (progressTrackers) => set({ progressTrackers }),

	addProgressTracker: (tracker) =>
		set((state) => ({
			progressTrackers: {
				...state.progressTrackers,
				[tracker.progress_id]: tracker,
			},
		})),

	updateProgressTracker: (update) =>
		set((state) => {
			const tracker = state.progressTrackers[update.progress_id];
			if (!tracker) return {};
			return {
				progressTrackers: {
					...state.progressTrackers,
					[update.progress_id]: { ...tracker, ...update },
				},
			};
		}),

	removeProgressTracker: (progressId) =>
		set((state) => {
			const { [progressId]: _, ...remainingTrackers } = state.progressTrackers;
			return { progressTrackers: remainingTrackers };
		}),

	toggleInfoBoxes: () =>
		set((state) => ({ showInfoBoxes: !state.showInfoBoxes })),

	setHoveredFrame: (frame) => set({ hoveredFrame: frame }),

	setScreenshotCapture: (fn) => set({ screenshotCapture: fn }),
	setPathtracerCapture: (fn) => set({ pathtracerCapture: fn }),

	requestPathtracingUpdate: () => set({ pathtracingNeedsUpdate: true }),
	clearPathtracingUpdate: () => set({ pathtracingNeedsUpdate: false }),
});
