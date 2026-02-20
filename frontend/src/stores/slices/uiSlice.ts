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

	toggleInfoBoxes: () =>
		set((state) => ({ showInfoBoxes: !state.showInfoBoxes })),

	setHoveredFrame: (frame) => set({ hoveredFrame: frame }),

	setScreenshotCapture: (fn) => set({ screenshotCapture: fn }),
	setPathtracerCapture: (fn) => set({ pathtracerCapture: fn }),

	requestPathtracingUpdate: () => set({ pathtracingNeedsUpdate: true }),
	clearPathtracingUpdate: () => set({ pathtracingNeedsUpdate: false }),
});
