import type { StateCreator } from "zustand";
import {
	deleteBookmark as deleteBookmarkAPI,
	setBookmark as setBookmarkAPI,
	updateFrameSelection as updateFrameSelectionAPI,
} from "../../myapi/client";
import type { AppState } from "../../store";

export interface PlaybackSlice {
	currentFrame: number;
	frameCount: number;
	skipFrames: number;
	frame_selection: number[] | null;
	frameSelectionEnabled: boolean;
	playing: boolean;
	bookmarks: Record<number, string> | null;
	fps: number | null;
	frameLoadTime: number | null;
	lastFrameChangeTime: number | null;
	synchronizedMode: boolean;

	setCurrentFrame: (frame: number) => void;
	setFrameCount: (count: number) => void;
	setSkipFrames: (skip: number) => void;
	setFrameSelection: (selection: number[] | null) => void;
	updateFrameSelection: (selection: number[] | null) => void;
	setFrameSelectionEnabled: (enabled: boolean) => void;
	setPlaying: (playing: boolean) => void;
	setBookmarks: (bookmark: Record<number, string> | null) => void;
	addBookmark: (frame: number, label?: string) => void;
	deleteBookmark: (frame: number) => void;
	setFps: (fps: number | null) => void;
	setFrameLoadTime: (time: number | null) => void;
	setLastFrameChangeTime: (time: number | null) => void;
	setSynchronizedMode: (enabled: boolean) => void;
}

export const createPlaybackSlice: StateCreator<
	AppState,
	[],
	[],
	PlaybackSlice
> = (set, get) => ({
	currentFrame: 0,
	frameCount: 0,
	skipFrames: 1,
	frame_selection: null,
	frameSelectionEnabled: false,
	playing: false,
	bookmarks: null,
	fps: null,
	frameLoadTime: null,
	lastFrameChangeTime: null,
	synchronizedMode: true,

	setCurrentFrame: (frame) => {
		if (typeof frame !== "number" || !Number.isFinite(frame)) {
			console.error("Invalid frame value:", frame);
			return;
		}
		set({ currentFrame: frame });
	},

	setFrameCount: (count) =>
		set((state) => {
			const newCurrentFrame =
				state.currentFrame >= count ? 0 : state.currentFrame;
			return { frameCount: count, currentFrame: newCurrentFrame };
		}),

	setSkipFrames: (skip) => set({ skipFrames: skip }),

	setFrameSelection: (frame_selection) => set({ frame_selection }),

	updateFrameSelection: (frame_selection) => {
		const roomId = get().roomId;
		if (!roomId) return;
		set({ frame_selection });
		updateFrameSelectionAPI(roomId, frame_selection || []).catch((error) => {
			console.error("Failed to update frame selection:", error);
		});
	},

	setFrameSelectionEnabled: (enabled) =>
		set({ frameSelectionEnabled: enabled }),

	setPlaying: (playing) => set({ playing }),

	setBookmarks: (bookmarks) => set({ bookmarks }),

	addBookmark: (frame, label) => {
		const roomId = get().roomId;
		if (!roomId) return;
		const bookmarkLabel = label || `Bookmark ${frame}`;
		set((state) => ({
			bookmarks: { ...(state.bookmarks || {}), [frame]: bookmarkLabel },
		}));
		setBookmarkAPI(roomId, frame, bookmarkLabel).catch((error) => {
			console.error(`Failed to set bookmark at frame ${frame}:`, error);
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
		const oldLabel = bookmarks[frame];
		set((state) => {
			if (!state.bookmarks) return state;
			const { [frame]: removed, ...rest } = state.bookmarks;
			return { bookmarks: rest };
		});
		deleteBookmarkAPI(roomId, frame).catch((error) => {
			console.error(`Failed to delete bookmark at frame ${frame}:`, error);
			set((state) => ({
				bookmarks: { ...(state.bookmarks || {}), [frame]: oldLabel },
			}));
		});
	},

	setFps: (fps) => set({ fps }),
	setFrameLoadTime: (time) => set({ frameLoadTime: time }),
	setLastFrameChangeTime: (time) => set({ lastFrameChangeTime: time }),
	setSynchronizedMode: (enabled) => set({ synchronizedMode: enabled }),
});
