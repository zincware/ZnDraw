import { create } from 'zustand';

interface AppState {
  // Connection & Room
  roomId: string | null;
  userId: string | null;
  isConnected: boolean;
  isLoading: boolean;
  currentFrame: number;
  frameCount: number;
  skipFrames: number;

  // Actions (functions to modify the state)
  setConnected: (status: boolean, room: string, user: string) => void;
  setCurrentFrame: (frame: number) => void;
  setFrameCount: (count: number) => void;
  setLoading: (loading: boolean) => void;
  setSkipFrames: (skip: number) => void;
}

export const useAppStore = create<AppState>((set) => ({
  // Initial State
  roomId: null,
  userId: null,
  isConnected: false,
  currentFrame: 0,
  frameCount: 0,
  isLoading: false,
  skipFrames: 1,

  // Actions
  setConnected: (status, room, user) => set({ isConnected: status, roomId: room, userId: user }),
  setCurrentFrame: (frame) => set({ currentFrame: frame }),
  setFrameCount: (count) => set((state) => {
    // If currentFrame is beyond the new frameCount, reset to 0
    const newCurrentFrame = state.currentFrame >= count ? 0 : state.currentFrame;
    return { frameCount: count, currentFrame: newCurrentFrame };
  }),
  setLoading: (loading) => set({ isLoading: loading }),
  setSkipFrames: (skip) => set({ skipFrames: skip }),
}));
