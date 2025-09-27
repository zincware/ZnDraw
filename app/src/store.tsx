import { create } from 'zustand';

interface AppState {
  // Connection & Room
  roomId: string | null;
  isConnected: boolean;
  isLoading: boolean;
  currentFrame: number;
  frameCount: number;
  skipFrames: number;
  isPresenter: boolean; // Is the current client the presenter?
  presenterSid: string | null;

  // Actions (functions to modify the state)
  setConnected: (status: boolean, room: string) => void;
  setCurrentFrame: (frame: number) => void;
  setFrameCount: (count: number) => void;
  setLoading: (loading: boolean) => void;
  setSkipFrames: (skip: number) => void;
  setPresenter: (status: boolean) => void;
  setPresenterSid: (sid: string | null) => void;
}

export const useAppStore = create<AppState>((set) => ({
  // Initial State
  roomId: null,
  isConnected: false,
  currentFrame: 0,
  frameCount: 0,
  isLoading: false,
  skipFrames: 1,
  isPresenter: false,
  presenterSid: null,

  // Actions
  setConnected: (status, room) => set({ isConnected: status, roomId: room }),
  setCurrentFrame: (frame) => set({ currentFrame: frame }),
  setFrameCount: (count) => set({ frameCount: count }),
  setLoading: (loading) => set({ isLoading: loading }),
  setSkipFrames: (skip) => set({ skipFrames: skip }),
  setPresenter: (status) => set({ isPresenter: status }),
  setPresenterSid: (sid) => set({ presenterSid: sid }),
}));
