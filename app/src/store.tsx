import { create } from 'zustand';
import { produce } from 'immer'; // check how this works

// Define the shape of your frame data cache
interface FrameData {
  positions?: any | number[][];
}

interface AppState {
  // Connection & Room
  roomId: string | null;
  isConnected: boolean;
  currentFrame: number;
  frameData?: FrameData | null;
  isLoadingFrame: boolean;
  frameCount: number;

  // Actions (functions to modify the state)
  setConnected: (status: boolean, room: string) => void;
  setCurrentFrame: (frame: number) => void;
  setFrameData: (frameIndex: number, data: FrameData) => void;
  setLoadingFrame: (isLoading: boolean) => void;
  setFrameCount: (count: number) => void;
}

export const useAppStore = create<AppState>((set) => ({
  // Initial State
  roomId: null,
  isConnected: false,
  currentFrame: 0,
  frameData: null,
  isLoadingFrame: false,
  frameCount: 0,

  // Actions
  setConnected: (status, room) => set({ isConnected: status, roomId: room }),
  setCurrentFrame: (frame) => set({ currentFrame: frame }),
  setFrameData: (frameIndex, data) => set(produce(draft => {
    draft.frameData = { ...draft.frameData, [frameIndex]: data };
  })),
  setLoadingFrame: (isLoading) => set({ isLoadingFrame: isLoading }),
  setFrameCount: (count) => set({ frameCount: count }),
}));
