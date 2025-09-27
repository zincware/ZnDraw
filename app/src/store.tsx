import { create } from 'zustand';
import { produce } from 'immer'; // check how this works
import * as THREE from 'three';

// Define the shape of your frame data cache
interface FrameData {
  positions?: THREE.Vector3[];
}

interface AppState {
  // Connection & Room
  roomId: string | null;
  isConnected: boolean;
  currentFrame: number;
  frameCount: number;

  // Actions (functions to modify the state)
  setConnected: (status: boolean, room: string) => void;
  setCurrentFrame: (frame: number) => void;
  setFrameCount: (count: number) => void;
}

export const useAppStore = create<AppState>((set) => ({
  // Initial State
  roomId: null,
  isConnected: false,
  currentFrame: 0,
  frameCount: 0,

  // Actions
  setConnected: (status, room) => set({ isConnected: status, roomId: room }),
  setCurrentFrame: (frame) => set({ currentFrame: frame }),
  setFrameCount: (count) => set({ frameCount: count }),
}));
