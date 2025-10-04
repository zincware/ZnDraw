import { create } from 'zustand';

// Defines the state for a single window
export interface WindowState {
  x: number;
  y: number;
  width: number;
  height: number;
  zIndex: number;
}

// Defines the entire store's state and actions
interface WindowManagerStore {
  openWindows: Record<string, WindowState>; // A dictionary of figureKey -> WindowState
  _zIndexCounter: number; // Private-like counter to ensure new windows are on top
  openWindow: (key: string) => void;
  closeWindow: (key: string) => void;
  updateWindowState: (key: string, updates: Partial<WindowState>) => void;
  bringToFront: (key: string) => void;
}

export const useWindowManagerStore = create<WindowManagerStore>((set, get) => ({
  openWindows: {},
  _zIndexCounter: 100, // Start z-index from a base number

  /**
   * Opens a new window for a given figure key or brings it to the front if already open.
   */
  openWindow: (key) => {
    const { openWindows, bringToFront } = get();
    if (openWindows[key]) {
      bringToFront(key); // If already open, just bring it to the front
      return;
    }

    const newZIndex = get()._zIndexCounter + 1;
    // Cascade new windows for better visibility
    const openCount = Object.keys(openWindows).length;
    
    set((state) => ({
      openWindows: {
        ...state.openWindows,
        [key]: {
          x: 50 + openCount * 20,
          y: 100 + openCount * 20,
          width: 400,
          height: 300,
          zIndex: newZIndex,
        },
      },
      _zIndexCounter: newZIndex,
    }));
  },

  /**
   * Closes a window by removing it from the state.
   */
  closeWindow: (key) => {
    set((state) => {
      const { [key]: _, ...remainingWindows } = state.openWindows;
      return { openWindows: remainingWindows };
    });
  },

  /**
   * Updates a window's state (position, size) after dragging or resizing.
   */
  updateWindowState: (key, updates) => {
    set((state) => ({
      openWindows: {
        ...state.openWindows,
        [key]: { ...state.openWindows[key], ...updates },
      },
    }));
  },

  /**
   * Brings a specific window to the top by giving it the highest z-index.
   */
  bringToFront: (key) => {
    const currentWindow = get().openWindows[key];
    if (!currentWindow) return;

    const newZIndex = get()._zIndexCounter + 1;
    
    // Only update if it's not already on top
    if (currentWindow.zIndex !== newZIndex -1) {
        set((state) => ({
          openWindows: {
            ...state.openWindows,
            [key]: { ...currentWindow, zIndex: newZIndex },
          },
          _zIndexCounter: newZIndex,
        }));
    }
  },
}));