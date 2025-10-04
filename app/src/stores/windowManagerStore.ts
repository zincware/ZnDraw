import { create } from 'zustand';

// A window now has its own stable ID and a figureKey property for its content
export interface WindowInstance {
  id: string;
  figureKey: string; 
  x: number;
  y: number;
  width: number;
  height: number;
  zIndex: number;
}

interface WindowManagerStore {
  openWindows: Record<string, WindowInstance>; // The key is now a unique windowId
  _zIndexCounter: number;
  openWindow: (figureKey: string) => void;
  closeWindow: (windowId: string) => void;
  changeFigureInWindow: (windowId: string, newFigureKey: string) => void; // New action
  updateWindowState: (windowId: string, updates: Partial<Omit<WindowInstance, 'id' | 'figureKey'>>) => void;
  bringToFront: (windowId: string) => void;
}

export const useWindowManagerStore = create<WindowManagerStore>((set, get) => ({
  openWindows: {},
  _zIndexCounter: 100,

  /**
   * Opens a new window displaying the given figure.
   */
  openWindow: (figureKey) => {
    const newZIndex = get()._zIndexCounter + 1;
    const openCount = Object.keys(get().openWindows).length;
    const windowId = `window-${Date.now()}`; // Generate a unique, stable ID

    set((state) => ({
      openWindows: {
        ...state.openWindows,
        [windowId]: {
          id: windowId,
          figureKey, // The figure to display initially
          x: 50 + openCount * 20,
          y: 100 + openCount * 20,
          width: 500,
          height: 400,
          zIndex: newZIndex,
        },
      },
      _zIndexCounter: newZIndex,
    }));
  },

  /**
   * Closes a window using its unique ID.
   */
  closeWindow: (windowId) => {
    set((state) => {
      const { [windowId]: _, ...remainingWindows } = state.openWindows;
      return { openWindows: remainingWindows };
    });
  },
  
  /**
   * Changes the figure displayed within an existing window.
   */
  changeFigureInWindow: (windowId, newFigureKey) => {
    set((state) => {
      if (!state.openWindows[windowId]) return state; // Safety check
      return {
        openWindows: {
          ...state.openWindows,
          [windowId]: {
            ...state.openWindows[windowId],
            figureKey: newFigureKey,
          },
        },
      };
    });
  },

  /**
   * Updates a window's state (position, size) after user interaction.
   */
  updateWindowState: (windowId, updates) => {
    set((state) => ({
      openWindows: {
        ...state.openWindows,
        [windowId]: { ...state.openWindows[windowId], ...updates },
      },
    }));
  },

  /**
   * Brings a specific window to the top.
   */
  bringToFront: (windowId) => {
    const currentWindow = get().openWindows[windowId];
    if (!currentWindow) return;

    const newZIndex = get()._zIndexCounter + 1;
    
    if (currentWindow.zIndex !== newZIndex - 1) {
        set((state) => ({
          openWindows: {
            ...state.openWindows,
            [windowId]: { ...currentWindow, zIndex: newZIndex },
          },
          _zIndexCounter: newZIndex,
        }));
    }
  },
}));