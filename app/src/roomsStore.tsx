import { create } from "zustand";
import { listRooms, getRoom, Room } from "./myapi/client";

interface RoomsState {
  // State - store both Map for fast lookups and array for rendering
  roomsMap: Map<string, Room>;
  roomsArray: Room[];
  loading: boolean;
  error: string | null;

  // Actions
  fetchRooms: () => Promise<void>;
  setRoom: (roomId: string, room: Room) => void;
  updateRoom: (roomId: string, updates: Partial<Room>) => void;
  removeRoom: (roomId: string) => void;
  getRoom: (roomId: string) => Room | undefined;
}

// Helper to convert Map to sorted array
function mapToSortedArray(map: Map<string, Room>): Room[] {
  return Array.from(map.values()).sort((a, b) => a.id.localeCompare(b.id));
}

export const useRoomsStore = create<RoomsState>((set, get) => ({
  // Initial state
  roomsMap: new Map(),
  roomsArray: [],
  loading: false,
  error: null,

  // Fetch all rooms from API (initial load)
  fetchRooms: async () => {
    set({ loading: true, error: null });
    try {
      const roomsList = await listRooms();
      const roomsMap = new Map<string, Room>();
      roomsList.forEach((room) => {
        roomsMap.set(room.id, room);
      });
      const roomsArray = mapToSortedArray(roomsMap);
      set({ roomsMap, roomsArray, loading: false });
    } catch (error) {
      console.error("Failed to fetch rooms:", error);
      set({ error: "Failed to fetch rooms", loading: false });
    }
  },

  // Set a complete room object
  setRoom: (roomId: string, room: Room) => {
    set((state) => {
      const newRoomsMap = new Map(state.roomsMap);
      newRoomsMap.set(roomId, room);
      const newRoomsArray = mapToSortedArray(newRoomsMap);
      return { roomsMap: newRoomsMap, roomsArray: newRoomsArray };
    });
  },

  // Update specific fields of a room
  updateRoom: (roomId: string, updates: Partial<Room>) => {
    set((state) => {
      const existingRoom = state.roomsMap.get(roomId);
      if (!existingRoom) {
        // Room doesn't exist yet, create it with updates
        const newRoom: Room = {
          id: roomId,
          frameCount: 0,
          locked: false,
          hidden: false,
          isDefault: false,
          ...updates,
        };
        const newRoomsMap = new Map(state.roomsMap);
        newRoomsMap.set(roomId, newRoom);
        const newRoomsArray = mapToSortedArray(newRoomsMap);
        return { roomsMap: newRoomsMap, roomsArray: newRoomsArray };
      }

      // Update existing room
      const updatedRoom = { ...existingRoom, ...updates };
      const newRoomsMap = new Map(state.roomsMap);
      newRoomsMap.set(roomId, updatedRoom);
      const newRoomsArray = mapToSortedArray(newRoomsMap);
      return { roomsMap: newRoomsMap, roomsArray: newRoomsArray };
    });
  },

  // Remove a room
  removeRoom: (roomId: string) => {
    set((state) => {
      const newRoomsMap = new Map(state.roomsMap);
      newRoomsMap.delete(roomId);
      const newRoomsArray = mapToSortedArray(newRoomsMap);
      return { roomsMap: newRoomsMap, roomsArray: newRoomsArray };
    });
  },

  // Get a specific room
  getRoom: (roomId: string) => {
    return get().roomsMap.get(roomId);
  },
}));
