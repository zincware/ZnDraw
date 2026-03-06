import type { StateCreator } from "zustand";
import type { GlobalSettings } from "../../myapi/client";
import type { UserInfo } from "../../utils/auth";
import type { AppState, InitializationError } from "../../store";

export interface ConnectionSlice {
	roomId: string | null;
	user: UserInfo | null;
	sessionId: string | null;
	cameraKey: string | null;
	isConnected: boolean;
	isLoading: boolean;
	initializationError: InitializationError | null;
	serverVersion: string | null;
	globalSettings: GlobalSettings | null;

	setRoomId: (roomId: string) => void;
	setUser: (user: UserInfo) => void;
	setSessionId: (sessionId: string | null) => void;
	setCameraKey: (cameraKey: string | null) => void;
	setConnected: (status: boolean) => void;
	setLoading: (loading: boolean) => void;
	setInitializationError: (error: InitializationError | null) => void;
	setServerVersion: (version: string | null) => void;
	setGlobalSettings: (settings: GlobalSettings | null) => void;
}

export const createConnectionSlice: StateCreator<
	AppState,
	[],
	[],
	ConnectionSlice
> = (set) => ({
	roomId: null,
	user: null,
	sessionId: null,
	cameraKey: null,
	isConnected: false,
	isLoading: false,
	initializationError: null,
	serverVersion: null,
	globalSettings: null,

	setRoomId: (roomId) => set({ roomId }),
	setUser: (user) => set({ user }),
	setSessionId: (sessionId) => set({ sessionId }),
	setCameraKey: (cameraKey) => set({ cameraKey }),
	setConnected: (status) => set({ isConnected: status }),
	setLoading: (loading) => set({ isLoading: loading }),
	setInitializationError: (error) => set({ initializationError: error }),
	setServerVersion: (version) => set({ serverVersion: version }),
	setGlobalSettings: (settings) => set({ globalSettings: settings }),
});
