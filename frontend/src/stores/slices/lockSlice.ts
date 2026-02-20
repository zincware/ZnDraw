import type { StateCreator } from "zustand";
import { acquireEditLock, releaseEditLock } from "../../myapi/client";
import type { AppState } from "../../store";

export interface LockSlice {
	superuserLock: boolean;
	userLock: string | null;
	userLockMessage: string | null;
	lockRenewalIntervalId: number | null;

	setSuperuserLock: (locked: boolean) => void;
	setUserLock: (email: string | null, message?: string | null) => void;
	startLockRenewal: () => void;
	stopLockRenewal: () => void;
	acquireLock: (msg: string) => Promise<boolean>;
	releaseLock: () => Promise<boolean>;
	hasLock: () => boolean;
}

export const createLockSlice: StateCreator<AppState, [], [], LockSlice> = (
	set,
	get,
) => ({
	superuserLock: false,
	userLock: null,
	userLockMessage: null,
	lockRenewalIntervalId: null,

	setSuperuserLock: (locked) => set({ superuserLock: locked }),

	setUserLock: (email, message = null) =>
		set({ userLock: email, userLockMessage: message ?? null }),

	startLockRenewal: () => {
		const { lockRenewalIntervalId, userLock, stopLockRenewal } = get();
		if (!userLock) {
			console.warn("[startLockRenewal] No lock to renew");
			return;
		}
		if (lockRenewalIntervalId !== null) {
			stopLockRenewal();
		}
		const refreshIntervalMs = 5000;
		const intervalId = window.setInterval(async () => {
			const state = get();
			const currentMessage = state.userLockMessage;
			if (!state.userLock || !state.roomId) {
				stopLockRenewal();
				return;
			}
			try {
				await acquireEditLock(state.roomId, currentMessage ?? undefined);
			} catch (error) {
				console.error("[startLockRenewal] Error refreshing lock:", error);
				set({ userLock: null, userLockMessage: null });
				get().showSnackbar("Lock lost - returning to view mode", "warning");
				set({ mode: "view" });
				get().stopLockRenewal();
			}
		}, refreshIntervalMs);
		set({ lockRenewalIntervalId: intervalId });
	},

	stopLockRenewal: () => {
		const { lockRenewalIntervalId } = get();
		if (lockRenewalIntervalId !== null) {
			window.clearInterval(lockRenewalIntervalId);
			set({ lockRenewalIntervalId: null });
		}
	},

	acquireLock: async (msg) => {
		const { roomId, userLock, user } = get();
		if (!roomId) {
			console.error("[acquireLock] No roomId");
			return false;
		}
		if (userLock) {
			try {
				await acquireEditLock(roomId, msg);
				set({ userLockMessage: msg });
				return true;
			} catch (error) {
				console.error("[acquireLock] Error refreshing lock:", error);
				return false;
			}
		}
		try {
			const response = await acquireEditLock(roomId, msg);
			if (!response.locked) return false;
			const currentUserEmail = user?.email ?? null;
			set({ userLock: currentUserEmail, userLockMessage: msg });
			get().startLockRenewal();
			return true;
		} catch (error) {
			console.error("[acquireLock] Error acquiring lock:", error);
			return false;
		}
	},

	releaseLock: async () => {
		const { roomId, userLock } = get();
		if (!userLock) return true;
		if (!roomId) {
			console.error("[releaseLock] No roomId");
			return false;
		}
		try {
			await releaseEditLock(roomId);
			set({ userLock: null, userLockMessage: null });
			get().stopLockRenewal();
			return true;
		} catch (error) {
			console.error("[releaseLock] Error releasing lock:", error);
			return false;
		}
	},

	hasLock: () => get().userLock !== null,
});
