import type { StateCreator } from "zustand";
import {
	acquireEditLock,
	getEditLockStatus,
	releaseEditLock,
} from "../../myapi/client";
import type { AppState } from "../../store";

export interface LockSlice {
	superuserLock: boolean;
	userLock: string | null;
	userLockMessage: string | null;
	lockToken: string | null;
	lockRenewalIntervalId: number | null;
	lockExpiryTimerId: number | null;

	setSuperuserLock: (locked: boolean) => void;
	setUserLock: (email: string | null, message?: string | null) => void;
	startLockRenewal: () => void;
	stopLockRenewal: () => void;
	acquireLock: (msg: string) => Promise<boolean>;
	releaseLock: () => Promise<boolean>;
	hasLock: () => boolean;
	startLockExpiryTimer: (ttl: number) => void;
	stopLockExpiryTimer: () => void;
}

export const createLockSlice: StateCreator<AppState, [], [], LockSlice> = (
	set,
	get,
) => ({
	superuserLock: false,
	userLock: null,
	userLockMessage: null,
	lockToken: null,
	lockRenewalIntervalId: null,
	lockExpiryTimerId: null,

	setSuperuserLock: (locked) => set({ superuserLock: locked }),

	setUserLock: (email, message = null) =>
		set({ userLock: email, userLockMessage: message ?? null }),

	startLockRenewal: () => {
		const { lockRenewalIntervalId, lockToken, stopLockRenewal } = get();
		if (!lockToken) {
			console.warn("[startLockRenewal] No lock token to renew");
			return;
		}
		if (lockRenewalIntervalId !== null) {
			stopLockRenewal();
		}
		const refreshIntervalMs = 5000;
		const intervalId = window.setInterval(async () => {
			const state = get();
			if (!state.lockToken || !state.roomId) {
				stopLockRenewal();
				return;
			}
			try {
				await acquireEditLock(
					state.roomId,
					state.userLockMessage ?? undefined,
					state.lockToken,
				);
			} catch (error: any) {
				console.error("[startLockRenewal] Error refreshing lock:", error);
				set({
					userLock: null,
					userLockMessage: null,
					lockToken: null,
				});
				const is409 = error?.response?.status === 409;
				get().showSnackbar(
					is409
						? "Lock expired - returning to view mode"
						: "Lock lost - returning to view mode",
					"warning",
				);
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
		const { roomId, lockToken, user } = get();
		if (!roomId) {
			console.error("[acquireLock] No roomId");
			return false;
		}
		if (lockToken) {
			// Already holding — refresh
			try {
				await acquireEditLock(roomId, msg, lockToken);
				set({ userLockMessage: msg });
				return true;
			} catch (error) {
				console.error("[acquireLock] Error refreshing lock:", error);
				return false;
			}
		}
		try {
			const response = await acquireEditLock(roomId, msg);
			if (!response.locked || !response.lock_token) return false;
			const currentUserEmail = user?.email ?? null;
			set({
				userLock: currentUserEmail,
				userLockMessage: msg,
				lockToken: response.lock_token,
			});
			get().startLockRenewal();
			return true;
		} catch (error) {
			console.error("[acquireLock] Error acquiring lock:", error);
			return false;
		}
	},

	releaseLock: async () => {
		const { roomId, lockToken } = get();
		if (!lockToken) return true;
		if (!roomId) {
			console.error("[releaseLock] No roomId");
			return false;
		}
		try {
			await releaseEditLock(roomId, lockToken);
			set({ userLock: null, userLockMessage: null, lockToken: null });
			get().stopLockRenewal();
			return true;
		} catch (error) {
			console.error("[releaseLock] Error releasing lock:", error);
			return false;
		}
	},

	hasLock: () => get().lockToken !== null,

	startLockExpiryTimer: (ttl: number) => {
		get().stopLockExpiryTimer();
		const timerId = window.setTimeout(
			async () => {
				// TTL expired — verify with server
				const { roomId } = get();
				if (!roomId) return;
				try {
					const status = await getEditLockStatus(roomId);
					if (!status.locked) {
						set({ userLock: null, userLockMessage: null });
					}
				} catch {
					// Ignore — will be caught on next interaction
				}
			},
			(ttl + 1) * 1000,
		);
		set({ lockExpiryTimerId: timerId });
	},

	stopLockExpiryTimer: () => {
		const { lockExpiryTimerId } = get();
		if (lockExpiryTimerId !== null) {
			window.clearTimeout(lockExpiryTimerId);
			set({ lockExpiryTimerId: null });
		}
	},
});
