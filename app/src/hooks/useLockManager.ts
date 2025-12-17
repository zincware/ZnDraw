import { useCallback } from "react";
import { useAppStore } from "../store";
import { acquireLock, releaseLock } from "../myapi/client";

/**
 * Hook for managing room locks via REST API
 *
 * Provides methods to acquire, release, and manage locks with the backend.
 * Locks are used to ensure exclusive access to resources (e.g., trajectory meta, geometry editing).
 *
 * Server controls TTL and refresh intervals - clients no longer specify TTL.
 * The acquire endpoint is idempotent - calling it when you already hold the lock
 * will refresh the TTL and return the existing token.
 */
export function useLockManager() {
	const roomId = useAppStore((state) => state.roomId);

	/**
	 * Acquire a lock for a specific target (idempotent)
	 *
	 * If the same session already holds the lock, this will refresh the TTL
	 * and return the existing token with refreshed: true.
	 *
	 * @param target - Lock target (e.g., "trajectory:meta", "geometry:editing")
	 * @param msg - Optional message to display to other users
	 * @returns Promise with lock acquisition result
	 */
	const acquireLock_ = useCallback(
		async (
			target: string,
			msg?: string,
		): Promise<{
			success: boolean;
			lockToken?: string;
			ttl?: number;
			refreshInterval?: number;
			refreshed?: boolean;
		}> => {
			if (!roomId) {
				console.error("Cannot acquire lock: no room ID");
				return { success: false };
			}

			try {
				const response = await acquireLock(roomId, target, msg);
				return response;
			} catch (error: any) {
				console.error("Failed to acquire lock:", error);
				// Check if error is 423 Locked
				if (error.response?.status === 423) {
					return {
						success: false,
						...error.response.data,
					};
				}
				return { success: false };
			}
		},
		[roomId],
	);

	/**
	 * Release a lock for a specific target
	 * @param target - Lock target to release
	 * @param lockToken - Lock token from acquire response
	 * @returns Promise<boolean> - true if lock released, false otherwise
	 */
	const releaseLock_ = useCallback(
		async (target: string, lockToken: string): Promise<boolean> => {
			if (!roomId) {
				console.error("Cannot release lock: no room ID");
				return false;
			}

			try {
				const response = await releaseLock(roomId, target, lockToken);
				return response.success || false;
			} catch (error) {
				console.error("Failed to release lock:", error);
				return false;
			}
		},
		[roomId],
	);

	return {
		acquireLock: acquireLock_,
		releaseLock: releaseLock_,
	};
}
