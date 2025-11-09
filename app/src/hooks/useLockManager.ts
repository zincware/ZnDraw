import { useCallback } from "react";
import { useAppStore } from "../store";
import { acquireLock, refreshLock, releaseLock } from "../myapi/client";

/**
 * Hook for managing room locks via REST API
 *
 * Provides methods to acquire, release, and manage locks with the backend.
 * Locks are used to ensure exclusive access to resources (e.g., trajectory meta, geometry editing).
 *
 * Server controls TTL and refresh intervals - clients no longer specify TTL.
 */
export function useLockManager() {
  const roomId = useAppStore((state) => state.roomId);

  /**
   * Acquire a lock for a specific target
   * @param target - Lock target (e.g., "trajectory:meta", "geometry:editing")
   * @param msg - Optional message to display to other users
   * @returns Promise<{success: boolean, lockToken?: string, ttl?: number, refreshInterval?: number}> - Lock acquisition result with server-provided lock token, TTL and refresh interval
   */
  const acquireLock_ = useCallback(
    async (
      target: string,
      msg?: string
    ): Promise<{ success: boolean; lockToken?: string; ttl?: number; refreshInterval?: number }> => {
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
    [roomId]
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
    [roomId]
  );

  /**
   * Refresh a lock to extend its TTL and optionally update message
   * @param target - Lock target to refresh
   * @param lockToken - Lock token from acquire response
   * @param msg - Optional updated message (if provided, updates the lock message)
   * @returns Promise<boolean> - true if lock refreshed, false otherwise
   */
  const refreshLock_ = useCallback(
    async (target: string, lockToken: string, msg?: string): Promise<boolean> => {
      if (!roomId) {
        console.error("Cannot refresh lock: no room ID");
        return false;
      }

      try {
        const response = await refreshLock(roomId, target, lockToken, msg);
        return response.success || false;
      } catch (error) {
        console.error("Failed to refresh lock:", error);
        return false;
      }
    },
    [roomId]
  );

  /**
   * Update the message associated with a lock
   * This is just an alias for refreshLock with a message parameter
   * @param target - Lock target
   * @param lockToken - Lock token from acquire response
   * @param msg - New message to display
   * @returns Promise<boolean> - true if message updated, false otherwise
   */
  const updateLockMessage = useCallback(
    async (target: string, lockToken: string, msg: string): Promise<boolean> => {
      return refreshLock_(target, lockToken, msg);
    },
    [refreshLock_]
  );

  return {
    acquireLock: acquireLock_,
    releaseLock: releaseLock_,
    refreshLock: refreshLock_,
    updateLockMessage,
  };
}
