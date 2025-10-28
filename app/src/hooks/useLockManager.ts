import { useCallback } from "react";
import { socket } from "../socket";

/**
 * Hook for managing room locks via socket.io
 *
 * Provides methods to acquire, release, and manage locks with the backend.
 * Locks are used to ensure exclusive access to resources (e.g., trajectory meta, geometry editing).
 */
export function useLockManager() {
  /**
   * Acquire a lock for a specific target
   * @param target - Lock target (e.g., "trajectory:meta", "geometry:editing")
   * @param ttl - Time to live in seconds (default 60, max 300)
   * @param msg - Optional message to display to other users
   * @returns Promise<boolean> - true if lock acquired, false if already held by another client
   */
  const acquireLock = useCallback(
    async (target: string, ttl: number = 60, msg?: string): Promise<boolean> => {
      return new Promise((resolve) => {
        socket.emit("lock:acquire", { target, ttl }, (response: any) => {
          if (response.success && msg) {
            // Update lock message after acquiring
            socket.emit("lock:msg", { target, msg });
          }
          resolve(response.success || false);
        });
      });
    },
    []
  );

  /**
   * Release a lock for a specific target
   * @param target - Lock target to release
   * @returns Promise<boolean> - true if lock released, false otherwise
   */
  const releaseLock = useCallback(
    async (target: string): Promise<boolean> => {
      return new Promise((resolve) => {
        socket.emit("lock:release", { target }, (response: any) => {
          resolve(response.success || false);
        });
      });
    },
    []
  );

  /**
   * Refresh a lock to extend its TTL
   * @param target - Lock target to refresh
   * @param ttl - New time to live in seconds (default 60, max 300)
   * @returns Promise<boolean> - true if lock refreshed, false otherwise
   */
  const refreshLock = useCallback(
    async (target: string, ttl: number = 60): Promise<boolean> => {
      return new Promise((resolve) => {
        socket.emit("lock:refresh", { target, ttl }, (response: any) => {
          resolve(response.success || false);
        });
      });
    },
    []
  );

  /**
   * Update the message associated with a lock
   * @param target - Lock target
   * @param msg - New message to display
   * @returns Promise<boolean> - true if message updated, false otherwise
   */
  const updateLockMessage = useCallback(
    async (target: string, msg: string): Promise<boolean> => {
      return new Promise((resolve) => {
        socket.emit("lock:msg", { target, msg }, (response: any) => {
          resolve(response.success || false);
        });
      });
    },
    []
  );

  return {
    acquireLock,
    releaseLock,
    refreshLock,
    updateLockMessage,
  };
}
