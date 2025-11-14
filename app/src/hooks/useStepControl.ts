import { useRef, useCallback, useEffect, useState, useMemo } from "react";
import { throttle } from "lodash";
import { socket } from "../socket";
import { useLockManager } from "./useLockManager";
import { useAppStore } from "../store";
import { updateStep } from "../myapi/client";

interface StepControlHook {
  setStep: (frame: number) => Promise<void>;
  remoteLocked: boolean;
}

/**
 * Unified hook for step/frame control with automatic lock management.
 *
 * Features:
 * - Single setStep() function for all use cases (clicks, keyboard, slider, playback)
 * - Auto-acquires lock on first step update
 * - Debounced lock release (500ms after last update)
 * - Throttled server updates (100ms) with trailing edge to ensure final value synced
 * - Automatic heartbeat while lock is held
 * - Filters lock:update events using sessionId (ignores own events)
 *
 * Usage:
 *   const { setStep, remoteLocked } = useStepControl();
 *   await setStep(42); // Updates local state, syncs to server, manages lock automatically
 *   <Slider disabled={remoteLocked} />
 *
 * remoteLocked logic:
 * - false if no one has the lock
 * - false if WE have the lock
 * - true ONLY if someone else has the lock
 */
export const useStepControl = (): StepControlHook => {
  const [remoteLocked, setRemoteLocked] = useState<boolean>(false);

  const lockTokenRef = useRef<string | null>(null);
  const releaseDebounceRef = useRef<NodeJS.Timeout | null>(null);
  const heartbeatIntervalRef = useRef<NodeJS.Timeout | null>(null);

  const { acquireLock, releaseLock, refreshLock } = useLockManager();
  const roomId = useAppStore((state) => state.roomId);
  const setCurrentFrame = useAppStore((state) => state.setCurrentFrame);

  // Start heartbeat to refresh lock every 30 seconds
  const startHeartbeat = useCallback(() => {
    if (heartbeatIntervalRef.current) return; // Already running

    heartbeatIntervalRef.current = setInterval(async () => {
      if (lockTokenRef.current && roomId) {
        try {
          const success = await refreshLock("step", lockTokenRef.current, "Updating");
          if (!success) {
            console.error("Step lock refresh failed. Releasing lock state.");
            lockTokenRef.current = null;
            if (heartbeatIntervalRef.current) {
              clearInterval(heartbeatIntervalRef.current);
              heartbeatIntervalRef.current = null;
            }
          }
        } catch (error) {
          console.error("Failed to refresh step lock:", error);
          lockTokenRef.current = null;
          if (heartbeatIntervalRef.current) {
            clearInterval(heartbeatIntervalRef.current);
            heartbeatIntervalRef.current = null;
          }
        }
      }
    }, 30000); // 30 seconds
  }, [roomId, refreshLock]);

  // Stop heartbeat
  const stopHeartbeat = useCallback(() => {
    if (heartbeatIntervalRef.current) {
      clearInterval(heartbeatIntervalRef.current);
      heartbeatIntervalRef.current = null;
    }
  }, []);

  // Throttled server update (100ms throttle with trailing edge)
  const throttledUpdate = useMemo(
    () =>
      throttle(
        async (frame: number) => {
          if (!roomId || !lockTokenRef.current) return;

          try {
            await updateStep(roomId, frame);
          } catch (error: any) {
            console.error("Failed to update step:", error);
            // If we get 423 (lock lost), release our local lock state
            if (error.response?.status === 423) {
              lockTokenRef.current = null;
              stopHeartbeat();
            }
          }
        },
        100, // Max 10 updates/second
        { leading: true, trailing: true }, // Execute first call immediately + last call after throttle
      ),
    [roomId, stopHeartbeat],
  );

  // Release lock and cleanup all state
  const releaseLockAndCleanup = useCallback(async () => {
    if (lockTokenRef.current && roomId) {
      // Flush any pending throttled update before releasing lock
      await throttledUpdate.flush();

      try {
        await releaseLock("step", lockTokenRef.current);
      } catch (error) {
        console.error("Failed to release step lock:", error);
      } finally {
        lockTokenRef.current = null;
        stopHeartbeat();
      }
    }
  }, [roomId, releaseLock, throttledUpdate, stopHeartbeat]);

  // Universal step setter
  const setStep = useCallback(
    async (frame: number) => {
      if (!roomId) {
        console.error("Cannot set step: no room ID");
        return;
      }

      // 1. Update local state immediately for responsive UI
      setCurrentFrame(frame);

      // 2. Acquire lock if we don't have it
      if (!lockTokenRef.current) {
        try {
          const response = await acquireLock("step", "Updating");

          if (!response.success || !response.lockToken) {
            console.warn("[setStep] Failed to acquire step lock - someone else has it");
            return;
          }

          lockTokenRef.current = response.lockToken;
          startHeartbeat();
        } catch (error) {
          console.error("[setStep] Failed to acquire step lock:", error);
          return;
        }
      }

      // 3. Update server (throttled)
      throttledUpdate(frame);

      // 4. Reset debounce timer - release lock 500ms after last update
      if (releaseDebounceRef.current) {
        clearTimeout(releaseDebounceRef.current);
      }

      releaseDebounceRef.current = setTimeout(() => {
        releaseLockAndCleanup();
      }, 500);
    },
    [roomId, setCurrentFrame, acquireLock, startHeartbeat, throttledUpdate, releaseLockAndCleanup],
  );

  // Handle lock:update events from other clients
  useEffect(() => {
    const onLockUpdate = (data: any) => {
      if (data.target !== "step") return;

      const mySessionId = useAppStore.getState().sessionId;

      // If this is OUR lock event, ignore it (we already have local state)
      if (data.sessionId && data.sessionId === mySessionId) {
        return;
      }

      // This is from another client - update remoteLocked state
      if (data.action === "released") {
        setRemoteLocked(false);
      } else if (data.action === "acquired" || data.action === "refreshed") {
        setRemoteLocked(true);
      }
    };

    socket.on("lock:update", onLockUpdate);

    return () => {
      socket.off("lock:update", onLockUpdate);
    };
  }, []);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      // Cancel any pending debounced release
      if (releaseDebounceRef.current) {
        clearTimeout(releaseDebounceRef.current);
      }

      // Flush throttled updates
      throttledUpdate.cancel();

      // Release lock if we're holding it
      if (lockTokenRef.current && roomId) {
        releaseLock("step", lockTokenRef.current).catch(console.error);
      }

      stopHeartbeat();
    };
  }, [roomId, releaseLock, throttledUpdate, stopHeartbeat]);

  return {
    setStep,
    remoteLocked,
  };
};
