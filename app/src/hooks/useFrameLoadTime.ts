import { useEffect, useRef } from "react";
import { useAppStore } from "../store";

/**
 * Hook to track frame load time when not playing.
 * Measures the time from when a frame changes until the scene is fully rendered.
 * Uses Zustand subscription for efficient change detection instead of RAF polling.
 */
export const useFrameLoadTime = () => {
  const {
    currentFrame,
    playing,
    setFrameLoadTime,
  } = useAppStore();

  const frameStartTimeRef = useRef<number | null>(null);
  const previousFrameRef = useRef<number | null>(null);
  const rafIdRef = useRef<number | null>(null);
  const unsubscribeRef = useRef<(() => void) | null>(null);

  useEffect(() => {
    // Only track load time when not playing
    if (playing) {
      setFrameLoadTime(null);
      frameStartTimeRef.current = null;
      previousFrameRef.current = null;
      // Cancel any pending RAF
      if (rafIdRef.current !== null) {
        cancelAnimationFrame(rafIdRef.current);
        rafIdRef.current = null;
      }
      // Unsubscribe from store changes
      if (unsubscribeRef.current !== null) {
        unsubscribeRef.current();
        unsubscribeRef.current = null;
      }
      return;
    }

    // If frame changed, start timing
    if (previousFrameRef.current !== currentFrame) {
      // Cancel any previous pending measurement
      if (rafIdRef.current !== null) {
        cancelAnimationFrame(rafIdRef.current);
        rafIdRef.current = null;
      }
      // Unsubscribe from previous frame's subscription
      if (unsubscribeRef.current !== null) {
        unsubscribeRef.current();
        unsubscribeRef.current = null;
      }

      const startTime = performance.now();
      frameStartTimeRef.current = startTime;
      previousFrameRef.current = currentFrame;

      // Helper to check if any active geometry is still fetching
      const checkIfFetching = () => {
        const state = useAppStore.getState();
        const { geometryFetchingStates, geometries } = state;

        // Early bailout: no fetching states at all
        if (Object.keys(geometryFetchingStates).length === 0) {
          return false;
        }

        // Check if any active geometry is fetching
        return Object.entries(geometryFetchingStates).some(([key, isFetching]) => {
          if (!isFetching) return false; // Early bailout for non-fetching
          const geometry = geometries[key];
          const isActive = geometry?.data?.active !== false;
          return isActive;
        });
      };

      // Helper to complete measurement (called when all fetching done)
      const completeMeasurement = () => {
        // Wait one more frame to ensure render is complete
        rafIdRef.current = requestAnimationFrame(() => {
          const loadTime = Math.round(performance.now() - startTime);
          setFrameLoadTime(loadTime);
          rafIdRef.current = null;
        });

        // Clean up subscription after measurement starts
        if (unsubscribeRef.current !== null) {
          unsubscribeRef.current();
          unsubscribeRef.current = null;
        }
      };

      // Check immediately if already done
      const initiallyFetching = checkIfFetching();
      if (!initiallyFetching) {
        completeMeasurement();
        return;
      }

      // Track previous fetching states to detect changes
      let previousFetchingStates = useAppStore.getState().geometryFetchingStates;

      // Subscribe to store changes (event-driven, not polling)
      // Only check when geometryFetchingStates actually changes
      unsubscribeRef.current = useAppStore.subscribe((state) => {
        const currentFetchingStates = state.geometryFetchingStates;

        // Only check if geometryFetchingStates reference changed
        if (currentFetchingStates !== previousFetchingStates) {
          previousFetchingStates = currentFetchingStates;

          const stillFetching = checkIfFetching();
          if (!stillFetching) {
            completeMeasurement();
          }
        }
      });
    }

    // Cleanup on unmount or dependencies change
    return () => {
      if (rafIdRef.current !== null) {
        cancelAnimationFrame(rafIdRef.current);
        rafIdRef.current = null;
      }
      if (unsubscribeRef.current !== null) {
        unsubscribeRef.current();
        unsubscribeRef.current = null;
      }
    };
  }, [currentFrame, playing, setFrameLoadTime]);
};
