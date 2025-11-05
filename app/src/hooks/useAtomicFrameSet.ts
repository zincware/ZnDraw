import { useCallback } from "react";
import { socket } from "../socket";
import { useAppStore } from "../store";

/**
 * Hook for setting the current frame with atomic synchronization.
 *
 * This ensures frame changes are:
 * - Updated locally in the store
 * - Broadcast to all other clients via the server
 * - Rejected if a presenter is currently active
 *
 * Use this for discrete frame jumps (clicks, keyboard nav, manual input).
 * Do NOT use during continuous updates (slider dragging, playback) - use
 * usePresenterToken for those cases.
 */
export const useAtomicFrameSet = () => {
  // Use individual selectors to prevent unnecessary re-renders
  const setCurrentFrame = useAppStore((state) => state.setCurrentFrame);

  return useCallback(
    (frame: number) => {
      // Update local state immediately for responsive UI
      setCurrentFrame(frame);

      // Sync with server and other clients
      socket.emit("set_frame_atomic", { frame }, (response: any) => {
        if (response && !response.success) {
          console.error(
            `Failed to set frame: ${response.error}${response.message ? ` - ${response.message}` : ""}`,
          );
        }
      });
    },
    [setCurrentFrame],
  );
};
