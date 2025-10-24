import { useRef, useCallback, useEffect, useState } from "react";
import { socket } from "../socket";

interface PresenterModeHook {
  requestPresenterMode: () => void;
  releasePresenterMode: () => void;
  presenterMode: PresenterMode;
}

type PresenterMode = "idle" | "requesting" | "presenting" | "locked";

export const usePresenterMode = (): PresenterModeHook => {
  const [presenterMode, setPresenterMode] = useState<PresenterMode>("idle");
  const presenterLockTimeout = useRef<NodeJS.Timeout | null>(null);

  // 1. Create a ref to hold the heartbeat interval ID
  const heartbeatInterval = useRef<NodeJS.Timeout | null>(null);

  const requestPresenterMode = useCallback(() => {
    // Prevent multiple requests while one is in flight
    if (presenterMode === "requesting") return;

    setPresenterMode("requesting");
    socket.emit("request_presenter_token", (response: { success: boolean }) => {
      if (!response.success) {
        console.warn("Presenter mode request denied");
      } else {
        setPresenterMode("presenting");
      }
    });
  }, [socket, presenterMode]);

  const releasePresenterMode = useCallback(() => {
    if (presenterMode === "presenting") {
      socket.emit("release_presenter_token");
      setPresenterMode("idle");
    }
  }, [presenterMode, socket]);

  // 2. Add a new useEffect to manage the heartbeat interval
  useEffect(() => {
    // If we are the presenter, start the heartbeat
    if (presenterMode === "presenting") {
      // Start a new interval
      heartbeatInterval.current = setInterval(() => {
        console.log("Sending presenter heartbeat...");
        // We can reuse the same event. The server logic now handles renewals.
        socket.emit(
          "request_presenter_token",
          (response: { success: boolean }) => {
            // If our heartbeat fails for any reason (e.g., server restarted, lock lost),
            // we should stop presenting.
            if (!response.success) {
              console.error(
                "Presenter heartbeat failed. Releasing presenter mode.",
              );
              setPresenterMode("idle"); // Or 'locked' if a presenter_update says so
            }
          },
        );
      }, 3000); // Send heartbeat every 3 seconds
    }

    return () => {
      if (heartbeatInterval.current) {
        clearInterval(heartbeatInterval.current);
        heartbeatInterval.current = null;
      }
    };
  }, [presenterMode, socket]);

  // Handle presenter events from others
  useEffect(() => {
    const onRoomUpdate = (data: any) => {
      // Only handle presenterSid updates
      if (!("presenterSid" in data)) {
        return;
      }

      console.log("Presenter update received:", data);

      if (presenterLockTimeout.current) {
        clearTimeout(presenterLockTimeout.current);
      }

      if (data.presenterSid === null) {
        setPresenterMode("idle");
      } else {
        setPresenterMode("locked");

        presenterLockTimeout.current = setTimeout(() => {
          console.warn(
            "No presenter update received for 5 seconds. Reverting to idle.",
          );
          setPresenterMode("idle");
        }, 5000);
      }
    };

    socket.on("room:update", onRoomUpdate);

    return () => {
      socket.off("room:update", onRoomUpdate);

      if (presenterLockTimeout.current) {
        clearTimeout(presenterLockTimeout.current);
      }

      if (presenterMode === "presenting") {
        socket.emit("release_presenter_token");
      }
    };
  }, [socket, presenterMode]);

  return {
    requestPresenterMode,
    releasePresenterMode,
    presenterMode,
  };
};
