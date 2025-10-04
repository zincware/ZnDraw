import { useState, useCallback } from "react";
import { socket } from "../socket";

export const usePresenterToken = () => {
  const [hasToken, setHasToken] = useState(false);

  const requestToken = useCallback(async (): Promise<boolean> => {
    return new Promise((resolve) => {
      socket.emit("request_presenter_token", (response: any) => {
        const success = response?.success === true;
        setHasToken(success);
        if (!success) {
          console.error("Failed to acquire presenter token:", response?.reason);
        }
        resolve(success);
      });
    });
  }, []);

  const releaseToken = useCallback(() => {
    if (hasToken) {
      socket.emit("release_presenter_token");
      setHasToken(false);
    }
  }, [hasToken]);

  const setFrame = useCallback(
    (frame: number) => {
      if (hasToken) {
        // Use continuous mode when we have the token
        socket.emit("set_frame_continuous", { frame });
      } else {
        // Use atomic mode without token
        socket.emit("set_frame_atomic", { frame }, (response: any) => {
          if (response && !response.success) {
            console.error(`Failed to set frame: ${response.error}`);
          }
        });
      }
    },
    [hasToken],
  );

  return {
    hasToken,
    requestToken,
    releaseToken,
    setFrame,
  };
};
