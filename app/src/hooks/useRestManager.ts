import { useEffect, useMemo, useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { useParams, useSearchParams } from "react-router-dom";
import { set, throttle } from "lodash";
import { useQueryClient } from "@tanstack/react-query";
import { joinRoom as joinRoomApi } from "../myapi/client";
import { convertBookmarkKeys } from "../utils/bookmarks";
import { ensureAuthenticated, getUsername, getUserRole } from "../utils/auth";
import { useLazyRoomData } from "./useLazyRoomData";

export const useRestJoinManager = () => {
  const {
    roomId: storeRoomId,
    userName: storeUserName,
    setRoomId,
    setUserName,
    setUserRole,
    setSessionId,
    setCurrentFrame,
    setFrameCount,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setFrameSelection,
    setBookmarks,
    setGeometries,
    setGeometryDefaults,
    setLockMetadata,
  } = useAppStore();
  const { roomId: room } = useParams<{
    roomId: string;
  }>();
  const [searchParams] = useSearchParams();
  const queryClient = useQueryClient();

  const abortControllerRef = useRef<AbortController | null>(null);

  // Lazily load room data after join succeeds
  // This hook fetches all room data in parallel when roomId is available
  const {
    isLoading: isLoadingRoomData,
    isReady: isRoomDataReady,
    hasErrors: roomDataHasErrors,
  } = useLazyRoomData(
    storeRoomId,
    !!storeRoomId, // Enable when roomId is set
  );

  // Log when room data is ready
  useEffect(() => {
    if (isRoomDataReady) {
      // Room data loaded successfully
    }
  }, [isRoomDataReady, storeRoomId]);

  const joinRoom = useCallback(async () => {
    if (!room) {
      return;
    }

    // Ensure user is authenticated before making REST API calls
    // This will auto-login if needed and ensure we have a username
    try {
      await ensureAuthenticated();

      // Initialize user role from localStorage
      const userRole = getUserRole();
      if (userRole) {
        setUserRole(userRole);
      }
    } catch (error) {
      console.error("Authentication failed:", error);
      return;
    }

    // Get username from auth system (not from URL)
    const userName = getUsername();
    if (!userName) {
      console.error("No username available after authentication");
      return;
    }

    // 2. Create a new AbortController for this specific request.
    const controller = new AbortController();
    abortControllerRef.current = controller; // Store it in the ref

    // Get template from query parameters
    const template = searchParams.get("template");

    // Build request body (userName not needed - backend uses JWT)
    const requestBody: { template?: string } = {};
    if (template) {
      requestBody.template = template;
    }

    try {
      const data = await joinRoomApi(room, requestBody, controller.signal);

      // Store session ID for this browser tab
      if (data.sessionId) {
        setSessionId(data.sessionId);
      }

      // Store essential fields only
      if (data.userName) {
        setUserName(data.userName);
      }
      setRoomId(room);

      // Room data will be fetched lazily via useLazyRoomData hook (see above)
    } catch (error: any) {
      // 4. Check if the error was due to the request being aborted.
      if (error instanceof Error && error.name === "AbortError") {
        // Fetch aborted on unmount or re-run
      } else if (error.response?.status === 401) {
        // Handle authentication error (stale token)
        const { logout, login, getUsername } = await import("../utils/auth");

        logout();

        try {
          await login();

          // Update store with new username
          const newUsername = getUsername();
          if (newUsername) {
            setUserName(newUsername);
          }

          // Retry joining the room
          joinRoom();
        } catch (loginError) {
          console.error("Re-login failed:", loginError);
        }
      } else {
        console.error("Error joining room:", error);
      }
    } finally {
      // Clean up the ref if the controller for this fetch is the current one
      if (abortControllerRef.current === controller) {
        abortControllerRef.current = null;
      }
    }
  }, [
    room,
    searchParams,
    queryClient,
    setRoomId,
    setUserName,
    setUserRole,
    setSessionId,
    setCurrentFrame,
    setFrameCount,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setFrameSelection,
    setBookmarks,
    setGeometries,
    setGeometryDefaults,
    setLockMetadata,
  ]);

  const throttledJoin = useMemo(
    () => throttle(joinRoom, 10000, { leading: true, trailing: false }),
    [joinRoom],
  );

  useEffect(() => {
    throttledJoin();

    return () => {
      abortControllerRef.current?.abort(); // Abort the fetch
      throttledJoin.cancel(); // Cancel any pending throttled execution
    };
  }, [room, throttledJoin]);

  return {
    isLoadingRoomData,
    isRoomDataReady,
    roomDataHasErrors,
  };
};
