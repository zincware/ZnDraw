import { useEffect, useMemo, useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { useParams, useSearchParams } from "react-router-dom";
import { set, throttle } from "lodash";
import { useQueryClient } from "@tanstack/react-query";
import { joinRoom as joinRoomApi } from "../myapi/client";
import { convertBookmarkKeys } from "../utils/bookmarks";
import { ensureAuthenticated, getUsername, getUserRole } from "../utils/auth";

export const useRestJoinManager = () => {
  const {
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

  const joinRoom = useCallback(async () => {
    if (!room) {
      return;
    }

    // Ensure user is authenticated before making REST API calls
    // This will auto-login if needed and ensure we have a username
    try {
      await ensureAuthenticated();
      console.log("Authentication verified for REST API");

      // Initialize user role from localStorage
      const userRole = getUserRole();
      if (userRole) {
        setUserRole(userRole);
        console.log("User role initialized:", userRole);
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

    console.log("Joining room via REST:", room, "as user:", userName);

    // Get template from query parameters
    const template = searchParams.get("template");

    // Build request body (userName not needed - backend uses JWT)
    const requestBody: { template?: string } = {};
    if (template) {
      requestBody.template = template;
    }

    try {
      const data = await joinRoomApi(room, requestBody, controller.signal);
      console.log("Join response data:", data);

      // Store session ID for this browser tab
      if (data.sessionId) {
        console.log("Setting sessionId from join:", data.sessionId);
        setSessionId(data.sessionId);
      }

      // Update Zustand store with room data
      // TODO: all of these should be fetched lazily instead of at join time
      if (typeof data.frameCount === "number") {
        setFrameCount(data.frameCount);
      }
      if (data.selections !== undefined) {
        console.log("Setting selections from join:", data.selections);
        setSelections(data.selections);
      }
      if (data.selectionGroups !== undefined) {
        console.log("Setting selection groups from join:", data.selectionGroups);
        setSelectionGroups(data.selectionGroups);
      }
      if (data.activeSelectionGroup !== undefined) {
        console.log("Setting active selection group from join:", data.activeSelectionGroup);
        setActiveSelectionGroup(data.activeSelectionGroup);
      }
      if (data.frame_selection !== undefined) {
        setFrameSelection(data.frame_selection);
      }
      if (data.step !== undefined && data.step !== null) {
        setCurrentFrame(data.step);
      }
      if (data.bookmarks !== undefined) {
        setBookmarks(convertBookmarkKeys(data.bookmarks));
      }
      if (data.userName) {
        setUserName(data.userName);
      }
      // IMPORTANT: Set defaults BEFORE geometries to avoid race condition
      // Geometries need defaults to render properly
      if (data.geometryDefaults) {
        setGeometryDefaults(data.geometryDefaults);
      }
      if (data.geometries) {
        setGeometries(data.geometries);
      }
      // Initialize lockMetadata if present in join response
      if (data.metadataLocked !== undefined) {
        if (data.metadataLocked) {
          setLockMetadata({
            locked: true,
            holder: undefined,
            userName: data.metadataLocked.userName,
            msg: data.metadataLocked.msg,
            timestamp: data.metadataLocked.timestamp,
          });
        } else {
          setLockMetadata({ locked: false });
        }
      }
      setRoomId(room);

      // Store settings in query cache for each category
      if (data.settings) {
        for (const [categoryName, categoryData] of Object.entries(
          data.settings,
        )) {
          const queryKey = [
            "extensionData",
            room,
            userName,
            "settings",
            categoryName,
          ];
          queryClient.setQueryData(queryKey, categoryData);
        }
      }
    } catch (error: any) {
      // 4. Check if the error was due to the request being aborted.
      if (error instanceof Error && error.name === "AbortError") {
        console.log("Fetch aborted on unmount or re-run.");
      } else if (error.response?.status === 401) {
        // Handle authentication error (stale token)
        console.log("Authentication error, clearing stale token and retrying...");
        const { logout, login, getUsername } = await import("../utils/auth");

        logout();

        try {
          await login();
          console.log("Re-login successful, retrying join room...");

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

  return {};
};
