import { useEffect, useMemo, useCallback, useRef } from "react";
import { useAppStore } from "../store";
import { useParams, useSearchParams } from "react-router-dom";
import { set, throttle } from "lodash";
import { useQueryClient } from "@tanstack/react-query";
import { joinRoom as joinRoomApi } from "../myapi/client";

export const useRestJoinManager = () => {
  const {
    setClientId,
    setRoomId,
    setUserId,
    setCurrentFrame,
    setFrameCount,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setFrameSelection,
    setBookmarks,
    setJoinToken,
    setGeometries,
  } = useAppStore();
  const { roomId: room, userId } = useParams<{
    roomId: string;
    userId: string;
  }>();
  const [searchParams] = useSearchParams();
  const queryClient = useQueryClient();

  const abortControllerRef = useRef<AbortController | null>(null);

  const joinRoom = useCallback(async () => {
    if (!room || !userId) {
      return;
    }

    // 2. Create a new AbortController for this specific request.
    const controller = new AbortController();
    abortControllerRef.current = controller; // Store it in the ref

    console.log("Joining room via REST:", room, userId);

    // Get template from query parameters
    const template = searchParams.get("template");

    // Build request body
    const requestBody: { userId: string; template?: string } = { userId };
    if (template) {
      requestBody.template = template;
    }

    try {
      const data = await joinRoomApi(room, requestBody, controller.signal);
      console.log("Join response data:", data);

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
        setBookmarks(data.bookmarks);
      }
      if (data.clientId) {
        setClientId(data.clientId);
      }
      if (data.joinToken) {
        setJoinToken(data.joinToken);
      }
      if (data.geometries) {
        setGeometries(data.geometries);
      }
      setRoomId(room);
      setUserId(userId);

      // Store settings in query cache for each category
      if (data.settings) {
        for (const [categoryName, categoryData] of Object.entries(
          data.settings,
        )) {
          const queryKey = [
            "extensionData",
            room,
            userId,
            "settings",
            categoryName,
          ];
          queryClient.setQueryData(queryKey, categoryData);
        }
      }
    } catch (error) {
      // 4. Check if the error was due to the request being aborted.
      if (error instanceof Error && error.name === "AbortError") {
        console.log("Fetch aborted on unmount or re-run.");
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
    userId,
    searchParams,
    queryClient,
    setClientId,
    setRoomId,
    setUserId,
    setCurrentFrame,
    setFrameCount,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setFrameSelection,
    setBookmarks,
    setJoinToken,
    setGeometries,
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
  }, [room, userId, throttledJoin]);

  return {};
};
