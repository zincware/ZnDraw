import { useEffect } from "react";
import { useQuery, useQueryClient } from "@tanstack/react-query";
import {
  getRoomInfo,
  getAllSelections,
  getFrameSelection,
  getCurrentStep,
  getAllBookmarks,
  getUserRoomSettings,
  listGeometries,
} from "../myapi/client";
import { useGeometrySchemas } from "./useGeometries";
import { useAppStore } from "../store";

// Query cache stale time constants (in milliseconds)
const ROOM_INFO_STALE_TIME = 30000;      // 30 seconds - room metadata changes infrequently
const SELECTIONS_STALE_TIME = 10000;     // 10 seconds - selections change moderately
const FRAME_SELECTION_STALE_TIME = 10000; // 10 seconds - frame selection changes moderately
const STEP_STALE_TIME = 5000;            // 5 seconds - current step/frame changes frequently
const BOOKMARKS_STALE_TIME = 30000;      // 30 seconds - bookmarks change infrequently
const SETTINGS_STALE_TIME = 60000;       // 60 seconds (1 minute) - user settings rarely change
const GEOMETRIES_STALE_TIME = 30000;     // 30 seconds - geometries change infrequently

/**
 * Hook to lazily load all room data after joining a room.
 *
 * This hook fetches room data in parallel and populates the Zustand store.
 * Data includes: frameCount, selections, geometries, bookmarks, settings, etc.
 *
 * User identity is handled via JWT authentication, so userName is not needed
 * in cache keys.
 *
 * @param roomId - The room to fetch data for
 * @param enabled - Whether to start fetching (typically after join succeeds)
 */
export const useLazyRoomData = (
  roomId: string | null,
  enabled: boolean,
) => {
  const queryClient = useQueryClient();
  const {
    setFrameCount,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setFrameSelection,
    setCurrentFrame,
    setBookmarks,
    setGeometries,
    setGeometryDefaults,
  } = useAppStore();

  // Fetch room info (includes frameCount)
  const {
    data: roomInfo,
    isLoading: isLoadingRoomInfo,
    error: roomInfoError,
  } = useQuery({
    queryKey: ["roomInfo", roomId],
    queryFn: () => getRoomInfo(roomId!),
    enabled: enabled && !!roomId,
    staleTime: ROOM_INFO_STALE_TIME,
    retry: 2,
  });

  // Fetch selections (includes selections, groups, activeGroup)
  const {
    data: selectionsData,
    isLoading: isLoadingSelections,
    error: selectionsError,
  } = useQuery({
    queryKey: ["selections", roomId],
    queryFn: () => getAllSelections(roomId!),
    enabled: enabled && !!roomId,
    staleTime: SELECTIONS_STALE_TIME,
    retry: 2,
  });

  // Fetch frame selection
  const {
    data: frameSelectionData,
    isLoading: isLoadingFrameSelection,
    error: frameSelectionError,
  } = useQuery({
    queryKey: ["frameSelection", roomId],
    queryFn: () => getFrameSelection(roomId!),
    enabled: enabled && !!roomId,
    staleTime: FRAME_SELECTION_STALE_TIME,
    retry: 2,
  });

  // Fetch current step/frame
  const {
    data: stepData,
    isLoading: isLoadingStep,
    error: stepError,
  } = useQuery({
    queryKey: ["currentStep", roomId],
    queryFn: () => getCurrentStep(roomId!),
    enabled: enabled && !!roomId,
    staleTime: STEP_STALE_TIME,
    retry: 2,
  });

  // Fetch bookmarks
  const {
    data: bookmarksData,
    isLoading: isLoadingBookmarks,
    error: bookmarksError,
  } = useQuery({
    queryKey: ["bookmarks", roomId],
    queryFn: () => getAllBookmarks(roomId!),
    enabled: enabled && !!roomId,
    staleTime: BOOKMARKS_STALE_TIME,
    retry: 2,
  });

  // Fetch user settings (JWT auth handles user identity)
  const {
    data: settingsData,
    isLoading: isLoadingSettings,
    error: settingsError,
  } = useQuery({
    queryKey: ["userSettings", roomId],
    queryFn: () => getUserRoomSettings(roomId!),
    enabled: enabled && !!roomId,
    staleTime: SETTINGS_STALE_TIME,
    retry: 2,
  });

  // Fetch geometries with full data
  const {
    data: geometriesData,
    isLoading: isLoadingGeometries,
    error: geometriesError,
  } = useQuery({
    queryKey: ["geometries", roomId],
    queryFn: () => listGeometries(roomId!),
    enabled: enabled && !!roomId,
    staleTime: GEOMETRIES_STALE_TIME,
    retry: 2,
  });

  // Fetch geometry schemas (defaults)
  const { data: geometrySchemasData, isLoading: isLoadingSchemas } = useGeometrySchemas(
    enabled && roomId ? roomId : null,
  );

  // Update Zustand store when data arrives
  useEffect(() => {
    if (roomInfo?.frameCount !== undefined) {
      setFrameCount(roomInfo.frameCount);
    }
  }, [roomInfo, setFrameCount]);

  useEffect(() => {
    if (selectionsData?.selections !== undefined) {
      setSelections(selectionsData.selections);
    }
    if (selectionsData?.groups !== undefined) {
      setSelectionGroups(selectionsData.groups);
    }
    if (selectionsData?.activeGroup !== undefined) {
      setActiveSelectionGroup(selectionsData.activeGroup);
    }
  }, [selectionsData, setSelections, setSelectionGroups, setActiveSelectionGroup]);

  useEffect(() => {
    if (frameSelectionData?.frameSelection !== undefined) {
      setFrameSelection(frameSelectionData.frameSelection);
    }
  }, [frameSelectionData, setFrameSelection]);

  useEffect(() => {
    if (stepData?.step !== undefined && stepData.step !== null) {
      setCurrentFrame(stepData.step);
    }
  }, [stepData, setCurrentFrame]);

  useEffect(() => {
    if (bookmarksData?.bookmarks !== undefined) {
      // Convert bookmark keys to numbers if needed
      const convertedBookmarks =
        bookmarksData.bookmarks === null
          ? null
          : Object.entries(bookmarksData.bookmarks).reduce(
              (acc, [key, value]) => {
                acc[Number(key)] = value;
                return acc;
              },
              {} as Record<number, string>,
            );
      setBookmarks(convertedBookmarks);
    }
  }, [bookmarksData, setBookmarks]);

  useEffect(() => {
    if (geometriesData?.geometries !== undefined) {
      setGeometries(geometriesData.geometries);
    }
  }, [geometriesData, setGeometries]);

  useEffect(() => {
    if (geometrySchemasData?.schemas !== undefined) {
      setGeometryDefaults(geometrySchemasData.schemas);
    }
  }, [geometrySchemasData, setGeometryDefaults]);

  // Settings are handled differently - they're stored in query cache per category
  // The useExtensionData hook will read from the cache when needed
  useEffect(() => {
    if (settingsData?.settings && roomId) {
      // Store settings in query cache for each category
      for (const [categoryName, categoryData] of Object.entries(
        settingsData.settings,
      )) {
        const queryKey = [
          "extensionData",
          roomId,
          "settings",
          categoryName,
        ];
        queryClient.setQueryData(queryKey, categoryData);
      }
    }
  }, [settingsData, roomId, queryClient]);

  // Aggregate loading state
  const isLoading =
    isLoadingRoomInfo ||
    isLoadingSelections ||
    isLoadingFrameSelection ||
    isLoadingStep ||
    isLoadingBookmarks ||
    isLoadingSettings ||
    isLoadingGeometries ||
    isLoadingSchemas;

  // Aggregate error state
  const hasErrors = !!(
    roomInfoError ||
    selectionsError ||
    frameSelectionError ||
    stepError ||
    bookmarksError ||
    settingsError ||
    geometriesError
  );

  // Log errors for debugging (only when enabled to avoid spam)
  useEffect(() => {
    if (!enabled) return;

    if (roomInfoError) {
      console.error("[useLazyRoomData] Failed to load room info:", roomInfoError);
    }
    if (selectionsError) {
      console.error("[useLazyRoomData] Failed to load selections:", selectionsError);
    }
    if (frameSelectionError) {
      console.error("[useLazyRoomData] Failed to load frame selection:", frameSelectionError);
    }
    if (stepError) {
      console.error("[useLazyRoomData] Failed to load current step:", stepError);
    }
    if (bookmarksError) {
      console.error("[useLazyRoomData] Failed to load bookmarks:", bookmarksError);
    }
    if (settingsError) {
      console.error("[useLazyRoomData] Failed to load settings:", settingsError);
    }
    if (geometriesError) {
      console.error("[useLazyRoomData] Failed to load geometries:", geometriesError);
    }
  }, [
    enabled,
    roomInfoError,
    selectionsError,
    frameSelectionError,
    stepError,
    bookmarksError,
    settingsError,
    geometriesError,
  ]);

  return {
    isLoading,
    isReady: !isLoading && !hasErrors && enabled,
    hasErrors,
  };
};
