import { useEffect } from "react";
import { socket } from "../socket";
import { useAppStore } from "../store";
import { useQueryClient } from "@tanstack/react-query";
import { useWindowManagerStore } from "../stores/windowManagerStore";
import { listGeometries, getGeometry, getAllSelections, getAllBookmarks } from "../myapi/client";
import { convertBookmarkKeys } from "../utils/bookmarks";

/**
 * Factory function for creating consistent invalidate handlers.
 * Ensures uniform error handling across all handlers.
 */
function createInvalidateHandler<T>(
  fetchFn: (roomId: string) => Promise<T>,
  updateStoreFn: (data: T) => void,
  eventName: string
) {
  return async (data: any, roomId: string | null) => {
    if (!roomId) return;
    try {
      console.log(`Received ${eventName} event:`, data);
      const response = await fetchFn(roomId);
      updateStoreFn(response);
      console.log(`Updated ${eventName} from server:`, response);
    } catch (error) {
      console.error(`Error fetching ${eventName}:`, error);
    }
  };
}

export const useSocketManager = () => {
  const {
    setConnected,
    setFrameCount,
    setCurrentFrame,
    setFrameSelection,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setBookmarks,
    roomId,
    userId,
    joinToken,
    setGeometries,
    updateGeometry,
    removeGeometry,
    setActiveCurveForDrawing,
  } = useAppStore();
  const queryClient = useQueryClient();
  const { openWindow } = useWindowManagerStore();

  useEffect(() => {
    if (!joinToken) {
      console.warn("No join token available, cannot connect socket.");
      return;
    }

    async function onConnect() {
      console.log("Socket connected and joining room:", roomId, userId);
      setConnected(true);
    }
    function onDisconnect() {
      console.log("Socket disconnected");
      setConnected(false);
    }
    function onLenUpdate(data: any) {
      if (data && typeof data.count === "number") {
        setFrameCount(data.count);
      } else {
        console.error("Invalid len_frames data:", data);
      }
    }

    function onFrameUpdate(data: any) {
      const { frame } = data;
      setCurrentFrame(frame);
    }

    function onInvalidate(data: any) {
      const { roomId, userId, category, extension } = data;
      queryClient.invalidateQueries({
        queryKey: ["extensionData", roomId, userId, category, extension],
      });
      console.log(
        `Invalidated extension data for user ${userId}, category ${category}, extension ${extension} in room ${roomId}`,
      );
    }

    function onSchemaInvalidate(data: any) {
      const { roomId, category } = data;
      queryClient.invalidateQueries({
        queryKey: ["schemas", roomId, category],
      });
      console.log(
        `Invalidated schemas for category ${category} in room ${roomId}`,
      );
    }

    function onQueueUpdate(data: any) {
      const {
        roomId,
        category,
        extension,
        queueLength,
        idleWorkers,
        progressingWorkers,
      } = data;

      // 1. Update the cached schema data with new queue length and worker counts
      queryClient.setQueryData(
        ["schemas", roomId, category],
        (oldData: any) => {
          if (!oldData || !oldData[extension]) return oldData;
          return {
            ...oldData,
            [extension]: {
              ...oldData[extension],
              queueLength,
              idleWorkers,
              progressingWorkers,
            },
          };
        },
      );

      // 2. Invalidate jobs list to refetch all jobs (covers created/started/completed/failed/deleted)
      queryClient.invalidateQueries({ queryKey: ["jobs", roomId] });

      console.log(
        `Queue updated for ${category}/${extension}: ${queueLength} queued, ${idleWorkers} idle, ${progressingWorkers} progressing`,
      );
    }

    function onFrameSelectionUpdate(data: any) {
      console.log(data);
      setFrameSelection(data["indices"] || null);
    }

    // Create handlers using factory for consistency
    const onBookmarksInvalidate = createInvalidateHandler(
      getAllBookmarks,
      (response) => {
        const bookmarksWithNumberKeys = convertBookmarkKeys(response.bookmarks);
        setBookmarks(bookmarksWithNumberKeys);
      },
      "bookmarks"
    );

    function onFramesInvalidate(data: any) {
      const { roomId, operation, affectedIndex, affectedFrom } = data;
      console.log("Invalidating frame cache:", {
        roomId,
        operation,
        affectedIndex,
        affectedFrom,
      });

      queryClient.invalidateQueries({
        predicate: (query) => {
          // Query keys are expected to be in the format: ['frame', roomId, frameIndex]
          const [type, qRoomId, frameIndex] = query.queryKey;

          // Only invalidate frame queries for this room
          if (type !== "frame" || qRoomId !== roomId) return false;

          // Ensure we are only dealing with queries for specific frames
          if (typeof frameIndex !== "number") return false;

          if (operation === "replace") {
            // Handles single-item replacement ONLY.
            return frameIndex === affectedIndex;
          } else if (
            operation === "delete" ||
            operation === "insert" ||
            operation === "bulk_replace"
          ) {
            // Handles any operation that shifts indices or modifies a range.
            // Invalidate all frames from the affected position onward.
            return frameIndex >= affectedFrom;
          }

          // Default to not invalidating if the operation is unknown.
          return false;
        },
      });
    }

    function onChatMessageNew(data: any) {
      queryClient.setQueryData(["chat", roomId], (oldData: any) => {
        if (!oldData) return oldData;
        const newPages = [...oldData.pages];
        const lastPageIndex = newPages.length - 1;
        if (lastPageIndex >= 0) {
          newPages[lastPageIndex] = {
            ...newPages[lastPageIndex],
            messages: [...newPages[lastPageIndex].messages, data],
            metadata: {
              ...newPages[lastPageIndex].metadata,
              totalCount: newPages[lastPageIndex].metadata.totalCount + 1,
            },
          };
        }
        return { ...oldData, pages: newPages };
      });
    }

    function onChatMessageUpdated(data: any) {
      queryClient.setQueryData(["chat", roomId], (oldData: any) => {
        if (!oldData) return oldData;
        const newPages = oldData.pages.map((page: any) => ({
          ...page,
          messages: page.messages.map((msg: any) =>
            msg.id === data.id ? data : msg,
          ),
        }));
        return { ...oldData, pages: newPages };
      });
    }

    async function onGeometriesInvalidate(data: any) {
      if (!roomId) return;

      try {
        const operation = data?.operation || "set"; // default to 'set' for backward compatibility

        if (operation === "delete") {
          // Handle geometry deletion
          const { key } = data;
          if (!key) {
            console.warn("Delete operation received without key");
            return;
          }

          console.log(`Deleting geometry: ${key}`);

          // Remove from the store
          removeGeometry(key);

          // Remove the specific geometry from the cache
          queryClient.removeQueries({
            queryKey: ["geometries", roomId, "detail", key],
          });

          // Invalidate the list to update the UI
          queryClient.invalidateQueries({
            queryKey: ["geometries", roomId, "list"],
          });

          console.log(`Removed geometry ${key} from store and cache`);
        } else if (operation === "set") {
          // Handle geometry creation/update
          if (data && data.key) {
            const { key } = data;
            console.log(`Updating specific geometry: ${key}`);

            // Invalidate the specific geometry detail query
            queryClient.invalidateQueries({
              queryKey: ["geometries", roomId, "detail", key],
            });

            // Also invalidate frame queries - geometry data may reference dynamic frame data
            // that needs to be refetched (e.g., "arrays.position")
            queryClient.invalidateQueries({
              queryKey: ["frame", roomId],
            });

            // Fetch only the updated geometry
            try {
              const response = await getGeometry(roomId, key);
              // Update only this specific geometry in the store
              updateGeometry(key, response.geometry);
              console.log(`Updated geometry ${key}:`, response.geometry);

              // Auto-select newly created curves ONLY if no curve is currently selected
              const currentActiveCurve = useAppStore.getState().activeCurveForDrawing;
              if (!currentActiveCurve && response.geometry.type === "Curve" && response.geometry.data?.active !== false) {
                console.log(`Auto-selecting newly created curve (no curve was selected): ${key}`);
                setActiveCurveForDrawing(key);
              }
            } catch (error) {
              console.error(`Error fetching geometry ${key}:`, error);
            }
          } else {
            // No specific key - refetch all geometries (fallback for backward compatibility)
            console.log("No specific key provided, refetching all geometries");

            // Invalidate React Query cache for geometries
            queryClient.invalidateQueries({
              queryKey: ["geometries", roomId, "list"],
            });

            // Fetch list of geometry keys first
            const listResponse = await listGeometries(roomId);
            const keys = listResponse.geometries || [];
            console.log("Received geometry keys:", keys);

            // Fetch all geometries in parallel
            const geometryPromises = keys.map(async (key: string) => {
              try {
                const response = await getGeometry(roomId, key);
                return { key, geometry: response.geometry };
              } catch (error) {
                console.error(`Error fetching geometry ${key}:`, error);
                return null;
              }
            });

            const geometries = await Promise.all(geometryPromises);

            // Build geometries object from results
            const geometriesObj: Record<string, any> = {};
            geometries.forEach((item) => {
              if (item && item.geometry) {
                geometriesObj[item.key] = item.geometry;
              }
            });

            console.log("Received updated geometries:", geometriesObj);
            setGeometries(geometriesObj);
          }
        }
      } catch (error) {
        console.error("Error handling geometry invalidation:", error);
      }
    }

    function onFiguresInvalidate(data: {
      key: string;
      operation?: "set" | "delete";
    }) {
      console.log(
        `Received invalidation for figure: ${data.key}, operation: ${data.operation}`,
      );
      if (!data.key) return;

      const operation = data.operation || "set"; // default to 'set' for backward compatibility

      if (operation === "delete") {
        // Step 1: Remove the figure data from the cache
        queryClient.removeQueries({
          queryKey: ["figures", roomId, "detail", data.key],
        });

        // Step 2: Close any open windows displaying this figure
        // Find all windows showing this figure key
        const windowsToClose = Object.entries(
          useWindowManagerStore.getState().openWindows,
        )
          .filter(([_, window]) => window.figureKey === data.key)
          .map(([windowId]) => windowId);

        windowsToClose.forEach((windowId) => {
          console.log(
            `Closing window ${windowId} for deleted figure: ${data.key}`,
          );
          useWindowManagerStore.getState().closeWindow(windowId);
        });

        // Step 3: Invalidate the figures list to update the UI
        queryClient.invalidateQueries({
          queryKey: ["figures", roomId, "list"],
        });
      } else if (operation === "set") {
        // Step 1: Invalidate the data so it's fresh when needed
        queryClient.invalidateQueries({
          queryKey: ["figures", roomId, "detail", data.key],
        });

        // Step 2: Invalidate the list in case this is a new figure
        queryClient.invalidateQueries({
          queryKey: ["figures", roomId, "list"],
        });

        // Step 3: Check if ANY window is already displaying this figure key
        // We need to check if any window.figureKey matches, not just if openWindows[data.key] exists
        const openWindowsState = useWindowManagerStore.getState().openWindows;
        const isWindowDisplayingFigure = Object.values(openWindowsState).some(
          (window) => window.figureKey === data.key,
        );

        const openWindowCount = Object.keys(openWindowsState).length;
        const MAX_AUTO_OPEN_WINDOWS = 5;

        if (!isWindowDisplayingFigure) {
          if (openWindowCount < MAX_AUTO_OPEN_WINDOWS) {
            console.log(
              `Popping up new window for figure: ${data.key} (${openWindowCount + 1}/${MAX_AUTO_OPEN_WINDOWS})`,
            );
            openWindow(data.key);
          } else {
            console.log(
              `Not opening window for figure: ${data.key} - maximum of ${MAX_AUTO_OPEN_WINDOWS} windows already open`,
            );
          }
        } else {
          console.log(
            `Window already open for figure: ${data.key}, not opening a new one`,
          );
        }
      }
    }

    const onSelectionsInvalidate = createInvalidateHandler(
      getAllSelections,
      (response) => {
        setSelections(response.selections);
        setSelectionGroups(response.groups);
        setActiveSelectionGroup(response.activeGroup);
      },
      "selections"
    );

    const onSelectionGroupsInvalidate = createInvalidateHandler(
      getAllSelections,
      (response) => {
        setSelectionGroups(response.groups);
        setActiveSelectionGroup(response.activeGroup);
      },
      "selection_groups"
    );

    function onConnectError(err: any) {
      console.error("Socket connection error:", err);
    }

    socket.on("disconnect", onDisconnect);
    socket.on("connect", onConnect);
    socket.on("connect_error", onConnectError);
    socket.on("len_frames", onLenUpdate);
    socket.on("frame_update", onFrameUpdate);
    socket.on("invalidate", onInvalidate);
    socket.on("invalidate:schema", onSchemaInvalidate);
    socket.on("queue:update", onQueueUpdate);
    socket.on("frame_selection:update", onFrameSelectionUpdate);
    socket.on("bookmarks:invalidate", onBookmarksInvalidate);
    socket.on("frames:invalidate", onFramesInvalidate);
    socket.on("chat:message:new", onChatMessageNew);
    socket.on("chat:message:updated", onChatMessageUpdated);
    socket.on("invalidate:geometry", onGeometriesInvalidate);
    socket.on("invalidate:figure", onFiguresInvalidate);
    socket.on("invalidate:selection", onSelectionsInvalidate);
    socket.on("invalidate:selection_groups", onSelectionGroupsInvalidate);

    socket.auth = { token: joinToken };
    socket.connect();

    return () => {
      socket.off("connect", onConnect);
      socket.off("disconnect", onDisconnect);
      socket.off("len_frames", onLenUpdate);
      socket.off("frame_update", onFrameUpdate);
      socket.off("invalidate", onInvalidate);
      socket.off("invalidate:schema", onSchemaInvalidate);
      socket.off("queue:update", onQueueUpdate);
      socket.off("frame_selection:update", onFrameSelectionUpdate);
      socket.off("bookmarks:invalidate", onBookmarksInvalidate);
      socket.off("frames:invalidate", onFramesInvalidate);
      socket.off("chat:message:new", onChatMessageNew);
      socket.off("chat:message:updated", onChatMessageUpdated);
      socket.off("invalidate:geometry", onGeometriesInvalidate);
      socket.off("invalidate:figure", onFiguresInvalidate);
      socket.off("invalidate:selection", onSelectionsInvalidate);
      socket.off("invalidate:selection_groups", onSelectionGroupsInvalidate);

      // Disconnect socket when component unmounts to ensure clean reconnection
      socket.disconnect();
    };
  }, [
    joinToken,
    roomId,
    userId,
    setConnected,
    setFrameCount,
    setCurrentFrame,
    queryClient,
    setBookmarks,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setFrameSelection,
    setGeometries,
    updateGeometry,
    removeGeometry,
    setActiveCurveForDrawing,
    openWindow,
  ]);
};
