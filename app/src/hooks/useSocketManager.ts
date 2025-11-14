import { useEffect, useRef } from "react";
import { socket } from "../socket";
import { useAppStore } from "../store";
import { useQueryClient } from "@tanstack/react-query";
import { useWindowManagerStore } from "../stores/windowManagerStore";
import { listGeometries, getGeometry, getAllSelections, getAllBookmarks, getServerVersion, getGlobalSettings } from "../myapi/client";
import { convertBookmarkKeys } from "../utils/bookmarks";
import { checkVersionCompatibility, getClientVersion } from "../utils/versionCompatibility";
import { useRoomsStore } from "../roomsStore";
import { ensureAuthenticated } from "../utils/auth";

const MAX_AUTH_RETRIES = 3;
const INITIAL_RETRY_DELAY = 1000; // 1 second


interface SocketManagerOptions {
  roomId?: string;  // Room ID when on /rooms/:roomId page
  isOverview?: boolean;  // True when on /rooms page
}

export const useSocketManager = (options: SocketManagerOptions = {}) => {
  const {
    setConnected,
    setFrameCount,
    setCurrentFrame,
    setFrameSelection,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setBookmarks,
    roomId: appStoreRoomId,
    userName,
    setUserName,
    setGeometries,
    updateGeometry,
    removeGeometry,
    setActiveCurveForDrawing,
    setServerVersion,
    setGlobalSettings,
    setLockMetadata,
    setProgressTrackers,
    addProgressTracker,
    updateProgressTracker,
    removeProgressTracker,
  } = useAppStore();
  const queryClient = useQueryClient();
  const { openWindow } = useWindowManagerStore();

  // Use provided roomId from options, fallback to appStore roomId
  const roomId = options.roomId || appStoreRoomId;
  const { isOverview = false } = options;

  // Track retry attempts for stale token recovery
  const authRetryCountRef = useRef(0);

  useEffect(() => {
    /**
     * Factory function for creating consistent invalidate handlers.
     * Ensures uniform error handling across all handlers.
     * Accesses roomId from the closure.
     */
    function createInvalidateHandler<T>(
      fetchFn: (roomId: string) => Promise<T>,
      updateStoreFn: (data: T) => void,
      eventName: string
    ) {
      return async (data: any) => {
        if (!roomId) return;
        try {
          const response = await fetchFn(roomId);
          updateStoreFn(response);
        } catch (error) {
          console.error(`Error fetching ${eventName}:`, error);
        }
      };
    }

    async function onConnect() {
      try {
        // Fetch server version and global settings
        const { version: serverVersion } = await getServerVersion();
        setServerVersion(serverVersion);

        const globalSettings = await getGlobalSettings();
        setGlobalSettings(globalSettings);

        // Check version compatibility
        const clientVersion = getClientVersion();

        const compatibility = checkVersionCompatibility(clientVersion, serverVersion);

        if (!compatibility.compatible) {
          console.error(compatibility.message);
          alert(compatibility.message);
          socket.disconnect();
          return;
        }

        if (compatibility.severity === 'warning') {
          console.warn(compatibility.message);
        }

        // Join appropriate room based on page
        if (isOverview) {
          socket.emit("join:overview");
        } else if (roomId) {
          socket.emit("join:room", { roomId });
          // Reset chat unread count when entering a room
          useAppStore.getState().resetChatUnread();
        }

        setConnected(true);
      } catch (error) {
        console.error("Error checking version compatibility:", error);
        // Still connect even if version check fails
        setConnected(true);
      }
    }
    function onDisconnect() {
      setConnected(false);
    }

    function onFrameUpdate(data: any) {
      const { frame } = data;
      setCurrentFrame(frame);
    }

    function onInvalidate(data: any) {
      const { roomId, userName, category, extension } = data;
      queryClient.invalidateQueries({
        queryKey: ["extensionData", roomId, userName, category, extension],
      });
    }

    function onSchemaInvalidate(data: any) {
      const { category } = data;
      queryClient.invalidateQueries({
        queryKey: ["schemas", appStoreRoomId, category],
      });
    }

    function onFrameSelectionUpdate(data: any) {
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
      const { roomId, operation, affectedIndex, affectedFrom, affectedKeys } = data;

      if (affectedKeys) {
        console.warn("Frame invalidation with affectedKeys is not yet implemented. Invalidating all frames.");
      }

      queryClient.invalidateQueries({
        predicate: (query) => {
          // Query keys are expected to be in the format: ['frame', roomId, frameIndex]
          // or ['frame-keys', roomId, frameIndex]
          const [type, qRoomId, frameIndex] = query.queryKey;

          // Only invalidate frame and frame-keys queries for this room
          if ((type !== "frame" && type !== "frame-keys") || qRoomId !== roomId) return false;

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

      // Increment unread count if chat is closed
      const { chatOpen, incrementChatUnread } = useAppStore.getState();
      if (!chatOpen) {
        incrementChatUnread();
      }
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
        } else if (operation === "set") {
          // Handle geometry creation/update
          if (data && data.key) {
            const { key } = data;

            // Invalidate the specific geometry detail query
            queryClient.invalidateQueries({
              queryKey: ["geometries", roomId, "detail", key],
            });

            // Invalidate the list to update the geometry grid
            queryClient.invalidateQueries({
              queryKey: ["geometries", roomId, "list"],
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

              // Auto-select newly created curves ONLY if no curve is currently selected
              const currentActiveCurve = useAppStore.getState().activeCurveForDrawing;
              if (!currentActiveCurve && response.geometry.type === "Curve" && response.geometry.data?.active !== false) {
                setActiveCurveForDrawing(key);
              }
            } catch (error) {
              console.error(`Error fetching geometry ${key}:`, error);
            }
          } else {
            // No specific key - refetch all geometries (fallback for backward compatibility)

            // Invalidate React Query cache for geometries
            queryClient.invalidateQueries({
              queryKey: ["geometries", roomId, "list"],
            });

            // Fetch list of geometry keys first
            const listResponse = await listGeometries(roomId);
            const keys = listResponse.geometries || [];

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
            openWindow(data.key);
          }
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

    function onRoomUpdate(data: any) {
      const { roomId: updateRoomId, created, ...updates } = data;

      // Handle frame count update
      if ("frameCount" in updates) {
        setFrameCount(updates.frameCount);
      }

      // Update Zustand rooms store
      if (created) {
        // New room created
        useRoomsStore.getState().setRoom(updateRoomId, {
          id: updateRoomId,
          frameCount: 0,
          locked: false,
          hidden: false,
          isDefault: false,
          ...updates,
        });
      } else {
        // Update existing room
        useRoomsStore.getState().updateRoom(updateRoomId, updates);
      }
    }

    function onRoomDelete(data: any) {
      const { roomId: deletedRoomId } = data;
      useRoomsStore.getState().removeRoom(deletedRoomId);
    }

    function onLockUpdate(data: any) {
      const { roomId: lockRoomId, target, action, holder, message, timestamp, sessionId } = data;

      // Handle trajectory:meta lock updates
      // in the UI (lock icon should show current user as holder)
      if (target === "trajectory:meta") {
        if (action === "acquired" || action === "refreshed") {
          setLockMetadata({
            locked: true,
            holder: holder,
            userName: holder,
            msg: message,
            timestamp: timestamp,
          });
        } else if (action === "released") {
          setLockMetadata({ locked: false });
        }
      }

      // Step lock updates are handled by useStepControl hook
      // (no action needed here)

      // Future: handle other lock targets
      // if (target === "geometry:selection") { ... }
    }

    function onProgressInitial(data: any) {
      const { progressTrackers } = data;
      if (progressTrackers) {
        setProgressTrackers(progressTrackers);
      }
    }

    function onProgressStarted(data: any) {
      const { progressId, roomId, description } = data;
      addProgressTracker(progressId, description, null, roomId);
    }

    function onProgressUpdate(data: any) {
      const { progressId, description, progress } = data;
      updateProgressTracker(progressId, description, progress);
    }

    function onProgressComplete(data: any) {
      const { progressId } = data;
      removeProgressTracker(progressId);
    }

    async function onConnectError(err: any) {
      console.error("Socket connection error:", err);

      // Handle "Client not registered" error (stale token in localStorage)
      // This happens when Redis is flushed but localStorage still has old token
      if (err.message && err.message.includes("Client not registered")) {
        // Check if we've exceeded max retries
        if (authRetryCountRef.current >= MAX_AUTH_RETRIES) {
          console.error(`Max authentication retries (${MAX_AUTH_RETRIES}) exceeded. Please refresh the page.`);
          // TODO: Show user-friendly error message
          return;
        }

        authRetryCountRef.current += 1;
        const retryAttempt = authRetryCountRef.current;

        // Calculate exponential backoff delay: 1s, 2s, 4s
        const delay = INITIAL_RETRY_DELAY * Math.pow(2, retryAttempt - 1);

        const { logout, login, getUsername } = await import("../utils/auth");

        // Clear stale authentication data
        logout();

        // Wait for exponential backoff delay
        await new Promise(resolve => setTimeout(resolve, delay));

        // Force a new login
        try {
          const loginData = await login();

          // Update store with new username
          const newUsername = getUsername();
          if (newUsername) {
            setUserName(newUsername);
          }

          // Reset retry count on successful login
          authRetryCountRef.current = 0;

          // Reconnect socket with new credentials
          socket.connect();
        } catch (loginError) {
          console.error(`Re-login failed (attempt ${retryAttempt}/${MAX_AUTH_RETRIES}):`, loginError);
          // Will retry on next connect_error event if under max retries
        }
      }
    }

    // Register event handlers FIRST before connecting
    socket.on("disconnect", onDisconnect);
    socket.on("connect", onConnect);
    socket.on("connect_error", onConnectError);
    socket.on("frame_update", onFrameUpdate);
    socket.on("invalidate", onInvalidate);
    socket.on("invalidate:schema", onSchemaInvalidate);
    socket.on("frame_selection:update", onFrameSelectionUpdate);
    socket.on("bookmarks:invalidate", onBookmarksInvalidate);
    socket.on("frames:invalidate", onFramesInvalidate);
    socket.on("chat:message:new", onChatMessageNew);
    socket.on("chat:message:updated", onChatMessageUpdated);
    socket.on("invalidate:geometry", onGeometriesInvalidate);
    socket.on("invalidate:figure", onFiguresInvalidate);
    socket.on("invalidate:selection", onSelectionsInvalidate);
    socket.on("invalidate:selection_groups", onSelectionGroupsInvalidate);
    socket.on("room:update", onRoomUpdate);
    socket.on("room:delete", onRoomDelete);
    socket.on("lock:update", onLockUpdate);
    socket.on("progress:initial", onProgressInitial);
    socket.on("progress:started", onProgressStarted);
    socket.on("progress:updated", onProgressUpdate);
    socket.on("progress:completed", onProgressComplete);

    // Ensure user is authenticated before connecting socket
    // This will auto-login with a server-generated username if no token exists
    ensureAuthenticated()
      .then(() => {
        // Force disconnect/reconnect when userName changes to ensure clean state
        if (socket.connected) {
          // Disconnect synchronously
          socket.disconnect();
          // Small delay to ensure disconnect completes
          setTimeout(() => {
            socket.connect();
          }, 100);
        } else {
          socket.connect();
        }
      })
      .catch((error) => {
        console.error("Authentication failed:", error);
      });

    return () => {
      // Leave room on cleanup
      if (isOverview) {
        socket.emit("leave:overview");
      } else if (roomId) {
        socket.emit("leave:room", { roomId });
      }

      socket.off("connect", onConnect);
      socket.off("disconnect", onDisconnect);
      socket.off("frame_update", onFrameUpdate);
      socket.off("invalidate", onInvalidate);
      socket.off("invalidate:schema", onSchemaInvalidate);
      socket.off("frame_selection:update", onFrameSelectionUpdate);
      socket.off("bookmarks:invalidate", onBookmarksInvalidate);
      socket.off("frames:invalidate", onFramesInvalidate);
      socket.off("chat:message:new", onChatMessageNew);
      socket.off("chat:message:updated", onChatMessageUpdated);
      socket.off("invalidate:geometry", onGeometriesInvalidate);
      socket.off("invalidate:figure", onFiguresInvalidate);
      socket.off("invalidate:selection", onSelectionsInvalidate);
      socket.off("invalidate:selection_groups", onSelectionGroupsInvalidate);
      socket.off("room:update", onRoomUpdate);
      socket.off("room:delete", onRoomDelete);
      socket.off("lock:update", onLockUpdate);
      socket.off("progress:initial", onProgressInitial);
      socket.off("progress:started", onProgressStarted);
      socket.off("progress:updated", onProgressUpdate);
      socket.off("progress:completed", onProgressComplete);

      // Disconnect socket when component unmounts to ensure clean reconnection
      socket.disconnect();
    };
  }, [
    roomId,
    userName,
    isOverview,
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
    setServerVersion,
    setGlobalSettings,
    setLockMetadata,
    setProgressTrackers,
    addProgressTracker,
    updateProgressTracker,
    removeProgressTracker,
    openWindow,
  ]);
};
