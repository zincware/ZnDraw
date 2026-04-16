import { useQueryClient } from "@tanstack/react-query";
import { useEffect } from "react";
import { connectWithAuth, socket } from "../socket";
import { useAppStore } from "../store";
import {
	createChatHandlers,
	createConnectionHandlers,
	createFigureHandlers,
	createFrameHandlers,
	createGeometryHandlers,
	createRoomHandlers,
	createSceneInvalidationHandlers,
	type HandlerContext,
} from "./socketHandlers";

interface SocketManagerOptions {
	roomId?: string; // Room ID when on /rooms/:roomId page
	isOverview?: boolean; // True when on /rooms page
}

export const useSocketManager = (options: SocketManagerOptions = {}) => {
	// Selector-based subscriptions to avoid re-renders on unrelated state changes
	const appStoreRoomId = useAppStore((state) => state.roomId);
	const setConnected = useAppStore((state) => state.setConnected);
	const setInitializationError = useAppStore(
		(state) => state.setInitializationError,
	);
	const setFrameCount = useAppStore((state) => state.setFrameCount);
	const setCurrentFrame = useAppStore((state) => state.setCurrentFrame);
	const setFrameSelection = useAppStore((state) => state.setFrameSelection);
	const setSelections = useAppStore((state) => state.setSelections);
	const setSelectionGroups = useAppStore((state) => state.setSelectionGroups);
	const setBookmarks = useAppStore((state) => state.setBookmarks);
	const setUser = useAppStore((state) => state.setUser);
	const setSessionId = useAppStore((state) => state.setSessionId);
	const setCameraKey = useAppStore((state) => state.setCameraKey);
	const setGeometries = useAppStore((state) => state.setGeometries);
	const setGeometrySchemas = useAppStore((state) => state.setGeometrySchemas);
	const setGeometryDefaults = useAppStore((state) => state.setGeometryDefaults);
	const updateGeometry = useAppStore((state) => state.updateGeometry);
	const removeGeometry = useAppStore((state) => state.removeGeometry);
	const setActiveCurveForDrawing = useAppStore(
		(state) => state.setActiveCurveForDrawing,
	);
	const setServerVersion = useAppStore((state) => state.setServerVersion);
	const setGlobalSettings = useAppStore((state) => state.setGlobalSettings);
	const setSuperuserLock = useAppStore((state) => state.setSuperuserLock);
	const setUserLock = useAppStore((state) => state.setUserLock);
	const setProgressTrackers = useAppStore((state) => state.setProgressTrackers);
	const addProgressTracker = useAppStore((state) => state.addProgressTracker);
	const updateProgressTracker = useAppStore(
		(state) => state.updateProgressTracker,
	);
	const removeProgressTracker = useAppStore(
		(state) => state.removeProgressTracker,
	);
	const queryClient = useQueryClient();

	// Use provided roomId from options, fallback to appStore roomId
	const roomId = options.roomId || appStoreRoomId;
	const { isOverview = false } = options;

	useEffect(() => {
		// Capture current room context for cleanup comparison
		const effectRoomId = roomId;
		const effectIsOverview = isOverview;

		// Guard against stale async callbacks (React StrictMode double-mount).
		let cancelled = false;

		// Build flat dependency context for all handler factories
		const ctx: HandlerContext = {
			roomId: roomId ?? undefined,
			appStoreRoomId,
			isOverview,
			isCancelled: () => cancelled,
			queryClient,
			setConnected,
			setInitializationError,
			setUser,
			setSessionId,
			setCameraKey,
			setServerVersion,
			setGlobalSettings,
			setFrameCount,
			setCurrentFrame,
			setFrameSelection,
			setBookmarks,
			setSelections,
			setSelectionGroups,
			setGeometries,
			setGeometrySchemas,
			setGeometryDefaults,
			updateGeometry,
			removeGeometry,
			setActiveCurveForDrawing,
			setSuperuserLock,
			setUserLock,
			setProgressTrackers,
			addProgressTracker,
			updateProgressTracker,
			removeProgressTracker,
		};

		// Create handler groups from factories
		const connection = createConnectionHandlers(ctx);
		const frames = createFrameHandlers(ctx);
		const geometries = createGeometryHandlers(ctx);
		const { handlers: chat, cleanup: chatCleanup } = createChatHandlers(ctx);
		const sceneInvalidation = createSceneInvalidationHandlers(ctx);
		const figures = createFigureHandlers(ctx);
		const room = createRoomHandlers(ctx);

		// Register event handlers FIRST before connecting
		socket.on("connect", connection.onConnect);
		socket.on("disconnect", connection.onDisconnect);
		socket.on("connect_error", connection.onConnectError);
		socket.on("frame_update", frames.onFrameUpdate);
		socket.on("active_camera_update", geometries.onActiveCameraUpdate);
		socket.on("invalidate", sceneInvalidation.onInvalidate);
		socket.on("schema_invalidate", sceneInvalidation.onSchemaInvalidate);
		socket.on("frame_selection_update", frames.onFrameSelectionUpdate);
		socket.on("bookmarks_invalidate", geometries.onBookmarksInvalidate);
		socket.on("frames_invalidate", frames.onFramesInvalidate);
		socket.on("message_new", chat.onChatMessageNew);
		socket.on("message_edited", chat.onChatMessageUpdated);
		socket.on("typing", chat.onTyping);
		socket.on("geometry_invalidate", geometries.onGeometriesInvalidate);
		socket.on(
			"default_camera_invalidate",
			geometries.onDefaultCameraInvalidate,
		);
		socket.on("figure_invalidate", figures.onFiguresInvalidate);
		socket.on("selection_invalidate", geometries.onSelectionsInvalidate);
		socket.on(
			"selection_groups_invalidate",
			geometries.onSelectionGroupsInvalidate,
		);
		socket.on("room_update", room.onRoomUpdate);
		socket.on("room_delete", room.onRoomDelete);
		socket.on("lock_update", room.onLockUpdate);
		socket.on("progress_start", room.onProgressStarted);
		socket.on("progress_update", room.onProgressUpdate);
		socket.on("progress_complete", room.onProgressComplete);

		// Single auth -> connect sequence
		if (socket.connected) {
			connection.onConnect(); // room switch -- re-join, don't reconnect
		} else {
			connectWithAuth()
				.then(({ user }) => {
					if (cancelled) return;
					setUser(user);
				})
				.catch((error) => {
					if (cancelled) return;
					console.error("Authentication failed:", error);
				});
		}

		return () => {
			cancelled = true;
			// Note: room_leave is NOT emitted during room switching
			// Backend automatically handles leaving old room in room_join handler
			const currentState = useAppStore.getState();
			const isNavigatingToOverview =
				!effectIsOverview && options.isOverview === true;
			const isLeavingRoom =
				effectRoomId !== currentState.roomId ||
				effectIsOverview !== options.isOverview;

			if (isNavigatingToOverview && effectRoomId) {
				// Explicitly leaving a room to go to overview
				socket.emit("room_leave", { room_id: effectRoomId });
			} else if (isOverview) {
				// Leaving overview system room
				socket.emit("room_leave", { room_id: "@overview" });
			}

			// Only clear sessionId/cameraKey if we're actually leaving a room
			if (isLeavingRoom) {
				setSessionId(null);
				setCameraKey(null);
			}

			chatCleanup();

			socket.off("connect", connection.onConnect);
			socket.off("disconnect", connection.onDisconnect);
			socket.off("connect_error", connection.onConnectError);
			socket.off("frame_update", frames.onFrameUpdate);
			socket.off("active_camera_update", geometries.onActiveCameraUpdate);
			socket.off("invalidate", sceneInvalidation.onInvalidate);
			socket.off("schema_invalidate", sceneInvalidation.onSchemaInvalidate);
			socket.off("frame_selection_update", frames.onFrameSelectionUpdate);
			socket.off("bookmarks_invalidate", geometries.onBookmarksInvalidate);
			socket.off("frames_invalidate", frames.onFramesInvalidate);
			socket.off("message_new", chat.onChatMessageNew);
			socket.off("message_edited", chat.onChatMessageUpdated);
			socket.off("typing", chat.onTyping);
			socket.off("geometry_invalidate", geometries.onGeometriesInvalidate);
			socket.off(
				"default_camera_invalidate",
				geometries.onDefaultCameraInvalidate,
			);
			socket.off("figure_invalidate", figures.onFiguresInvalidate);
			socket.off("selection_invalidate", geometries.onSelectionsInvalidate);
			socket.off(
				"selection_groups_invalidate",
				geometries.onSelectionGroupsInvalidate,
			);
			socket.off("room_update", room.onRoomUpdate);
			socket.off("room_delete", room.onRoomDelete);
			socket.off("lock_update", room.onLockUpdate);
			socket.off("progress_start", room.onProgressStarted);
			socket.off("progress_update", room.onProgressUpdate);
			socket.off("progress_complete", room.onProgressComplete);

			// REMOVED: socket.disconnect()
			// We now maintain persistent Socket.IO connections across navigation
			// Connection will disconnect on: tab close, network error, or explicit logout
		};
	}, [
		roomId,
		isOverview,
		setConnected,
		setInitializationError,
		setFrameCount,
		setCurrentFrame,
		queryClient,
		setBookmarks,
		setSelections,
		setSelectionGroups,
		setFrameSelection,
		setUser,
		setSessionId,
		setGeometries,
		setGeometrySchemas,
		setGeometryDefaults,
		updateGeometry,
		removeGeometry,
		setActiveCurveForDrawing,
		setServerVersion,
		setGlobalSettings,
		setSuperuserLock,
		setUserLock,
		setProgressTrackers,
		addProgressTracker,
		updateProgressTracker,
		removeProgressTracker,
		setCameraKey,
	]);
};
