import { useQueryClient } from "@tanstack/react-query";
import { useEffect } from "react";
import {
	createRoom,
	getAllBookmarks,
	getEditLockStatus,
	getFrameSelection,
	getGeometry,
	getGlobalSettings,
	getServerVersion,
	getSettings,
	listGeometries,
	listSelectionGroups,
} from "../myapi/client";
import { useRoomsStore } from "../roomsStore";
import { connectWithAuth, socket } from "../socket";
import { useAppStore } from "../store";
import { useWindowManagerStore } from "../stores/windowManagerStore";
import type { MessageEditedEvent, MessageNewEvent } from "../types/chat";
import { logout } from "../utils/auth";
import { setLastVisitedRoom } from "../utils/roomTracking";
import {
	checkVersionCompatibility,
	getClientVersion,
} from "../utils/versionCompatibility";

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
	const { openWindow } = useWindowManagerStore();

	// Use provided roomId from options, fallback to appStore roomId
	const roomId = options.roomId || appStoreRoomId;
	const { isOverview = false } = options;

	useEffect(() => {
		// Capture current room context for cleanup comparison
		const effectRoomId = roomId;
		const effectIsOverview = isOverview;

		/**
		 * Factory function for creating consistent invalidate handlers.
		 * Ensures uniform error handling across all handlers.
		 * Accesses roomId from the closure.
		 */
		function createInvalidateHandler<T>(
			fetchFn: (roomId: string) => Promise<T>,
			updateStoreFn: (data: T) => void,
			eventName: string,
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

		// Guard against stale async callbacks (React StrictMode double-mount).
		let cancelled = false;

		const INITIAL_RETRY_DELAY = 1000;
		const MAX_RETRY_DELAY = 30000;
		let retryDelay = INITIAL_RETRY_DELAY;

		async function onConnect() {
			retryDelay = INITIAL_RETRY_DELAY;
			// Clear any previous initialization error on new connection
			setInitializationError(null);

			try {
				// Fetch server version and global settings
				const { version: serverVersion } = await getServerVersion();
				setServerVersion(serverVersion);

				const globalSettings = await getGlobalSettings();
				setGlobalSettings(globalSettings);

				// Check version compatibility
				const clientVersion = getClientVersion();

				const compatibility = checkVersionCompatibility(
					clientVersion,
					serverVersion,
				);

				if (!compatibility.compatible) {
					console.error(compatibility.message);
					useAppStore.getState().showSnackbar(compatibility.message, "error");
					socket.disconnect();
					return;
				}

				if (compatibility.severity === "warning") {
					console.warn(compatibility.message);
				}

				// Join appropriate room based on page
				if (isOverview) {
					// Join @overview system room for real-time room updates
					socket.emit(
						"room_join",
						{ room_id: "@overview", client_type: "frontend" },
						(response: any) => {
							if (response.status) {
								console.error("Failed to join overview room:", response);
								setInitializationError({
									message: "Failed to join overview",
									details: response.detail || "Server rejected the connection",
								});
								return;
							}
							// Overview joined successfully
							setConnected(true);
						},
					);
				} else if (roomId) {
					// Handle room_join response
					const handleJoinResponse = async (response: any) => {
						// Error responses have a `status` field (RFC 9457 problem JSON).
						// Success responses are plain RoomJoinResponse (no status field).
						if (response.status) {
							console.error("Failed to join room:", response);
							setInitializationError({
								message: "Failed to join room",
								details: response.detail || "Server rejected the connection",
							});
							return;
						}

						// Reset chat unread count when entering a room
						useAppStore.getState().resetChatUnread();

						// Extract minimal data from socket response (backend uses snake_case)
						const sessionId = response.session_id;
						const step = response.step;
						const frameCount = response.frame_count;
						setSessionId(sessionId);
						setCameraKey(response.camera_key ?? null);

						// Track room visit for localStorage persistence
						setLastVisitedRoom(roomId);

						// Set minimal state from socket response
						setFrameCount(frameCount);
						setCurrentFrame(step);

						// Superuser lock from join response (SQL room.locked)
						console.debug("[RoomJoin] locked:", response.locked);
						setSuperuserLock(response.locked ?? false);

						// Set progress trackers from join response
						if (response.progress_trackers) {
							setProgressTrackers(response.progress_trackers);
						}

						// === CRITICAL: Fetch render-blocking data via REST ===
						// Both calls run in parallel - we need both to render the scene
						try {
							const [geometriesResponse] = await Promise.all([
								listGeometries(roomId),
								// Use fetchQuery to fetch AND cache settings in one step
								queryClient.fetchQuery({
									queryKey: ["settings", roomId, sessionId],
									queryFn: () => getSettings(roomId, sessionId),
								}),
							]);

							// Set geometries and type metadata (schemas + defaults)
							const geos = geometriesResponse.items || {};
							setGeometries(geos);
							setGeometrySchemas(geometriesResponse.types?.schemas || {});
							setGeometryDefaults(geometriesResponse.types?.defaults || {});

							// Extract per-geometry selections (cameras don't have selections)
							const selectionsFromGeos: Record<string, number[]> = {};
							for (const [key, geo] of Object.entries(geos)) {
								if (geo.type === "Camera") continue;
								selectionsFromGeos[key] = geo.selection ?? [];
							}
							setSelections(selectionsFromGeos);

							// Clear any previous error and mark as connected
							setInitializationError(null);
							setConnected(true);
							startHeartbeat();
						} catch (error) {
							console.error("Error fetching critical data:", error);
							setInitializationError({
								message: "Failed to load scene data",
								details:
									error instanceof Error
										? error.message
										: "Could not fetch critical data from server",
							});
							return; // Don't proceed - UI will show error state
						}

						// === SECONDARY: Fetch non-blocking data after setConnected ===
						// These can fail without breaking the UI
						try {
							const [
								groupsResponse,
								bookmarksResponse,
								frameSelResponse,
								editLockResponse,
							] = await Promise.all([
								listSelectionGroups(roomId),
								getAllBookmarks(roomId),
								getFrameSelection(roomId),
								getEditLockStatus(roomId),
							]);

							// Set selection groups
							setSelectionGroups(groupsResponse.items || {});

							// Set bookmarks
							setBookmarks(bookmarksResponse.items || {});

							// Set frame selection
							setFrameSelection(frameSelResponse.frame_selection);

							// Set edit lock status (for edit/drawing modes)
							if (editLockResponse.locked) {
								const mySessionId = useAppStore.getState().sessionId;
								if (editLockResponse.sid === mySessionId) {
									// We hold this lock (page reload case)
									useAppStore.setState({
										lockToken: editLockResponse.lock_token ?? null,
										userLock: editLockResponse.user_id ?? null,
										userLockMessage: editLockResponse.msg ?? null,
									});
									useAppStore.getState().startLockRenewal();
								} else {
									setUserLock(
										editLockResponse.user_id ?? null,
										editLockResponse.msg ?? null,
									);
									if (editLockResponse.ttl) {
										useAppStore
											.getState()
											.startLockExpiryTimer(editLockResponse.ttl);
									}
								}
							}
						} catch (error) {
							console.error("Error fetching secondary data:", error);
							// Non-critical - don't fail initialization
						}
					};

					// Emit room_join event (matches RoomJoin model → room_join)
					socket.emit(
						"room_join",
						{ room_id: roomId, client_type: "frontend" },
						async (response: any) => {
							// Handle 404 - room doesn't exist, create it via REST API
							// Backend returns RFC 9457 problem details with status field
							if (response.status === 404) {
								// Check for explicit copy_from in URL, otherwise let server decide
								// Server will use: copy_from > server default > @empty
								const urlCopyFrom = new URLSearchParams(
									window.location.search,
								).get("copy_from");
								try {
									await createRoom({
										room_id: roomId,
										copy_from: urlCopyFrom ?? undefined,
									});
								} catch (error: any) {
									// 409 Conflict = room created by another client, continue
									if (error.response?.status !== 409) {
										console.error("Failed to create room:", error);
										setInitializationError({
											message: "Failed to create room",
											details:
												error instanceof Error
													? error.message
													: "Could not create room on server",
										});
										return;
									}
								}
								// Retry join after room creation
								socket.emit(
									"room_join",
									{ room_id: roomId, client_type: "frontend" },
									handleJoinResponse,
								);
								return;
							}
							handleJoinResponse(response);
						},
					);
				} else {
					setConnected(true);
				}
			} catch (error) {
				console.error("Error checking version compatibility:", error);
				// Still connect even if version check fails
				setConnected(true);
			}
		}
		function onDisconnect() {
			setConnected(false);
			// NOTE: Do NOT clear sessionId here - onDisconnect fires during temporary
			// reconnects (e.g., when userName changes). SessionId is only cleared
			// in the cleanup function when actually leaving a room.
		}

		function onFrameUpdate(data: any) {
			const { frame } = data;
			// During local playback, ignore server echoes of our own step updates.
			// Without session_id in the event we cannot distinguish self from remote,
			// so all frame_update events are ignored while playing.
			if (useAppStore.getState().playing) return;
			setCurrentFrame(frame);
		}

		function onActiveCameraUpdate(data: { active_camera: string }) {
			const { setAttachedCameraKey } = useAppStore.getState();
			setAttachedCameraKey(data.active_camera);
		}

		function onInvalidate(data: any) {
			const { roomId, userName, category, extension, sessionId } = data;
			// Invalidate extension data queries (for modifiers, analysis, selections)
			queryClient.invalidateQueries({
				queryKey: ["extensionData", roomId, userName, category, extension],
			});
			// Invalidate settings queries (per-session)
			if (category === "settings" && sessionId) {
				queryClient.invalidateQueries({
					queryKey: ["settings", roomId, sessionId],
				});
			}
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
				setBookmarks(response.items || {});
			},
			"bookmarks",
		);

		function onFramesInvalidate(data: {
			room_id: string;
			action: "add" | "delete" | "modify" | "clear";
			indices?: number[];
			count?: number | null;
			reason?: string | null;
		}) {
			const { room_id: eventRoomId, action, indices, count, reason } = data;

			// Update frameCount if provided (new total frame count)
			if (count != null) {
				setFrameCount(count);
			}

			// Notify user when a mounted source disconnects
			if (action === "clear" && reason === "provider_disconnected") {
				useAppStore.getState().showSnackbar("Source disconnected", "warning");
			}

			// Invalidate React Query cache based on action
			queryClient.invalidateQueries({
				predicate: (query) => {
					// Query keys: ['frame', roomId, frameIndex, key]
					// or ['metadata', roomId, frameIndex]
					const [type, qRoomId, frameIndex] = query.queryKey;

					// Only invalidate frame and metadata queries for this room
					if (
						(type !== "frame" && type !== "metadata") ||
						qRoomId !== eventRoomId
					)
						return false;

					// Ensure we are only dealing with queries for specific frames
					if (typeof frameIndex !== "number") return false;

					if (action === "modify") {
						if (indices) {
							return indices.includes(frameIndex);
						}
						// No specific indices — invalidate all frame queries
						return true;
					}
					if (action === "delete" && indices) {
						// Invalidate from first deleted index onward (indices shift)
						const minDeleted = Math.min(...indices);
						return frameIndex >= minDeleted;
					}
					if (action === "add") {
						// Refresh queries with stale null data (e.g., cached 404s
						// from before frames existed)
						return query.state.data === null;
					}
					if (action === "clear") {
						return true;
					}

					return false;
				},
			});
		}

		function onChatMessageNew(data: MessageNewEvent) {
			queryClient.setQueryData(["chat", roomId], (oldData: any) => {
				if (!oldData) return oldData;
				const newPages = [...oldData.pages];
				const lastPageIndex = newPages.length - 1;

				if (lastPageIndex >= 0) {
					newPages[lastPageIndex] = {
						...newPages[lastPageIndex],
						items: [...newPages[lastPageIndex].items, data],
						metadata: {
							...newPages[lastPageIndex].metadata,
							total_count: newPages[lastPageIndex].metadata.total_count + 1,
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

		function onChatMessageUpdated(data: MessageEditedEvent) {
			queryClient.setQueryData(["chat", roomId], (oldData: any) => {
				if (!oldData) return oldData;
				const newPages = oldData.pages.map((page: any) => ({
					...page,
					items: page.items.map((msg: any) =>
						msg.id === data.id ? { ...msg, ...data } : msg,
					),
				}));
				return { ...oldData, pages: newPages };
			});
		}

		const typingTimeouts = new Map<string, ReturnType<typeof setTimeout>>();

		function onTyping(data: {
			user_id: string;
			email: string;
			is_typing: boolean;
		}) {
			const { addTypingUser, removeTypingUser } = useAppStore.getState();
			const email = data.email;

			// Clear previous timeout for this user
			const prev = typingTimeouts.get(email);
			if (prev) clearTimeout(prev);

			if (data.is_typing) {
				addTypingUser(email);
				// Auto-remove after 5s in case typing_stop is missed
				typingTimeouts.set(
					email,
					setTimeout(() => {
						removeTypingUser(email);
						typingTimeouts.delete(email);
					}, 5000),
				);
			} else {
				removeTypingUser(email);
				typingTimeouts.delete(email);
			}
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

						// Fetch only the updated geometry
						try {
							const response = await getGeometry(roomId, key);
							// Update only this specific geometry in the store
							updateGeometry(key, response.geometry);

							// Auto-select newly created curves ONLY if no curve is currently selected
							const currentActiveCurve =
								useAppStore.getState().activeCurveForDrawing;
							if (
								!currentActiveCurve &&
								response.geometry.type === "Curve" &&
								response.geometry.data?.active !== false
							) {
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
						const keys = Object.keys(listResponse.items || {});

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
						geometries.forEach(
							(item: { key: string; geometry: any } | null) => {
								if (item && item.geometry) {
									geometriesObj[item.key] = item.geometry;
								}
							},
						);

						setGeometries(geometriesObj);
					}
				}
			} catch (error) {
				console.error("Error handling geometry invalidation:", error);
			}
		}

		function onDefaultCameraInvalidate(data: {
			room_id: string;
			default_camera: string | null;
		}) {
			queryClient.setQueryData(["defaultCamera", data.room_id], {
				default_camera: data.default_camera,
			});
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
			listGeometries,
			(response) => {
				const geos = response.items || {};
				const selectionsFromGeos: Record<string, number[]> = {};
				for (const [key, geo] of Object.entries(geos)) {
					if (geo.type === "Camera") continue;
					selectionsFromGeos[key] = geo.selection ?? [];
				}
				setSelections(selectionsFromGeos);
			},
			"selections",
		);

		const onSelectionGroupsInvalidate = createInvalidateHandler(
			listSelectionGroups,
			(response) => {
				setSelectionGroups(response.items || {});
			},
			"selection_groups",
		);

		function onRoomUpdate(data: any) {
			console.debug("[RoomUpdate] received:", {
				data,
				currentRoomId: roomId,
			});

			// Update in-room state if this event is for the current room
			if (data.id === roomId) {
				if (data.frame_count != null) {
					setFrameCount(data.frame_count);
				}
				if (data.locked != null) {
					console.debug("[RoomUpdate] lock change:", data.locked);
					setSuperuserLock(data.locked);
				}
			}

			// Upsert into rooms store (full snapshot — always safe to overwrite)
			useRoomsStore.getState().setRoom(data.id, data);
		}

		function onRoomDelete(data: any) {
			const { room_id: deletedRoomId } = data;
			useRoomsStore.getState().removeRoom(deletedRoomId);
		}

		function onLockUpdate(data: any) {
			const { action, user_id, sid, msg, ttl } = data;
			const mySessionId = useAppStore.getState().sessionId;

			if (action === "acquired" || action === "refreshed") {
				// If this is our own session, lockSlice already set the state
				if (sid === mySessionId) return;
				setUserLock(user_id ?? null, msg ?? null);
				// Start TTL countdown to verify expiry
				if (ttl && action === "acquired") {
					useAppStore.getState().startLockExpiryTimer(ttl);
				}
			} else if (action === "released") {
				// If we held the lock and it was released (e.g. disconnect cleanup)
				const currentLockToken = useAppStore.getState().lockToken;
				if (currentLockToken && sid === mySessionId) {
					useAppStore.getState().stopLockRenewal();
					useAppStore.setState({
						lockToken: null,
						userLock: null,
						userLockMessage: null,
						mode: "view",
					});
					useAppStore
						.getState()
						.showSnackbar("Lock released (session ended)", "info");
				} else {
					setUserLock(null, null);
				}
				useAppStore.getState().stopLockExpiryTimer();
			}
		}

		function onProgressStarted(data: any) {
			addProgressTracker({
				progress_id: data.progress_id,
				description: data.description,
				n: 0,
				total: null,
				elapsed: 0,
				unit: data.unit ?? "it",
			});
		}

		function onProgressUpdate(data: any) {
			updateProgressTracker(data);
		}

		function onProgressComplete(data: any) {
			removeProgressTracker(data.progress_id);
		}

		async function onConnectError(err: Error) {
			if (cancelled) return;
			console.error("Socket connection error:", err.message);

			// socket.active === true means temporary network failure — auto-reconnects.
			if (socket.active) return;

			// Server rejected the connection (auth error).
			// Clear stale token and retry with backoff — guest auth should
			// always succeed unless the server is down.
			logout();
			await new Promise((r) => setTimeout(r, retryDelay));
			retryDelay = Math.min(retryDelay * 2, MAX_RETRY_DELAY);

			if (cancelled) return;
			try {
				const { user } = await connectWithAuth();
				setUser(user);
			} catch {
				// connectWithAuth failed (server down?) — next connect_error will retry
			}
		}

		// Register event handlers FIRST before connecting
		socket.on("disconnect", onDisconnect);
		socket.on("connect", onConnect);
		socket.on("connect_error", onConnectError);
		socket.on("frame_update", onFrameUpdate);
		socket.on("active_camera_update", onActiveCameraUpdate);
		socket.on("invalidate", onInvalidate);
		socket.on("schema_invalidate", onSchemaInvalidate);
		socket.on("frame_selection_update", onFrameSelectionUpdate);
		socket.on("bookmarks_invalidate", onBookmarksInvalidate);
		socket.on("frames_invalidate", onFramesInvalidate);
		socket.on("message_new", onChatMessageNew);
		socket.on("message_edited", onChatMessageUpdated);
		socket.on("typing", onTyping);
		socket.on("geometry_invalidate", onGeometriesInvalidate);
		socket.on("default_camera_invalidate", onDefaultCameraInvalidate);
		socket.on("figure_invalidate", onFiguresInvalidate);
		socket.on("selection_invalidate", onSelectionsInvalidate);
		socket.on("selection_groups_invalidate", onSelectionGroupsInvalidate);
		socket.on("room_update", onRoomUpdate);
		socket.on("room_delete", onRoomDelete);
		socket.on("lock_update", onLockUpdate);
		socket.on("progress_start", onProgressStarted);
		socket.on("progress_update", onProgressUpdate);
		socket.on("progress_complete", onProgressComplete);

		// Heartbeat interval — keeps presence + session camera alive
		const HEARTBEAT_INTERVAL = 30_000; // 30s (TTL is 60s)
		let heartbeatTimer: ReturnType<typeof setInterval> | null = null;

		function startHeartbeat() {
			stopHeartbeat();
			if (!roomId || isOverview) return;
			heartbeatTimer = setInterval(() => {
				if (socket.connected && roomId) {
					socket.emit("heartbeat", { room_id: roomId });
				}
			}, HEARTBEAT_INTERVAL);
		}

		function stopHeartbeat() {
			if (heartbeatTimer) {
				clearInterval(heartbeatTimer);
				heartbeatTimer = null;
			}
		}

		// Single auth → connect sequence
		if (socket.connected) {
			onConnect(); // room switch — re-join, don't reconnect
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
			stopHeartbeat();
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

			socket.off("connect", onConnect);
			socket.off("disconnect", onDisconnect);
			socket.off("frame_update", onFrameUpdate);
			socket.off("active_camera_update", onActiveCameraUpdate);
			socket.off("invalidate", onInvalidate);
			socket.off("schema_invalidate", onSchemaInvalidate);
			socket.off("frame_selection_update", onFrameSelectionUpdate);
			socket.off("bookmarks_invalidate", onBookmarksInvalidate);
			socket.off("frames_invalidate", onFramesInvalidate);
			socket.off("message_new", onChatMessageNew);
			socket.off("message_edited", onChatMessageUpdated);
			socket.off("typing", onTyping);
			// Clean up typing timeouts
			for (const timeout of typingTimeouts.values()) clearTimeout(timeout);
			typingTimeouts.clear();
			socket.off("geometry_invalidate", onGeometriesInvalidate);
			socket.off("default_camera_invalidate", onDefaultCameraInvalidate);
			socket.off("figure_invalidate", onFiguresInvalidate);
			socket.off("selection_invalidate", onSelectionsInvalidate);
			socket.off("selection_groups_invalidate", onSelectionGroupsInvalidate);
			socket.off("room_update", onRoomUpdate);
			socket.off("room_delete", onRoomDelete);
			socket.off("lock_update", onLockUpdate);
			socket.off("progress_start", onProgressStarted);
			socket.off("progress_update", onProgressUpdate);
			socket.off("progress_complete", onProgressComplete);

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
		openWindow,
		setCameraKey,
	]);
};
