import { useEffect, useRef } from "react";
import { socket } from "../socket";
import { useAppStore } from "../store";
import { useQueryClient } from "@tanstack/react-query";
import { useWindowManagerStore } from "../stores/windowManagerStore";
import {
	listGeometries,
	getGeometry,
	getAllSelections,
	getAllBookmarks,
	getFrameSelection,
	getSettings,
	getServerVersion,
	getGlobalSettings,
	createRoom,
	getLockStatus,
} from "../myapi/client";
import {
	checkVersionCompatibility,
	getClientVersion,
} from "../utils/versionCompatibility";
import { useRoomsStore } from "../roomsStore";
import { ensureAuthenticated } from "../utils/auth";
import { setLastVisitedRoom } from "../utils/roomTracking";

const MAX_AUTH_RETRIES = 3;
const INITIAL_RETRY_DELAY = 1000; // 1 second

interface SocketManagerOptions {
	roomId?: string; // Room ID when on /rooms/:roomId page
	isOverview?: boolean; // True when on /rooms page
}

export const useSocketManager = (options: SocketManagerOptions = {}) => {
	const {
		setConnected,
		setInitializationError,
		setFrameCount,
		setCurrentFrame,
		setFrameSelection,
		setSelections,
		setSelectionGroups,
		setActiveSelectionGroup,
		setBookmarks,
		roomId: appStoreRoomId,
		setUserName,
		setSessionId,
		setGeometries,
		setGeometrySchemas,
		setGeometryDefaults,
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
		setPlaying,
	} = useAppStore();
	const queryClient = useQueryClient();
	const { openWindow } = useWindowManagerStore();

	// Use provided roomId from options, fallback to appStore roomId
	const roomId = options.roomId || appStoreRoomId;
	const { isOverview = false } = options;

	// Track retry attempts for stale token recovery
	const authRetryCountRef = useRef(0);
	// Guard against concurrent auth recovery (multiple connect_error events)
	const isRecoveringAuthRef = useRef(false);

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

		async function onConnect() {
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
					alert(compatibility.message);
					socket.disconnect();
					return;
				}

				if (compatibility.severity === "warning") {
					console.warn(compatibility.message);
				}

				// Join appropriate room based on page
				if (isOverview) {
					socket.emit("overview:join");
					setConnected(true);
				} else if (roomId) {
					// Handle room:join response
					const handleJoinResponse = async (response: any) => {
						if (response.status !== "ok") {
							console.error("Failed to join room:", response.message);
							setInitializationError({
								message: "Failed to join room",
								details: response.message || "Server rejected the connection",
							});
							return;
						}

						// Reset chat unread count when entering a room
						useAppStore.getState().resetChatUnread();

						// Extract minimal data from socket response
						const { sessionId, step, frameCount, locked } = response;
						setSessionId(sessionId);

						// Track room visit for localStorage persistence
						setLastVisitedRoom(roomId);

						// Set minimal state from socket response
						setFrameCount(frameCount);
						setCurrentFrame(step);

						// Note: `locked` from join response is the ADMIN room lock,
						// not trajectory:meta lock. Admin lock is fetched via roomDetail query.
						// lockMetadata (trajectory:meta) is updated via lock:update socket events.
						// Do NOT set lockMetadata here - it caused stale "Someone: using this room".

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
							setGeometries(geometriesResponse.geometries || {});
							setGeometrySchemas(geometriesResponse.types?.schemas || {});
							setGeometryDefaults(geometriesResponse.types?.defaults || {});

							// Clear any previous error and mark as connected
							setInitializationError(null);
							setConnected(true);
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
								selectionsResponse,
								bookmarksResponse,
								frameSelResponse,
								lockStatusResponse,
							] = await Promise.all([
								getAllSelections(roomId),
								getAllBookmarks(roomId),
								getFrameSelection(roomId),
								getLockStatus(roomId, "trajectory:meta"),
							]);

							// Set selections
							setSelections(selectionsResponse.selections || {});
							setSelectionGroups(selectionsResponse.groups || {});
							setActiveSelectionGroup(selectionsResponse.activeGroup || null);

							// Set bookmarks
							setBookmarks(bookmarksResponse.bookmarks || {});

							// Set frame selection
							setFrameSelection(frameSelResponse.frameSelection);

							// Set trajectory:meta lock status (for edit/drawing modes)
							if (lockStatusResponse.locked) {
								setLockMetadata({
									locked: true,
									holder: lockStatusResponse.holder,
									userName: lockStatusResponse.metadata?.userName,
									msg: lockStatusResponse.metadata?.msg,
									timestamp: lockStatusResponse.metadata?.timestamp,
								});
							}
						} catch (error) {
							console.error("Error fetching secondary data:", error);
							// Non-critical - don't fail initialization
						}
					};

					// Emit room:join event
					socket.emit(
						"room:join",
						{ roomId, clientType: "frontend" },
						async (response: any) => {
							// Handle 404 - room doesn't exist, create it via REST API
							if (response.code === 404) {
								// Check for explicit template in URL, otherwise let server decide
								// Server will use: copyFrom > template > default_room > "empty"
								const urlTemplate = new URLSearchParams(
									window.location.search,
								).get("template");
								try {
									await createRoom({
										roomId,
										template: urlTemplate ?? undefined,
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
									"room:join",
									{ roomId, clientType: "frontend" },
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
			// Stop playback when another client (e.g., Python) changes the frame
			setPlaying(false);
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
				setBookmarks(response.bookmarks || {});
			},
			"bookmarks",
		);

		function onFramesInvalidate(data: any) {
			const { roomId, operation, affectedIndex, affectedFrom, affectedKeys } =
				data;

			if (affectedKeys) {
				console.warn(
					"Frame invalidation with affectedKeys is not yet implemented. Invalidating all frames.",
				);
			}

			queryClient.invalidateQueries({
				predicate: (query) => {
					// Query keys are expected to be in the format: ['frame', roomId, frameIndex]
					// or ['frame-keys', roomId, frameIndex]
					const [type, qRoomId, frameIndex] = query.queryKey;

					// Only invalidate frame and frame-keys queries for this room
					if ((type !== "frame" && type !== "frame-keys") || qRoomId !== roomId)
						return false;

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
						const keys = Object.keys(listResponse.geometries || {});

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
			"selections",
		);

		const onSelectionGroupsInvalidate = createInvalidateHandler(
			getAllSelections,
			(response) => {
				setSelectionGroups(response.groups);
				setActiveSelectionGroup(response.activeGroup);
			},
			"selection_groups",
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
			const {
				roomId: lockRoomId,
				target,
				action,
				holder,
				message,
				timestamp,
				sessionId,
			} = data;

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

			// Handle authentication errors (stale token in localStorage)
			// This happens when Redis is flushed but localStorage still has old token
			// Handles both "Client not registered" and "User not found" errors
			const isAuthError =
				err.message &&
				(err.message.includes("Client not registered") ||
					err.message.includes("User not found"));

			if (isAuthError) {
				// Prevent concurrent auth recovery - socket.io may emit multiple
				// connect_error events before the first recovery completes
				if (isRecoveringAuthRef.current) {
					console.log(
						"[onConnectError] Auth recovery already in progress, skipping",
					);
					return;
				}
				isRecoveringAuthRef.current = true;

				try {
					// Check if we've exceeded max retries
					if (authRetryCountRef.current >= MAX_AUTH_RETRIES) {
						console.error(
							`Max authentication retries (${MAX_AUTH_RETRIES}) exceeded. Please refresh the page.`,
						);
						// TODO: Show user-friendly error message
						return;
					}

					authRetryCountRef.current += 1;
					const retryAttempt = authRetryCountRef.current;

					// Calculate exponential backoff delay: 1s, 2s, 4s
					const delay = INITIAL_RETRY_DELAY * Math.pow(2, retryAttempt - 1);

					const { logout, login, getUsernameFromToken } = await import(
						"../utils/auth"
					);

					// Clear stale authentication data
					logout();

					// Wait for exponential backoff delay
					await new Promise((resolve) => setTimeout(resolve, delay));

					// Force a new login
					try {
						await login();

						// Update store with new username
						const newUsername = getUsernameFromToken();
						if (newUsername) {
							setUserName(newUsername);
						}

						// Reset retry count on successful login
						authRetryCountRef.current = 0;

						// Reconnect socket with new credentials
						socket.connect();
					} catch (loginError) {
						console.error(
							`Re-login failed (attempt ${retryAttempt}/${MAX_AUTH_RETRIES}):`,
							loginError,
						);
						// Will retry on next connect_error event if under max retries
					}
				} finally {
					isRecoveringAuthRef.current = false;
				}
			}
		}

		// Register event handlers FIRST before connecting
		socket.on("disconnect", onDisconnect);
		socket.on("connect", onConnect);
		socket.on("connect_error", onConnectError);
		socket.on("frame:update", onFrameUpdate);
		socket.on("active-camera:update", onActiveCameraUpdate);
		socket.on("invalidate", onInvalidate);
		socket.on("schema:invalidate", onSchemaInvalidate);
		socket.on("frame-selection:update", onFrameSelectionUpdate);
		socket.on("bookmarks:invalidate", onBookmarksInvalidate);
		socket.on("frames:invalidate", onFramesInvalidate);
		socket.on("chat:new", onChatMessageNew);
		socket.on("chat:update", onChatMessageUpdated);
		socket.on("geometry:invalidate", onGeometriesInvalidate);
		socket.on("figure:invalidate", onFiguresInvalidate);
		socket.on("selection:invalidate", onSelectionsInvalidate);
		socket.on("selection-groups:invalidate", onSelectionGroupsInvalidate);
		socket.on("room:update", onRoomUpdate);
		socket.on("room:delete", onRoomDelete);
		socket.on("lock:update", onLockUpdate);
		socket.on("progress:init", onProgressInitial);
		socket.on("progress:start", onProgressStarted);
		socket.on("progress:update", onProgressUpdate);
		socket.on("progress:complete", onProgressComplete);

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
				socket.emit("overview:leave");
			} else if (roomId) {
				socket.emit("room:leave", { roomId });
			}

			// Only clear sessionId if we're actually leaving a room
			// Compare the room context from when effect was created with current store values
			const currentState = useAppStore.getState();
			const isLeavingRoom =
				effectRoomId !== currentState.roomId ||
				effectIsOverview !== options.isOverview;

			if (isLeavingRoom) {
				setSessionId(null);
			}

			socket.off("connect", onConnect);
			socket.off("disconnect", onDisconnect);
			socket.off("frame:update", onFrameUpdate);
			socket.off("active-camera:update", onActiveCameraUpdate);
			socket.off("invalidate", onInvalidate);
			socket.off("schema:invalidate", onSchemaInvalidate);
			socket.off("frame-selection:update", onFrameSelectionUpdate);
			socket.off("bookmarks:invalidate", onBookmarksInvalidate);
			socket.off("frames:invalidate", onFramesInvalidate);
			socket.off("chat:new", onChatMessageNew);
			socket.off("chat:update", onChatMessageUpdated);
			socket.off("geometry:invalidate", onGeometriesInvalidate);
			socket.off("figure:invalidate", onFiguresInvalidate);
			socket.off("selection:invalidate", onSelectionsInvalidate);
			socket.off("selection-groups:invalidate", onSelectionGroupsInvalidate);
			socket.off("room:update", onRoomUpdate);
			socket.off("room:delete", onRoomDelete);
			socket.off("lock:update", onLockUpdate);
			socket.off("progress:init", onProgressInitial);
			socket.off("progress:start", onProgressStarted);
			socket.off("progress:update", onProgressUpdate);
			socket.off("progress:complete", onProgressComplete);

			// Disconnect socket when component unmounts to ensure clean reconnection
			socket.disconnect();
		};
	}, [
		roomId,
		// NOTE: userName intentionally NOT included - socket auth uses token from localStorage
		// Having userName here caused infinite re-runs when onConnectError created new users
		isOverview,
		setConnected,
		setInitializationError,
		setFrameCount,
		setCurrentFrame,
		setPlaying,
		queryClient,
		setBookmarks,
		setSelections,
		setSelectionGroups,
		setActiveSelectionGroup,
		setFrameSelection,
		setSessionId,
		setGeometries,
		setGeometrySchemas,
		setGeometryDefaults,
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
