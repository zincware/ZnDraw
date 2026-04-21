import {
	createRoom,
	type GeometryData,
	getEditLockStatus,
	getFrameSelection,
	getGlobalSettings,
	getServerVersion,
	listGeometries,
	listSelectionGroups,
	getAllBookmarks,
} from "../../myapi/client";
import { connectWithAuth, socket } from "../../socket";
import { useAppStore } from "../../store";
import { logout } from "../../utils/auth";
import { setLastVisitedRoom } from "../../utils/roomTracking";
import {
	checkVersionCompatibility,
	getClientVersion,
} from "../../utils/versionCompatibility";
import type { HandlerContext } from "./types";

// --- Typed event interfaces ---

export interface RoomJoinResponse {
	session_id: string;
	step: number;
	frame_count: number;
	camera_key: string | null;
	locked: boolean;
	progress_trackers?: Record<string, import("../../store").Progress>;
}

export interface RoomJoinError {
	status: number;
	detail?: string;
}

// --- Factory ---

const INITIAL_RETRY_DELAY = 1000;
const MAX_RETRY_DELAY = 30000;

export function createConnectionHandlers(ctx: HandlerContext) {
	let retryDelay = INITIAL_RETRY_DELAY;

	/**
	 * Process a successful room_join response: set session state,
	 * fetch render-blocking geometry data, then fetch secondary data.
	 */
	async function handleRoomJoin(response: RoomJoinResponse) {
		// Reset chat unread count when entering a room
		useAppStore.getState().resetChatUnread();

		// Extract minimal data from socket response (backend uses snake_case)
		const sessionId = response.session_id;
		const step = response.step;
		const frameCount = response.frame_count;
		ctx.setSessionId(sessionId);
		ctx.setCameraKey(response.camera_key ?? null);

		// Track room visit for localStorage persistence
		setLastVisitedRoom(ctx.roomId!);

		// Set minimal state from socket response
		ctx.setFrameCount(frameCount);
		ctx.setCurrentFrame(step);

		// Superuser lock from join response (SQL room.locked)
		console.debug("[RoomJoin] locked:", response.locked);
		ctx.setSuperuserLock(response.locked ?? false);

		// Set progress trackers from join response
		if (response.progress_trackers) {
			ctx.setProgressTrackers(response.progress_trackers);
		}

		// === CRITICAL: Fetch render-blocking data via REST ===
		// Both calls run in parallel - we need both to render the scene
		try {
			const geometriesResponse = await listGeometries(ctx.roomId!);

			// Set geometries and type metadata (schemas + defaults)
			const geos = geometriesResponse.items || {};
			ctx.setGeometries(geos);
			ctx.setGeometrySchemas(geometriesResponse.types?.schemas || {});
			ctx.setGeometryDefaults(geometriesResponse.types?.defaults || {});

			// Extract per-geometry selections (cameras don't have selections)
			const selectionsFromGeos: Record<string, number[]> = {};
			for (const [key, geo] of Object.entries(geos) as [
				string,
				GeometryData,
			][]) {
				if (geo.type === "Camera") continue;
				selectionsFromGeos[key] = geo.selection ?? [];
			}
			ctx.setSelections(selectionsFromGeos);

			// Clear any previous error and mark as connected
			ctx.setInitializationError(null);
			ctx.setConnected(true);
		} catch (error) {
			console.error("Error fetching critical data:", error);
			ctx.setInitializationError({
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
				listSelectionGroups(ctx.roomId!),
				getAllBookmarks(ctx.roomId!),
				getFrameSelection(ctx.roomId!),
				getEditLockStatus(ctx.roomId!),
			]);

			// Set selection groups
			ctx.setSelectionGroups(groupsResponse.items || {});

			// Set bookmarks
			ctx.setBookmarks(bookmarksResponse.items || {});

			// Set frame selection
			ctx.setFrameSelection(frameSelResponse.frame_selection);

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
					ctx.setUserLock(
						editLockResponse.user_id ?? null,
						editLockResponse.msg ?? null,
					);
					if (editLockResponse.ttl) {
						useAppStore.getState().startLockExpiryTimer(editLockResponse.ttl);
					}
				}
			}
		} catch (error) {
			console.error("Error fetching secondary data:", error);
			// Non-critical - don't fail initialization
		}
	}

	async function onConnect() {
		retryDelay = INITIAL_RETRY_DELAY;
		// Clear any previous initialization error on new connection
		ctx.setInitializationError(null);

		try {
			// Fetch server version and global settings
			const { version: serverVersion } = await getServerVersion();
			ctx.setServerVersion(serverVersion);

			const globalSettings = await getGlobalSettings();
			ctx.setGlobalSettings(globalSettings);

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

			// Join the specific room if one is set; otherwise connect is
			// sufficient — rooms:feed is auto-joined server-side.
			if (ctx.roomId) {
				socket.emit(
					"room_join",
					{ room_id: ctx.roomId, client_type: "frontend" },
					async (response: RoomJoinResponse | RoomJoinError) => {
						// Handle 404 - room doesn't exist, create it via REST API
						if ("status" in response && response.status === 404) {
							const urlCopyFrom = new URLSearchParams(
								window.location.search,
							).get("copy_from");
							try {
								await createRoom({
									room_id: ctx.roomId!,
									copy_from: urlCopyFrom ?? undefined,
								});
							} catch (error: any) {
								// 409 Conflict = room created by another client, continue
								if (error.response?.status !== 409) {
									console.error("Failed to create room:", error);
									ctx.setInitializationError({
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
								{ room_id: ctx.roomId!, client_type: "frontend" },
								(retryResponse: RoomJoinResponse | RoomJoinError) => {
									if ("status" in retryResponse) {
										console.error(
											"Failed to join room:",
											retryResponse,
										);
										ctx.setInitializationError({
											message: "Failed to join room",
											details:
												retryResponse.detail ||
												"Server rejected the connection",
										});
										return;
									}
									handleRoomJoin(retryResponse);
								},
							);
							return;
						}

						if ("status" in response) {
							console.error("Failed to join room:", response);
							ctx.setInitializationError({
								message: "Failed to join room",
								details:
									response.detail || "Server rejected the connection",
							});
							return;
						}

						handleRoomJoin(response);
					},
				);
			} else {
				ctx.setConnected(true);
			}
		} catch (error) {
			console.error("Error checking version compatibility:", error);
			// Still connect even if version check fails
			ctx.setConnected(true);
		}
	}

	function onDisconnect() {
		ctx.setConnected(false);
		// NOTE: Do NOT clear sessionId here - onDisconnect fires during temporary
		// reconnects (e.g., when userName changes). SessionId is only cleared
		// in the cleanup function when actually leaving a room.
	}

	async function onConnectError(err: Error) {
		if (ctx.isCancelled()) return;
		console.error("Socket connection error:", err.message);

		// socket.active === true means temporary network failure -- auto-reconnects.
		if (socket.active) return;

		// Server rejected the connection (auth error).
		// Clear stale token and retry with backoff -- guest auth should
		// always succeed unless the server is down.
		logout();
		await new Promise((r) => setTimeout(r, retryDelay));
		retryDelay = Math.min(retryDelay * 2, MAX_RETRY_DELAY);

		if (ctx.isCancelled()) return;
		try {
			const { user } = await connectWithAuth();
			ctx.setUser(user);
		} catch {
			// connectWithAuth failed (server down?) -- next connect_error will retry
		}
	}

	return { onConnect, onDisconnect, onConnectError };
}
