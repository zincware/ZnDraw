import type { Room } from "../../myapi/client";
import { useRoomsStore } from "../../roomsStore";
import { useAppStore } from "../../store";
import type { HandlerContext } from "./types";

// --- Typed event interfaces ---

export interface RoomUpdateEvent {
	id: string;
	frame_count?: number | null;
	locked?: boolean | null;
	[key: string]: unknown;
}

export interface RoomDeleteEvent {
	room_id: string;
}

export interface LockUpdateEvent {
	action: "acquired" | "refreshed" | "released";
	user_id?: string | null;
	sid?: string;
	msg?: string | null;
	ttl?: number;
}

export interface ProgressStartedEvent {
	progress_id: string;
	description: string;
	unit?: string;
}

export interface ProgressUpdateEvent {
	progress_id: string;
	n?: number;
	total?: number | null;
	elapsed?: number;
}

export interface ProgressCompleteEvent {
	progress_id: string;
}

// --- Factory ---

export function createRoomHandlers(ctx: HandlerContext) {
	function onRoomUpdate(data: RoomUpdateEvent) {
		console.debug("[RoomUpdate] received:", {
			data,
			currentRoomId: ctx.roomId,
		});

		// Update in-room state if this event is for the current room
		if (data.id === ctx.roomId) {
			if (data.frame_count != null) {
				ctx.setFrameCount(data.frame_count);
			}
			if (data.locked != null) {
				console.debug("[RoomUpdate] lock change:", data.locked);
				ctx.setSuperuserLock(data.locked);
			}
		}

		// Upsert into rooms store (full snapshot -- always safe to overwrite)
		// Server sends complete Room objects; the event type is permissive for
		// partial reads above, so cast to Room for the store API.
		useRoomsStore.getState().setRoom(data.id, data as Room);
	}

	function onRoomDelete(data: RoomDeleteEvent) {
		const { room_id: deletedRoomId } = data;
		useRoomsStore.getState().removeRoom(deletedRoomId);
	}

	function onLockUpdate(data: LockUpdateEvent) {
		const { action, user_id, sid, msg, ttl } = data;
		const mySessionId = useAppStore.getState().sessionId;

		if (action === "acquired" || action === "refreshed") {
			// If this is our own session, lockSlice already set the state
			if (sid === mySessionId) return;
			ctx.setUserLock(user_id ?? null, msg ?? null);
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
				ctx.setUserLock(null, null);
			}
			useAppStore.getState().stopLockExpiryTimer();
		}
	}

	function onProgressStarted(data: ProgressStartedEvent) {
		ctx.addProgressTracker({
			progress_id: data.progress_id,
			description: data.description,
			n: 0,
			total: null,
			elapsed: 0,
			unit: data.unit ?? "it",
		});
	}

	function onProgressUpdate(data: ProgressUpdateEvent) {
		ctx.updateProgressTracker(data);
	}

	function onProgressComplete(data: ProgressCompleteEvent) {
		ctx.removeProgressTracker(data.progress_id);
	}

	return {
		onRoomUpdate,
		onRoomDelete,
		onLockUpdate,
		onProgressStarted,
		onProgressUpdate,
		onProgressComplete,
	};
}
