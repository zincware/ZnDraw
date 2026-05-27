import type { Room } from "../../myapi/client";
import { useRoomsStore } from "../../roomsStore";
import { useAppStore } from "../../store";
import type { HandlerContext } from "./types";

// --- Typed event interfaces ---

export interface RoomUpdateEvent {
	room_id: string;
	frame_count?: number | null;
	[key: string]: unknown;
}

export interface RoomDeleteEvent {
	room_id: string;
}

export interface LockUpdateEvent {
	action: "acquired" | "refreshed" | "released";
	user_id?: string | null;
	display_name?: string | null;
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
		const composedRoomId =
			ctx.ownerId && ctx.roomName
				? `${ctx.ownerId}/${ctx.roomName}`
				: undefined;
		console.debug("[RoomUpdate] received:", {
			data,
			currentRoomId: composedRoomId,
		});

		// Update in-room state if this event is for the current room
		if (data.room_id === composedRoomId) {
			if (data.frame_count != null) {
				ctx.setFrameCount(data.frame_count);
			}
		}

		// Upsert into rooms store (full snapshot -- always safe to overwrite)
		useRoomsStore.getState().setRoom(data.room_id, data as unknown as Room);
	}

	function onRoomDelete(data: RoomDeleteEvent) {
		const { room_id: deletedRoomId } = data;
		useRoomsStore.getState().removeRoom(deletedRoomId);
	}

	function onLockUpdate(data: LockUpdateEvent) {
		const { action, display_name, sid, msg, ttl } = data;
		const mySessionId = useAppStore.getState().sessionId;

		if (action === "acquired" || action === "refreshed") {
			// If this is our own session, lockSlice already set the state
			if (sid === mySessionId) return;
			ctx.setUserLock(display_name ?? null, msg ?? null);
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
