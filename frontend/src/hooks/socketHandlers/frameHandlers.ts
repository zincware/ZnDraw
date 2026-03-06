import { useAppStore } from "../../store";
import type { HandlerContext } from "./types";

// --- Typed event interfaces ---

export interface FrameUpdateEvent {
	frame: number;
}

export interface FramesInvalidateEvent {
	room_id: string;
	action: "add" | "delete" | "modify" | "clear";
	indices?: number[];
	count?: number | null;
	reason?: string | null;
}

export interface FrameSelectionUpdateEvent {
	indices: number[] | null;
}

// --- Factory ---

export function createFrameHandlers(ctx: HandlerContext) {
	function onFrameUpdate(data: FrameUpdateEvent) {
		const { frame } = data;
		// During local playback, ignore server echoes of our own step updates.
		// Without session_id in the event we cannot distinguish self from remote,
		// so all frame_update events are ignored while playing.
		if (useAppStore.getState().playing) return;
		ctx.setCurrentFrame(frame);
	}

	function onFramesInvalidate(data: FramesInvalidateEvent) {
		const { room_id: eventRoomId, action, indices, count, reason } = data;

		// Update frameCount if provided (new total frame count)
		if (count != null) {
			ctx.setFrameCount(count);
		}

		// Notify user when a mounted source disconnects
		if (action === "clear" && reason === "provider_disconnected") {
			useAppStore.getState().showSnackbar("Source disconnected", "warning");
		}

		// Invalidate React Query cache based on action
		ctx.queryClient.invalidateQueries({
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
					// No specific indices -- invalidate all frame queries
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

	function onFrameSelectionUpdate(data: FrameSelectionUpdateEvent) {
		ctx.setFrameSelection(data.indices || null);
	}

	return { onFrameUpdate, onFramesInvalidate, onFrameSelectionUpdate };
}
