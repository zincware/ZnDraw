/**
 * Factory for creating consistent invalidate handlers.
 *
 * Fetches data from the server and updates the store.
 * Uses a `getRoomId` getter so the handler always reads the current roomId
 * at event-fire time (same behavior as the original closure).
 *
 * If the incoming event payload carries a `room_address` and it does not
 * match the current room, the handler is a no-op. This prevents stale
 * events from a previously-joined room mutating the current room's cache
 * during navigation races.
 */
export function createInvalidateHandler<T>(
	fetchFn: (roomId: string) => Promise<T>,
	updateStoreFn: (data: T) => void,
	eventName: string,
	getRoomId: () => string | undefined,
): (data: unknown) => Promise<void> {
	return async (data) => {
		const roomId = getRoomId();
		if (!roomId) return;
		if (data && typeof data === "object" && "room_address" in data) {
			const evtAddr = (data as { room_address?: unknown }).room_address;
			if (typeof evtAddr === "string" && evtAddr !== roomId) return;
		}
		try {
			const response = await fetchFn(roomId);
			updateStoreFn(response);
		} catch (error) {
			console.error(`Error fetching ${eventName}:`, error);
		}
	};
}
