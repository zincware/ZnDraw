/**
 * Factory for creating consistent invalidate handlers.
 *
 * Fetches data from the server and updates the store.
 * Uses a `getRoomId` getter so the handler always reads the current roomId
 * at event-fire time (same behavior as the original closure).
 */
export function createInvalidateHandler<T>(
	fetchFn: (roomId: string) => Promise<T>,
	updateStoreFn: (data: T) => void,
	eventName: string,
	getRoomId: () => string | undefined,
): (data: unknown) => Promise<void> {
	return async () => {
		const roomId = getRoomId();
		if (!roomId) return;
		try {
			const response = await fetchFn(roomId);
			updateStoreFn(response);
		} catch (error) {
			console.error(`Error fetching ${eventName}:`, error);
		}
	};
}
