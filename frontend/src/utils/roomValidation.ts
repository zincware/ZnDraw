/**
 * Validate a room ID for creation/duplication.
 * Returns an error message string, or null if valid.
 */
export function validateRoomId(
	roomIdToCheck: string,
	existingIds: string[],
): string | null {
	if (!roomIdToCheck) {
		return null; // Empty means auto-generate, which is valid
	}

	if (existingIds.includes(roomIdToCheck)) {
		return "A room with this ID already exists";
	}

	if (!/^[a-zA-Z0-9_-]+$/.test(roomIdToCheck)) {
		return "Room ID can only contain letters, numbers, hyphens, and underscores";
	}

	return null;
}
