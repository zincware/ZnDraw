/**
 * Simple room tracking for localStorage-based room persistence.
 *
 * Stores only the last visited room ID.
 */

const LAST_ROOM_KEY = "zndraw_last_room";

/**
 * Get the last visited room ID from localStorage.
 *
 * Returns
 * -------
 * string | null
 *     The room ID or null if not set.
 */
export function getLastVisitedRoom(): string | null {
	return localStorage.getItem(LAST_ROOM_KEY);
}

/**
 * Set the last visited room ID in localStorage.
 *
 * Parameters
 * ----------
 * roomId : string
 *     The room ID to store.
 */
export function setLastVisitedRoom(roomId: string): void {
	localStorage.setItem(LAST_ROOM_KEY, roomId);
}
