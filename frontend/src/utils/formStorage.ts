/**
 * Form data persistence utilities using localStorage.
 *
 * Stores form data per room/user/job with a 7-day TTL.
 */

const FORM_STORAGE_PREFIX = "zndraw:form";
const FORM_TTL_MS = 7 * 24 * 60 * 60 * 1000; // 7 days

interface StoredFormData {
	data: Record<string, unknown>;
	timestamp: number;
}

const getStorageKey = (roomId: string, userName: string, jobName: string) =>
	`${FORM_STORAGE_PREFIX}:${roomId}:${userName}:${jobName}`;

/** Save form data to localStorage with timestamp. */
export const saveFormData = (
	roomId: string,
	userName: string,
	jobName: string,
	data: Record<string, unknown>,
): void => {
	const stored: StoredFormData = { data, timestamp: Date.now() };
	localStorage.setItem(
		getStorageKey(roomId, userName, jobName),
		JSON.stringify(stored),
	);
};

/** Load form data from localStorage, returning null if expired or not found. */
export const loadFormData = (
	roomId: string,
	userName: string,
	jobName: string,
): Record<string, unknown> | null => {
	const raw = localStorage.getItem(getStorageKey(roomId, userName, jobName));
	if (!raw) return null;

	try {
		const stored: StoredFormData = JSON.parse(raw);

		if (Date.now() - stored.timestamp > FORM_TTL_MS) {
			localStorage.removeItem(getStorageKey(roomId, userName, jobName));
			return null;
		}

		return stored.data;
	} catch {
		localStorage.removeItem(getStorageKey(roomId, userName, jobName));
		return null;
	}
};

/** Clear form data for a specific job. */
export const clearFormData = (
	roomId: string,
	userName: string,
	jobName: string,
): void => {
	localStorage.removeItem(getStorageKey(roomId, userName, jobName));
};
