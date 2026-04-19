import { isAxiosError } from "axios";

/**
 * Coerce a thrown error (typically from axios) into a flat error string.
 *
 * Handles FastAPI's three common detail shapes:
 *   - Plain string: returned as-is.
 *   - Validation array: joined into a single semicolon-separated string.
 *   - Missing/unknown: falls back to the axios message, then the fallback.
 */
export function extractDetail(err: unknown, fallback: string): string {
	if (!isAxiosError(err)) return fallback;
	const detail = err.response?.data?.detail;
	if (typeof detail === "string") return detail;
	if (Array.isArray(detail)) {
		return detail
			.map((d) => (typeof d === "object" && d?.msg ? String(d.msg) : String(d)))
			.join("; ");
	}
	return err.message ?? fallback;
}
