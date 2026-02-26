import type { QueryClient } from "@tanstack/react-query";
import { type FrameResponse, getFrames } from "../myapi/client";

/** Max retries for retriable 404s (provider computing frames). */
const MAX_RETRIES = 15;

interface PendingBatch {
	keys: Set<string>;
	promise: Promise<FrameResponse | null>;
	controller: AbortController;
}

/** In-flight batches keyed by "roomId:frame". */
const batchPromises = new Map<string, PendingBatch>();

/**
 * Check if an error is a retriable 404 (has Retry-After header).
 * Returns the delay in ms, or null if not retriable.
 */
function getRetryDelay(error: any): number | null {
	if (error?.response?.status !== 404) return null;
	const retryAfter = error.response.headers?.["retry-after"];
	if (!retryAfter) return null;
	const seconds = Number.parseInt(retryAfter, 10);
	return Number.isFinite(seconds) ? seconds * 1000 : 2000;
}

/** Sleep that respects AbortSignal. */
function abortableSleep(ms: number, signal: AbortSignal): Promise<void> {
	return new Promise((resolve, reject) => {
		if (signal.aborted) {
			reject(new DOMException("Aborted", "AbortError"));
			return;
		}
		const onAbort = () => {
			clearTimeout(timer);
			reject(new DOMException("Aborted", "AbortError"));
		};
		const timer = setTimeout(() => {
			signal.removeEventListener("abort", onAbort);
			resolve();
		}, ms);
		signal.addEventListener("abort", onAbort, { once: true });
	});
}

/**
 * Fetch frames with automatic retry for provider-computed frames.
 *
 * When the server returns 404 with Retry-After (provider still computing),
 * retries up to MAX_RETRIES times, respecting the server's delay.
 */
async function fetchWithRetry(
	roomId: string,
	frame: number,
	keys: string[],
	signal: AbortSignal,
): Promise<FrameResponse | null> {
	for (let attempt = 0; attempt <= MAX_RETRIES; attempt++) {
		try {
			return await getFrames(roomId, frame, keys, signal);
		} catch (error: any) {
			const delay = getRetryDelay(error);
			if (delay === null || attempt === MAX_RETRIES) {
				throw error;
			}
			await abortableSleep(delay, signal);
		}
	}
	return null; // unreachable, but satisfies TS
}

/**
 * Coalesces per-key frame queries into a single HTTP request.
 *
 * The HTTP request is deferred via setTimeout(0), which fires AFTER
 * React's useEffect flush. This ensures all components' queryFn calls
 * (triggered from useEffect) register their keys before the fetch starts.
 * Results are distributed to individual TanStack Query cache entries
 * via setQueryData.
 *
 * When the batch fires it is "sealed" — removed from the map so any
 * later-arriving queryFn creates a fresh batch instead of joining a
 * stale one.
 *
 * When a new frame is requested for the same room, any in-flight
 * batch for a previous frame is aborted (fast scrubbing support).
 *
 * Retriable 404s (provider-computed frames with Retry-After) are
 * retried internally — the promise stays pending so TanStack Query's
 * isFetching remains true.
 */
export function getFrameBatched(
	roomId: string,
	frame: number,
	key: string,
	queryClient: QueryClient,
): Promise<FrameResponse | null> {
	const batchKey = `${roomId}:${frame}`;

	// Abort all stale batches for this room (fast scrubbing support)
	const prefix = `${roomId}:`;
	for (const [bk, pending] of batchPromises) {
		if (bk.startsWith(prefix) && bk !== batchKey) {
			pending.controller.abort();
		}
	}

	let batch = batchPromises.get(batchKey);
	if (!batch) {
		const keys = new Set<string>();
		const controller = new AbortController();

		// Defer fetch so queryFn calls in the same event-loop turn
		// can register their keys. setTimeout(0) fires AFTER React's
		// useEffect flush, ensuring all components' queries join.
		const promise = new Promise<FrameResponse | null>((resolve, reject) => {
			setTimeout(() => {
				// Seal: late arrivals will create a new batch
				batchPromises.delete(batchKey);

				const toFetch = [...keys];

				fetchWithRetry(roomId, frame, toFetch, controller.signal)
					.then((data) => {
						for (const k of toFetch) {
							queryClient.setQueryData(
								["frame", roomId, frame, k],
								data?.[k] !== undefined ? { [k]: data[k] } : null,
							);
						}
						resolve(data);
					})
					.catch((error) => {
						reject(error);
					});
			}, 0);
		});

		batch = { keys, promise, controller };
		batchPromises.set(batchKey, batch);
	}

	// Add caller's key — will be included when the deferred timeout fires
	batch.keys.add(key);

	// All callers share the same promise, each extracts their key
	return batch.promise.then((data) =>
		data?.[key] !== undefined ? ({ [key]: data[key] } as FrameResponse) : null,
	);
}

/**
 * Prefetch a frame's data into the TanStack Query cache.
 *
 * Used by the prefetcher for ahead-of-current frames during playback.
 * Does NOT abort stale frames (prefetches for multiple frames can coexist).
 * Skips keys already in cache.
 */
export async function prefetchFrame(
	roomId: string,
	frame: number,
	allKeys: string[],
	queryClient: QueryClient,
	signal?: AbortSignal,
): Promise<void> {
	const uncached = allKeys.filter(
		(k) => queryClient.getQueryData(["frame", roomId, frame, k]) === undefined,
	);
	if (uncached.length === 0) return;

	try {
		const data = await getFrames(roomId, frame, uncached, signal);
		for (const key of uncached) {
			queryClient.setQueryData(
				["frame", roomId, frame, key],
				data?.[key] !== undefined ? { [key]: data[key] } : null,
			);
		}
	} catch {
		// Silently ignore prefetch errors (abort, network, 404, etc.)
	}
}
