import { useQueryClient } from "@tanstack/react-query";
import { useEffect, useRef } from "react";
import { useAppStore } from "../store";
import { prefetchFrame } from "./useFrameBatch";

const PREFETCH_AHEAD = 3;

/**
 * Batch-prefetches upcoming frames during playback.
 *
 * Uses a Zustand store subscription (not a React effect dependency on
 * currentFrame) so that a single AbortController persists across frame
 * ticks. Prefetches are NOT aborted when the frame advances — only when
 * playback stops or the component unmounts.
 *
 * For each upcoming frame, collects all needed keys from active geometries
 * and delegates to prefetchFrame() which filters already-cached keys and
 * fetches the remainder in ONE HTTP request per frame.
 */
export function useFramePrefetcher(): void {
	const queryClient = useQueryClient();
	const playing = useAppStore((s) => s.playing);
	const roomId = useAppStore((s) => s.roomId);

	const inFlightRef = useRef(new Set<number>());

	useEffect(() => {
		if (!playing || !roomId) return;

		const controller = new AbortController();

		function prefetchAround(currentFrame: number) {
			const state = useAppStore.getState();
			const rid = state.roomId;
			if (!rid) return;

			const allKeys = [...new Set(Object.values(state.neededFrameKeys).flat())];
			if (allKeys.length === 0) return;

			const skip = state.skipFrames;
			const fc = state.frameCount;
			const sel = state.frameSelectionEnabled ? state.frame_selection : null;

			// Compute upcoming frames
			const upcoming: number[] = [];
			if (sel && sel.length > 0) {
				const sorted = [...sel].sort((a, b) => a - b);
				const idx = sorted.indexOf(currentFrame);
				if (idx === -1) return;
				for (let i = 1; i <= PREFETCH_AHEAD; i++) {
					const nextIdx = idx + i * skip;
					if (nextIdx >= sorted.length) break;
					upcoming.push(sorted[nextIdx]);
				}
			} else {
				for (let i = 1; i <= PREFETCH_AHEAD; i++) {
					const next = currentFrame + i * skip;
					if (next >= fc) break;
					upcoming.push(next);
				}
			}

			for (const frame of upcoming) {
				if (inFlightRef.current.has(frame)) continue;
				inFlightRef.current.add(frame);
				prefetchFrame(
					rid,
					frame,
					allKeys,
					queryClient,
					controller.signal,
				).finally(() => inFlightRef.current.delete(frame));
			}
		}

		// Prefetch for current position immediately
		prefetchAround(useAppStore.getState().currentFrame);

		// Subscribe to frame changes — fires on every store update
		const unsubscribe = useAppStore.subscribe((state, prevState) => {
			if (state.currentFrame !== prevState.currentFrame) {
				prefetchAround(state.currentFrame);
			}
		});

		return () => {
			unsubscribe();
			controller.abort();
			inFlightRef.current.clear();
		};
	}, [playing, roomId, queryClient]);
}
