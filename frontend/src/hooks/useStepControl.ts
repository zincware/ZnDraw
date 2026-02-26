import { throttle } from "lodash";
import { useCallback, useEffect, useMemo, useState } from "react";
import { updateStep } from "../myapi/client";
import { socket } from "../socket";
import { useAppStore } from "../store";

interface StepControlHook {
	setStep: (frame: number) => void;
	remoteLocked: boolean;
}

/**
 * Unified hook for step/frame control.
 *
 * Features:
 * - Single setStep() function for all use cases (clicks, keyboard, slider, playback)
 * - Throttled server updates (100ms) with trailing edge to ensure final value synced
 * - No automatic lock acquisition - server handles via @check_lock decorator
 * - remoteLocked state indicates when another user has the step lock
 *
 * Usage:
 *   const { setStep, remoteLocked } = useStepControl();
 *   setStep(42); // Updates local state, syncs to server
 *   <Slider disabled={remoteLocked} />
 *
 * remoteLocked logic:
 * - false if no one has the lock
 * - false if WE have the lock
 * - true ONLY if someone else has the lock
 */
export const useStepControl = (): StepControlHook => {
	const [remoteLocked, setRemoteLocked] = useState<boolean>(false);

	const roomId = useAppStore((state) => state.roomId);
	const setCurrentFrame = useAppStore((state) => state.setCurrentFrame);

	// Throttled server update (100ms throttle with trailing edge)
	const throttledUpdate = useMemo(
		() =>
			throttle(
				async (frame: number) => {
					if (!roomId) return;

					try {
						await updateStep(roomId, frame);
					} catch (error: any) {
						if (error.response?.status === 423) {
							// Someone else has the step lock
							useAppStore
								.getState()
								.showSnackbar(
									"Step update blocked - another user is controlling playback",
									"warning",
								);
						} else {
							console.error("Failed to update step:", error);
						}
					}
				},
				100, // Max 10 updates/second
				{ leading: true, trailing: true }, // Execute first call immediately + last call after throttle
			),
		[roomId],
	);

	// Universal step setter
	const setStep = useCallback(
		(frame: number) => {
			if (!roomId) {
				console.error("Cannot set step: no room ID");
				return;
			}

			// 1. Update local state immediately for responsive UI
			setCurrentFrame(frame);

			// 2. Update server (throttled) - no lock needed, server checks via @check_lock
			throttledUpdate(frame);
		},
		[roomId, setCurrentFrame, throttledUpdate],
	);

	// Track remote lock state from room-level edit lock
	// If someone else holds the edit lock, step updates will be rejected (423)
	useEffect(() => {
		const onLockUpdate = (data: any) => {
			const { action } = data;
			// Check if we hold the lock ourselves
			const hasOurLock = useAppStore.getState().hasLock();

			if (action === "released") {
				setRemoteLocked(false);
			} else if (!hasOurLock) {
				// Someone else acquired/refreshed the lock
				setRemoteLocked(true);
			}
		};

		socket.on("lock_update", onLockUpdate);

		return () => {
			socket.off("lock_update", onLockUpdate);
		};
	}, []);

	// Cleanup on unmount
	useEffect(() => {
		return () => {
			throttledUpdate.cancel();
		};
	}, [throttledUpdate]);

	return {
		setStep,
		remoteLocked,
	};
};
