import { useCallback, useEffect } from "react";
import { useAppStore } from "../store";
import { isInputFocused } from "../utils/keyboard";
import { useStepControl } from "./useStepControl";

export const useKeyboardShortcuts = () => {
	const currentFrame = useAppStore((state) => state.currentFrame);
	const frameCount = useAppStore((state) => state.frameCount);
	const frame_selection = useAppStore((state) => state.frame_selection);
	const frameSelectionEnabled = useAppStore(
		(state) => state.frameSelectionEnabled,
	);
	const playing = useAppStore((state) => state.playing);
	const setPlaying = useAppStore((state) => state.setPlaying);
	const skipFrames = useAppStore((state) => state.skipFrames);
	const synchronizedMode = useAppStore((state) => state.synchronizedMode);
	const getIsFetching = useAppStore((state) => state.getIsFetching);
	const addBookmark = useAppStore((state) => state.addBookmark);
	const bookmarks = useAppStore((state) => state.bookmarks);
	const updateFrameSelection = useAppStore(
		(state) => state.updateFrameSelection,
	);
	const setFps = useAppStore((state) => state.setFps);
	const setLastFrameChangeTime = useAppStore(
		(state) => state.setLastFrameChangeTime,
	);

	const { setStep } = useStepControl();

	const getNavigableFrames = useCallback((): number[] => {
		if (
			frameSelectionEnabled &&
			frame_selection &&
			frame_selection.length > 0
		) {
			return [...frame_selection].sort((a, b) => a - b);
		}
		return Array.from({ length: frameCount }, (_, i) => i);
	}, [frameSelectionEnabled, frame_selection, frameCount]);

	const handlePreviousFrame = useCallback(() => {
		const { playing, currentFrame } = useAppStore.getState();
		if (playing) setPlaying(false);

		const navigableFrames = getNavigableFrames();
		if (navigableFrames.length === 0) return;

		const currentIndex = navigableFrames.indexOf(currentFrame);
		if (currentIndex > 0) {
			setStep(navigableFrames[currentIndex - 1]);
		} else {
			setStep(navigableFrames[navigableFrames.length - 1]);
		}
	}, [getNavigableFrames, setStep, setPlaying]);

	const handleNextFrame = useCallback(() => {
		const { playing, currentFrame } = useAppStore.getState();
		if (playing) setPlaying(false);

		const navigableFrames = getNavigableFrames();
		if (navigableFrames.length === 0) return;

		const currentIndex = navigableFrames.indexOf(currentFrame);
		if (currentIndex < navigableFrames.length - 1) {
			setStep(navigableFrames[currentIndex + 1]);
		} else {
			setStep(navigableFrames[0]);
		}
	}, [getNavigableFrames, setStep, setPlaying]);

	const handleTogglePlayback = useCallback(() => {
		const { playing, currentFrame } = useAppStore.getState();
		const navigableFrames = getNavigableFrames();
		const isAtLastFrame =
			currentFrame === navigableFrames[navigableFrames.length - 1];

		if (!playing) {
			if (isAtLastFrame) {
				setStep(navigableFrames[0]);
			}
			setPlaying(true);
		} else {
			setPlaying(false);
		}
	}, [getNavigableFrames, setStep, setPlaying]);

	const handleAddBookmark = useCallback(() => {
		addBookmark(useAppStore.getState().currentFrame);
	}, [addBookmark]);

	const handlePreviousBookmark = useCallback(() => {
		if (!bookmarks) return;
		const { playing, currentFrame } = useAppStore.getState();
		if (playing) setPlaying(false);

		const bookmarkFrames = Object.keys(bookmarks)
			.map(Number)
			.sort((a, b) => a - b);
		if (bookmarkFrames.length === 0) return;

		const prevBookmark = [...bookmarkFrames]
			.reverse()
			.find((frame) => frame < currentFrame);
		if (prevBookmark !== undefined) {
			setStep(prevBookmark);
		} else {
			setStep(bookmarkFrames[0]);
		}
	}, [bookmarks, setStep, setPlaying]);

	const handleNextBookmark = useCallback(() => {
		if (!bookmarks) return;
		const { playing, currentFrame } = useAppStore.getState();
		if (playing) setPlaying(false);

		const bookmarkFrames = Object.keys(bookmarks)
			.map(Number)
			.sort((a, b) => a - b);
		if (bookmarkFrames.length === 0) return;

		const nextBookmark = bookmarkFrames.find((frame) => frame > currentFrame);
		if (nextBookmark !== undefined) {
			setStep(nextBookmark);
		} else {
			setStep(bookmarkFrames[0]);
		}
	}, [bookmarks, setStep, setPlaying]);

	// Handle automatic frame advancement during playback
	useEffect(() => {
		if (!playing) return;

		const intervalId = setInterval(() => {
			if (synchronizedMode && getIsFetching()) return;

			const frame = useAppStore.getState().currentFrame;
			const navigableFrames = getNavigableFrames();
			const currentIndex = navigableFrames.indexOf(frame);
			const lastIndex = navigableFrames.length - 1;
			const nextIndex = currentIndex + skipFrames;

			if (nextIndex < lastIndex) {
				setStep(navigableFrames[nextIndex]);
			} else if (currentIndex < lastIndex) {
				setStep(navigableFrames[lastIndex]);
			} else {
				setPlaying(false);
			}
		}, 33);

		return () => clearInterval(intervalId);
	}, [
		playing,
		getNavigableFrames,
		setStep,
		skipFrames,
		setPlaying,
		synchronizedMode,
		getIsFetching,
	]);

	// Track FPS during playback
	useEffect(() => {
		if (!playing) {
			// Reset FPS when not playing
			setFps(null);
			setLastFrameChangeTime(null);
			return;
		}

		// Calculate FPS based on time since last frame change
		const now = performance.now();
		const lastFrameTime = useAppStore.getState().lastFrameChangeTime;

		if (lastFrameTime !== null) {
			const delta = now - lastFrameTime;
			if (delta > 0) {
				const instantFps = 1000 / delta;
				// Use exponential moving average to smooth FPS (80% old, 20% new)
				const currentFps = useAppStore.getState().fps;
				const smoothedFps =
					currentFps === null
						? instantFps
						: 0.8 * currentFps + 0.2 * instantFps;
				setFps(smoothedFps);
			}
		}

		setLastFrameChangeTime(now);
	}, [playing, currentFrame, setFps, setLastFrameChangeTime]); // Only depend on playing and currentFrame changes

	useEffect(() => {
		const handleKeyDown = (event: KeyboardEvent) => {
			if (isInputFocused(event.target)) return;

			// Handle Shift + Arrow keys for bookmark navigation
			if (event.shiftKey && event.key === "ArrowLeft") {
				event.preventDefault();
				handlePreviousBookmark();
				return;
			}

			if (event.shiftKey && event.key === "ArrowRight") {
				event.preventDefault();
				handleNextBookmark();
				return;
			}

			switch (event.key) {
				case "ArrowLeft":
					event.preventDefault();
					handlePreviousFrame();
					break;

				case "ArrowRight":
					event.preventDefault();
					handleNextFrame();
					break;

				case " ":
					event.preventDefault();
					handleTogglePlayback();
					break;

				case "b":
				case "B":
					event.preventDefault();
					handleAddBookmark();
					break;

				case "Escape":
					event.preventDefault();
					updateFrameSelection(null);
					break;

				default:
					break;
			}
		};

		window.addEventListener("keydown", handleKeyDown);

		return () => {
			window.removeEventListener("keydown", handleKeyDown);
		};
	}, [
		handlePreviousFrame,
		handleNextFrame,
		handleTogglePlayback,
		handleAddBookmark,
		handlePreviousBookmark,
		handleNextBookmark,
		updateFrameSelection,
	]);
};
