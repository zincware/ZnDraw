import { useEffect, useCallback } from "react";
import { useAppStore } from "../store";
import { useStepControl } from "./useStepControl";

export const useKeyboardShortcuts = () => {
  const {
    currentFrame,
    frameCount,
    frame_selection,
    frameSelectionEnabled,
    playing,
    setPlaying,
    skipFrames,
    synchronizedMode,
    getIsFetching,
    addBookmark,
    bookmarks,
    setFps,
    setLastFrameChangeTime,
  } = useAppStore();

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

  const handlePreviousFrame = useCallback(async () => {
    // Pause playback when manually navigating
    if (playing) {
      setPlaying(false);
    }

    const navigableFrames = getNavigableFrames();
    const currentIndex = navigableFrames.indexOf(currentFrame);

    if (currentIndex > 0) {
      setStep(navigableFrames[currentIndex - 1]);
    } else {
      // Wrap to last frame
      setStep(navigableFrames[navigableFrames.length - 1]);
    }
  }, [currentFrame, getNavigableFrames, setStep, playing, setPlaying]);

  const handleNextFrame = useCallback(async () => {
    // Pause playback when manually navigating
    if (playing) {
      setPlaying(false);
    }

    const navigableFrames = getNavigableFrames();
    const currentIndex = navigableFrames.indexOf(currentFrame);

    if (currentIndex < navigableFrames.length - 1) {
      setStep(navigableFrames[currentIndex + 1]);
    } else {
      // Wrap to first frame
      setStep(navigableFrames[0]);
    }
  }, [
    currentFrame,
    getNavigableFrames,
    setStep,
    playing,
    setPlaying,
  ]);

  const handleTogglePlayback = useCallback(async () => {
    const navigableFrames = getNavigableFrames();
    const isAtLastFrame =
      currentFrame === navigableFrames[navigableFrames.length - 1];

    if (!playing) {
      // Starting playback
      if (isAtLastFrame) {
        // If at last frame, jump to first before playing
        setStep(navigableFrames[0]);
      }
      setPlaying(true);
    } else {
      // Stopping playback
      setPlaying(false);
    }
  }, [
    playing,
    currentFrame,
    getNavigableFrames,
    setStep,
    setPlaying,
  ]);

  // Add bookmark (B)
  const handleAddBookmark = useCallback(() => {
    addBookmark(currentFrame);
  }, [addBookmark, currentFrame]);

  // Jump to previous bookmark (Shift + ArrowLeft)
  const handlePreviousBookmark = useCallback(async () => {
    if (!bookmarks) return;

    // Pause playback when manually navigating
    if (playing) {
      setPlaying(false);
    }

    const bookmarkFrames = Object.keys(bookmarks)
      .map(Number)
      .sort((a, b) => a - b);

    if (bookmarkFrames.length === 0) return;

    // Find the nearest bookmark before current frame
    const prevBookmark = bookmarkFrames
      .reverse()
      .find((frame) => frame < currentFrame);

    if (prevBookmark !== undefined) {
      setStep(prevBookmark);
    } else {
      // Wrap to last bookmark
      setStep(bookmarkFrames[0]);
    }
  }, [bookmarks, currentFrame, setStep, playing, setPlaying]);

  // Jump to next bookmark (Shift + ArrowRight)
  const handleNextBookmark = useCallback(async () => {
    if (!bookmarks) return;

    // Pause playback when manually navigating
    if (playing) {
      setPlaying(false);
    }

    const bookmarkFrames = Object.keys(bookmarks)
      .map(Number)
      .sort((a, b) => a - b);

    if (bookmarkFrames.length === 0) return;

    // Find the nearest bookmark after current frame
    const nextBookmark = bookmarkFrames.find((frame) => frame > currentFrame);

    if (nextBookmark !== undefined) {
      setStep(nextBookmark);
    } else {
      // Wrap to first bookmark
      setStep(bookmarkFrames[0]);
    }
  }, [bookmarks, currentFrame, setStep, playing, setPlaying]);

  // Handle automatic frame advancement during playback
  useEffect(() => {
    if (!playing) return;

    const intervalId = setInterval(() => {
      // In synchronized mode, wait for all active geometries to finish fetching
      if (synchronizedMode && getIsFetching()) {
        return; // Skip frame advancement until geometries are ready
      }

      const navigableFrames = getNavigableFrames();
      const currentIndex = navigableFrames.indexOf(currentFrame);
      const lastIndex = navigableFrames.length - 1;

      // Advance by skipFrames amount
      const nextIndex = currentIndex + skipFrames;

      if (nextIndex < lastIndex) {
        setStep(navigableFrames[nextIndex]);
      } else if (currentIndex < lastIndex) {
        // If next skip would overshoot, go to last frame
        setStep(navigableFrames[lastIndex]);
      } else {
        // Already at last frame, stop playing
        setPlaying(false);
      }
    }, 33); // approximately 30 fps

    return () => clearInterval(intervalId);
  }, [
    playing,
    currentFrame,
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
        const smoothedFps = currentFps === null ? instantFps : 0.8 * currentFps + 0.2 * instantFps;
        setFps(smoothedFps);
      }
    }

    setLastFrameChangeTime(now);
  }, [playing, currentFrame, setFps, setLastFrameChangeTime]); // Only depend on playing and currentFrame changes


  useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      // Check if focus is on an input element
      const target = event.target as HTMLElement;
      const isInputFocused =
        target.tagName === "INPUT" ||
        target.tagName === "TEXTAREA" ||
        target.isContentEditable;

      // Don't trigger shortcuts when typing in forms
      if (isInputFocused) {
        return;
      }

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
  ]);
};
