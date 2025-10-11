import { useEffect, useCallback } from "react";
import { useAppStore } from "../store";
import { usePresenterToken } from "./usePresenterToken";
import { useAtomicFrameSet } from "./useAtomicFrameSet";

export const useKeyboardShortcuts = () => {
  const {
    currentFrame,
    setCurrentFrame,
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
  } = useAppStore();

  const goToFrameAtomic = useAtomicFrameSet();
  const { requestToken, releaseToken, setFrame } = usePresenterToken(goToFrameAtomic);

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

  // For playback - uses continuous mode with presenter token
  const goToFrameContinuous = useCallback(
    (frame: number) => {
      setCurrentFrame(frame);
      setFrame(frame);
    },
    [setCurrentFrame, setFrame],
  );

  const handlePreviousFrame = useCallback(() => {
    const navigableFrames = getNavigableFrames();
    const currentIndex = navigableFrames.indexOf(currentFrame);

    if (currentIndex > 0) {
      goToFrameAtomic(navigableFrames[currentIndex - 1]);
    } else {
      // Wrap to last frame
      goToFrameAtomic(navigableFrames[navigableFrames.length - 1]);
    }
  }, [currentFrame, getNavigableFrames, goToFrameAtomic]);

  const handleNextFrame = useCallback(() => {
    const navigableFrames = getNavigableFrames();
    const currentIndex = navigableFrames.indexOf(currentFrame);

    if (currentIndex < navigableFrames.length - 1) {
      goToFrameAtomic(navigableFrames[currentIndex + 1]);
    } else {
      // Wrap to first frame
      goToFrameAtomic(navigableFrames[0]);
      // Stop playback when reaching the last frame
      if (playing) {
        releaseToken();
        setPlaying(false);
      }
    }
  }, [
    currentFrame,
    getNavigableFrames,
    goToFrameAtomic,
    playing,
    setPlaying,
    releaseToken,
  ]);

  const handleTogglePlayback = useCallback(async () => {
    const navigableFrames = getNavigableFrames();
    const isAtLastFrame =
      currentFrame === navigableFrames[navigableFrames.length - 1];

    if (!playing) {
      // Starting playback - request presenter token
      const success = await requestToken();
      if (success) {
        // If at last frame, jump to first before playing
        if (isAtLastFrame) {
          goToFrameContinuous(navigableFrames[0]);
        }
        setPlaying(true);
      }
    } else {
      // Stopping playback - release presenter token
      releaseToken();
      setPlaying(false);
    }
  }, [
    playing,
    setPlaying,
    currentFrame,
    getNavigableFrames,
    requestToken,
    releaseToken,
    goToFrameContinuous,
  ]);

  // Add bookmark (B)
  const handleAddBookmark = useCallback(() => {
    addBookmark(currentFrame);
  }, [addBookmark, currentFrame]);

  // Jump to previous bookmark (Shift + ArrowLeft)
  const handlePreviousBookmark = useCallback(() => {
    if (!bookmarks) return;

    const bookmarkFrames = Object.keys(bookmarks)
      .map(Number)
      .sort((a, b) => a - b);

    if (bookmarkFrames.length === 0) return;

    // Find the nearest bookmark before current frame
    const prevBookmark = bookmarkFrames
      .reverse()
      .find((frame) => frame < currentFrame);

    if (prevBookmark !== undefined) {
      goToFrameAtomic(prevBookmark);
    } else {
      // Wrap to last bookmark
      goToFrameAtomic(bookmarkFrames[0]);
    }
  }, [bookmarks, currentFrame, goToFrameAtomic]);

  // Jump to next bookmark (Shift + ArrowRight)
  const handleNextBookmark = useCallback(() => {
    if (!bookmarks) return;

    const bookmarkFrames = Object.keys(bookmarks)
      .map(Number)
      .sort((a, b) => a - b);

    if (bookmarkFrames.length === 0) return;

    // Find the nearest bookmark after current frame
    const nextBookmark = bookmarkFrames.find((frame) => frame > currentFrame);

    if (nextBookmark !== undefined) {
      goToFrameAtomic(nextBookmark);
    } else {
      // Wrap to first bookmark
      goToFrameAtomic(bookmarkFrames[0]);
    }
  }, [bookmarks, currentFrame, goToFrameAtomic]);

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
        goToFrameContinuous(navigableFrames[nextIndex]);
      } else if (currentIndex < lastIndex) {
        // If next skip would overshoot, go to last frame
        goToFrameContinuous(navigableFrames[lastIndex]);
      } else {
        // Already at last frame, stop playing and release presenter token
        releaseToken();
        setPlaying(false);
      }
    }, 33); // approximately 30 fps

    return () => clearInterval(intervalId);
  }, [
    playing,
    currentFrame,
    getNavigableFrames,
    goToFrameContinuous,
    setPlaying,
    skipFrames,
    releaseToken,
    synchronizedMode,
    getIsFetching,
  ]);

  // Cleanup presenter token when component unmounts while playing
  useEffect(() => {
    return () => {
      if (playing) {
        releaseToken();
      }
    };
  }, []);

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
