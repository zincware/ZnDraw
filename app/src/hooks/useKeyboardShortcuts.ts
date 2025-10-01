import { useEffect, useCallback } from 'react';
import { useAppStore } from '../store';
import { socket } from '../socket';

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
  } = useAppStore();

  const getNavigableFrames = useCallback((): number[] => {
    if (frameSelectionEnabled && frame_selection && frame_selection.length > 0) {
      return [...frame_selection].sort((a, b) => a - b);
    }
    return Array.from({ length: frameCount }, (_, i) => i);
  }, [frameSelectionEnabled, frame_selection, frameCount]);

  const goToFrame = useCallback((frame: number) => {
    setCurrentFrame(frame);
    socket.emit('set_frame_atomic', { frame }, (response: any) => {
      if (response && !response.success) {
        console.error(`Failed to set frame: ${response.error}`);
      }
    });
  }, [setCurrentFrame]);

  const handlePreviousFrame = useCallback(() => {
    const navigableFrames = getNavigableFrames();
    const currentIndex = navigableFrames.indexOf(currentFrame);

    if (currentIndex > 0) {
      goToFrame(navigableFrames[currentIndex - 1]);
    } else {
      // Wrap to last frame
      goToFrame(navigableFrames[navigableFrames.length - 1]);
    }
  }, [currentFrame, getNavigableFrames, goToFrame]);

  const handleNextFrame = useCallback(() => {
    const navigableFrames = getNavigableFrames();
    const currentIndex = navigableFrames.indexOf(currentFrame);

    if (currentIndex < navigableFrames.length - 1) {
      goToFrame(navigableFrames[currentIndex + 1]);
    } else {
      // Wrap to first frame
      goToFrame(navigableFrames[0]);
      // Stop playback when reaching the last frame
      if (playing) {
        setPlaying(false);
      }
    }
  }, [currentFrame, getNavigableFrames, goToFrame, playing, setPlaying]);

  const handleTogglePlayback = useCallback(() => {
    const navigableFrames = getNavigableFrames();
    const isAtLastFrame = currentFrame === navigableFrames[navigableFrames.length - 1];

    // If at last frame and not playing, jump to first and start playing
    if (isAtLastFrame && !playing) {
      goToFrame(navigableFrames[0]);
      setPlaying(true);
    } else {
      // Otherwise just toggle playback
      setPlaying(!playing);
    }
  }, [playing, setPlaying, currentFrame, getNavigableFrames, goToFrame]);

  // Handle automatic frame advancement during playback
  useEffect(() => {
    if (!playing) return;

    const intervalId = setInterval(() => {
      const navigableFrames = getNavigableFrames();
      const currentIndex = navigableFrames.indexOf(currentFrame);
      const lastIndex = navigableFrames.length - 1;

      // Advance by skipFrames amount
      const nextIndex = currentIndex + skipFrames;

      if (nextIndex < lastIndex) {
        goToFrame(navigableFrames[nextIndex]);
      } else if (currentIndex < lastIndex) {
        // If next skip would overshoot, go to last frame
        goToFrame(navigableFrames[lastIndex]);
      } else {
        // Already at last frame, stop playing
        setPlaying(false);
      }
    }, 33); // approximately 30 fps

    return () => clearInterval(intervalId);
  }, [playing, currentFrame, getNavigableFrames, goToFrame, setPlaying, skipFrames]);

  useEffect(() => {
    const handleKeyDown = (event: KeyboardEvent) => {
      // Check if focus is on an input element
      const target = event.target as HTMLElement;
      const isInputFocused =
        target.tagName === 'INPUT' ||
        target.tagName === 'TEXTAREA' ||
        target.isContentEditable;

      // Don't trigger shortcuts when typing in forms
      if (isInputFocused) {
        return;
      }

      switch (event.key) {
        case 'ArrowLeft':
          event.preventDefault();
          handlePreviousFrame();
          break;

        case 'ArrowRight':
          event.preventDefault();
          handleNextFrame();
          break;

        case ' ':
          event.preventDefault();
          handleTogglePlayback();
          break;

        default:
          break;
      }
    };

    window.addEventListener('keydown', handleKeyDown);

    return () => {
      window.removeEventListener('keydown', handleKeyDown);
    };
  }, [handlePreviousFrame, handleNextFrame, handleTogglePlayback]);
};
