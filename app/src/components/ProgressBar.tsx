import {
  Box,
  Slider,
  TextField,
  Typography,
  CircularProgress,
  IconButton,
  Tooltip,
} from "@mui/material";
import { useState, useMemo, useRef, useEffect } from "react";
import WifiIcon from "@mui/icons-material/Wifi";
import CloudDownloadIcon from "@mui/icons-material/CloudDownload";
import FilterListIcon from "@mui/icons-material/FilterList";
import SyncIcon from "@mui/icons-material/Sync";
import SyncDisabledIcon from "@mui/icons-material/SyncDisabled";
import { useAppStore } from "../store";
import { throttle } from "lodash";
import { socket } from "../socket";
import { usePresenterMode } from "../hooks/usePresenterMode";
import { useAtomicFrameSet } from "../hooks/useAtomicFrameSet";
import SelectionLayer from "./SelectionLayer";
import BookmarkLayer from "./BookmarkLayer";
import SelectionTrackOverlay from "./SelectionTrackOverlay";

const FrameProgressBar = () => {
  const [isEditing, setIsEditing] = useState(false);
  const [inputValue, setInputValue] = useState("0");
  const sliderContainerRef = useRef<HTMLDivElement>(null);
  const [sliderWidth, setSliderWidth] = useState(0);

  // Refs to manage the scrubbing logic without causing re-renders
  const scrubTimerRef = useRef<NodeJS.Timeout | null>(null);
  const isScrubbingRef = useRef(false);

  const {
    currentFrame,
    setCurrentFrame,
    frameCount,
    isConnected,
    isLoading,
    skipFrames,
    setSkipFrames,
    frame_selection,
    bookmarks,
    deleteBookmark,
    frameSelectionEnabled,
    setFrameSelectionEnabled,
    synchronizedMode,
    setSynchronizedMode,
    getIsFetching,
  } = useAppStore();

  const { requestPresenterMode, releasePresenterMode, presenterMode } =
    usePresenterMode();
  const setFrameAtomic = useAtomicFrameSet();

  // Find nearest selected frame when filtering is enabled
  const findNearestSelectedFrame = (targetFrame: number): number => {
    if (
      !frameSelectionEnabled ||
      !frame_selection ||
      frame_selection.length === 0
    ) {
      return targetFrame;
    }

    let nearest = frame_selection[0];
    let minDistance = Math.abs(targetFrame - nearest);

    for (const frame of frame_selection) {
      const distance = Math.abs(targetFrame - frame);
      if (distance < minDistance) {
        minDistance = distance;
        nearest = frame;
      }
    }

    return nearest;
  };

  // Emits CONTINUOUS updates while dragging
  const throttledFrameUpdate = useMemo(
    () =>
      throttle((frame: number) => {
        // Only send if we are confirmed to be the presenter
        if (presenterMode === "presenting") {
          socket.emit("set_frame_continuous", { frame });
        }
      }, 100),
    [presenterMode],
  );

  // This is the primary handler for the slider
  const handleSliderChange = (_e: Event, newFrame: number | number[]) => {
    // Apply snap-to-selected-frame logic if filtering is enabled
    const snappedFrame = findNearestSelectedFrame(newFrame as number);
    setCurrentFrame(snappedFrame);

    // If we're already in scrubbing mode, just send throttled updates
    if (isScrubbingRef.current) {
      throttledFrameUpdate(snappedFrame);
      return;
    }

    // If a timer is running, it means this is the 2nd change event -> it's a scrub!
    if (scrubTimerRef.current) {
      clearTimeout(scrubTimerRef.current);
      scrubTimerRef.current = null;
      isScrubbingRef.current = true;
      requestPresenterMode(); // Request presenter mode
    }
    // This is the first change event. Treat it as an atomic jump for now.
    else {
      socket.emit(
        "set_frame_atomic",
        { frame: snappedFrame },
        (response: any) => {
          if (response && !response.success) {
            console.error(
              `Failed to set frame: ${response.error} - ${response.message}`,
            );
          }
        },
      );
      // Start a timer. If it completes, the interaction was just a click.
      scrubTimerRef.current = setTimeout(() => {
        scrubTimerRef.current = null;
      }, 150); // A short delay to detect a drag
    }
  };

  // This handler cleans up after the interaction ends
  const handleScrubEnd = () => {
    // Always clear any pending timer
    if (scrubTimerRef.current) {
      clearTimeout(scrubTimerRef.current);
      scrubTimerRef.current = null;
    }

    // If we were in scrubbing mode, release the token
    if (isScrubbingRef.current) {
      throttledFrameUpdate.flush(); // Send the final frame
      releasePresenterMode(); // Release presenter mode
    }

    // Reset the scrubbing flag for the next interaction
    isScrubbingRef.current = false;
  };

  // Measure slider container width
  useEffect(() => {
    const updateWidth = () => {
      if (sliderContainerRef.current) {
        setSliderWidth(sliderContainerRef.current.offsetWidth);
      }
    };

    updateWidth();
    window.addEventListener("resize", updateWidth);
    return () => window.removeEventListener("resize", updateWidth);
  }, []);

  // Ensure timers are cleaned up if the component unmounts
  useEffect(() => {
    return () => {
      if (scrubTimerRef.current) {
        clearTimeout(scrubTimerRef.current);
      }
      // Release presenter mode if still active
      if (isScrubbingRef.current) {
        releasePresenterMode();
      }
    };
  }, [releasePresenterMode]);

  const waitingAnimation = {
    "@keyframes pulse": {
      "0%": { opacity: 1 },
      "50%": { opacity: 0.5 },
      "100%": { opacity: 1 },
    },
    animation: "pulse 1.5s ease-in-out infinite",
  };

  const compactTextFieldSx = {
    "& .MuiInputBase-input": {
      textAlign: "center",
      fontSize: "0.875rem",
      padding: "4px 8px",
    },
    "& .MuiInputLabel-root": {
      fontSize: "0.75rem",
    },
  };

  const handleDisplayClick = () => {
    setIsEditing(true);
    setInputValue(currentFrame.toString());
  };

  const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setInputValue(event.target.value);
  };

  const handleInputSubmit = () => {
    const newFrame = parseInt(inputValue, 10);
    if (!isNaN(newFrame) && newFrame >= 0 && newFrame <= frameCount - 1) {
      setFrameAtomic(newFrame);
    } else {
      setInputValue(currentFrame.toString());
    }
    setIsEditing(false);
  };

  const handleInputKeyPress = (event: React.KeyboardEvent) => {
    if (event.key === "Enter") {
      handleInputSubmit();
    } else if (event.key === "Escape") {
      setInputValue(currentFrame.toString());
      setIsEditing(false);
    }
  };

  const handleInputBlur = () => {
    handleInputSubmit();
  };

  const handleSkipFramesChange = (
    event: React.ChangeEvent<HTMLInputElement>,
  ) => {
    const value = parseInt(event.target.value, 10);
    if (!isNaN(value) && value > 0) {
      setSkipFrames(value);
    }
  };

  const renderConnectionStatus = () => {
    if (isLoading) {
      return (
        <CloudDownloadIcon
          sx={{ color: "#1976d2", fontSize: 22, ...waitingAnimation }}
        />
      );
    } else if (isConnected) {
      return <WifiIcon sx={{ color: "#4caf50", fontSize: 22 }} />;
    } else {
      return (
        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          <CircularProgress size={20} />
        </Box>
      );
    }
  };

  return (
    <Box
      sx={{
        padding: "5px 10px",
        backgroundColor: "background.paper",
        borderTop: "1px solid",
        borderColor: "divider",
        boxShadow: "0 -2px 8px rgba(0,0,0,0.1)",
        display: "flex",
        alignItems: "center",
        gap: 4,
        flexShrink: 0,
        zIndex: (theme) => theme.zIndex.drawer + 1,
      }}
    >
      <Box
        sx={{ minWidth: "100px", cursor: "pointer" }}
        onClick={handleDisplayClick}
      >
        {isEditing ? (
          <TextField
            value={inputValue}
            onChange={handleInputChange}
            onKeyDown={handleInputKeyPress}
            onBlur={handleInputBlur}
            size="small"
            variant="outlined"
            autoFocus
            sx={{
              width: "80px",
              "& .MuiInputBase-input": {
                textAlign: "center",
                fontSize: "0.875rem",
              },
            }}
          />
        ) : (
          <Typography
            variant="body2"
            sx={{
              textAlign: "center",
              padding: "8px",
              borderRadius: "4px",
              "&:hover": {
                backgroundColor: "action.hover",
              },
            }}
          >
            {currentFrame} / {frameCount - 1}
          </Typography>
        )}
      </Box>
      {/* Timeline annotations container with layers */}
      <Box ref={sliderContainerRef} sx={{ flexGrow: 1, position: "relative" }}>
        {/* Bookmark Layer - Top tier */}
        <BookmarkLayer
          bookmarks={bookmarks}
          frameCount={frameCount}
          currentFrame={currentFrame}
          containerWidth={sliderWidth}
          onBookmarkClick={setFrameAtomic}
          onBookmarkDelete={deleteBookmark}
        />

        {/* Selection Track Overlay - Always shown to replace MUI track/rail */}
        <SelectionTrackOverlay
          selectedFrames={frame_selection}
          frameCount={frameCount}
          containerWidth={sliderWidth}
          enabled={frameSelectionEnabled}
          currentFrame={currentFrame}
          disabled={presenterMode === "locked"}
        />

        {/* Selection Markers - Shows individual selected frames when filtering is disabled */}
        {!frameSelectionEnabled && (
          <SelectionLayer
            selectedFrames={frame_selection}
            frameCount={frameCount}
            currentFrame={currentFrame}
            containerWidth={sliderWidth}
          />
        )}

        {/* Slider - Bottom tier */}
        <Slider
          orientation="horizontal"
          value={currentFrame}
          onChange={handleSliderChange}
          onMouseUp={handleScrubEnd}
          onMouseDown={() => {
            isScrubbingRef.current = false;
          }}
          disabled={presenterMode === "locked"}
          max={frameCount - 1}
          aria-label="Frame Progress"
          valueLabelDisplay="auto"
          sx={{
            position: "relative",
            top: "2px", // Shift slider down to align with overlay
            "& .MuiSlider-thumb": {
              transition: "none !important",
              zIndex: 2,
            },
            "& .MuiSlider-track": {
              transition: "none !important",
              display: "none", // Always hide - use custom overlay
            },
            "& .MuiSlider-rail": {
              transition: "none !important",
              display: "none", // Always hide - use custom overlay
            },
          }}
        />
      </Box>
      <Box
        sx={{
          display: "flex",
          alignItems: "center",
          gap: 1,
          minWidth: "120px",
        }}
      >
        <TextField
          variant="outlined"
          size="small"
          label="Skip"
          value={skipFrames}
          onChange={handleSkipFramesChange}
          type="number"
          slotProps={{
            htmlInput: { min: 1 },
          }}
          sx={{ ...compactTextFieldSx, width: 80 }}
        />
        {frame_selection && frame_selection.length > 0 && (
          <Tooltip
            title={
              frameSelectionEnabled
                ? "Disable frame selection filter"
                : "Enable frame selection filter"
            }
          >
            <IconButton
              size="small"
              onClick={() => setFrameSelectionEnabled(!frameSelectionEnabled)}
              sx={{
                color: frameSelectionEnabled
                  ? "primary.main"
                  : "text.secondary",
                "&:hover": {
                  backgroundColor: "action.hover",
                },
              }}
            >
              <FilterListIcon fontSize="small" />
            </IconButton>
          </Tooltip>
        )}
        <Tooltip
          title={
            synchronizedMode
              ? "Synchronized mode enabled - playback waits for all active geometries to finish loading"
              : "Synchronized mode disabled - playback uses fixed FPS"
          }
        >
          <IconButton
            size="small"
            onClick={() => setSynchronizedMode(!synchronizedMode)}
            sx={{
              color: synchronizedMode
                ? "primary.main"
                : "text.secondary",
              ...(getIsFetching() && synchronizedMode ? waitingAnimation : {}),
              "&:hover": {
                backgroundColor: "action.hover",
              },
            }}
          >
            {synchronizedMode ? <SyncIcon fontSize="small" /> : <SyncDisabledIcon fontSize="small" />}
          </IconButton>
        </Tooltip>
        {renderConnectionStatus()}
      </Box>
    </Box>
  );
};

export default FrameProgressBar;
