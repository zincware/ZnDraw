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
import { socket } from "../socket";
import { useStepControl } from "../hooks/useStepControl";
import SelectionLayer from "./SelectionLayer";
import BookmarkLayer from "./BookmarkLayer";
import SelectionTrackOverlay from "./SelectionTrackOverlay";
import { LAYOUT_CONSTANTS } from "../constants/layout";

const FrameProgressBar = () => {
  const [isEditing, setIsEditing] = useState(false);
  const [inputValue, setInputValue] = useState("0");
  const sliderContainerRef = useRef<HTMLDivElement>(null);
  const [sliderWidth, setSliderWidth] = useState(0);

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
    addBookmark,
    deleteBookmark,
    frameSelectionEnabled,
    setFrameSelectionEnabled,
    synchronizedMode,
    setSynchronizedMode,
    getIsFetching,
  } = useAppStore();

  const { setStep, remoteLocked } = useStepControl();

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

  // This is the primary handler for the slider
  const handleSliderChange = (_e: Event, newFrame: number | number[]) => {
    // Apply snap-to-selected-frame logic if filtering is enabled
    const snappedFrame = findNearestSelectedFrame(newFrame as number);
    setStep(snappedFrame);
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


  // Auto-disable frame selection filter when frame selection is reset
  useEffect(() => {
    if (frameSelectionEnabled && (!frame_selection || frame_selection.length === 0)) {
      setFrameSelectionEnabled(false);
    }
  }, [frame_selection, frameSelectionEnabled, setFrameSelectionEnabled]);

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
      setStep(newFrame);
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
        height: `${LAYOUT_CONSTANTS.PROGRESSBAR_HEIGHT}px`,
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
          onBookmarkClick={setStep}
          onBookmarkDelete={deleteBookmark}
          onBookmarkEdit={addBookmark}
        />

        {/* Selection Track Overlay - Always shown to replace MUI track/rail */}
        <SelectionTrackOverlay
          selectedFrames={frame_selection}
          frameCount={frameCount}
          containerWidth={sliderWidth}
          enabled={frameSelectionEnabled}
          currentFrame={currentFrame}
          disabled={remoteLocked}
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
          disabled={remoteLocked}
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
