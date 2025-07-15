import React, { useState, useEffect, useRef } from "react";
import {
    Slider,
    TextField,
    InputAdornment,
    IconButton,
    Box,
    Typography,
    Tooltip,
    CircularProgress,
} from "@mui/material";
import { styled } from "@mui/system";
import { FaRegBookmark } from "react-icons/fa";
import WifiIcon from '@mui/icons-material/Wifi';
import VisibilityIcon from '@mui/icons-material/Visibility';
import VisibilityOffIcon from '@mui/icons-material/VisibilityOff';

interface IndicesState {
    active: boolean;
    indices: number[];
}

interface FrameProgressBarProps {
    length: number;
    step: number;
    setStep: (step: number) => void;
    selectedFrames: IndicesState;
    setSelectedFrames: (selectedFrames: IndicesState) => void;
    bookmarks: { [key: number]: string };
    setBookmarks: (bookmarks: { [key: number]: string }) => void;
    connected: boolean;
    frameRate: number;
    setFrameRate: (rate: number) => void;
    isFrameRendering: boolean;
}

const ProgressBarContainer = styled(Box)(() => ({
    width: "100%",
    padding: "12px 24px",
    backgroundColor: "#f4f6f8",
    borderRadius: "8px",
    boxShadow: "0px 3px 6px rgba(0,0,0,0.16)",
    display: "flex",
    alignItems: "center",
    gap: "24px",
    position: "relative",
    flexWrap: 'wrap',
}));

const SliderControlsWrapper = styled(Box)(() => ({
    position: "relative",
    flexGrow: 1,
    height: 48,
    display: "flex",
    alignItems: "center",
    paddingTop: 18,
    paddingBottom: 8,
    minWidth: 250,
}));

const BookmarkIndicator = styled("div")<{ left: string }>(({ left }) => ({
    position: "absolute",
    left: left,
    top: 0,
    zIndex: 2,
    transform: 'translateX(-50%)',
    "& svg": {
        color: "#ff9800",
        fontSize: "1.1em",
        cursor: "pointer",
        transition: "transform 0.2s ease-in-out",
        "&:hover": {
            transform: "scale(1.3)",
        },
    },
}));

const SelectedFrameHighlight = styled("div")<{ left: string; width: string }>(
    ({ left, width }) => ({
        position: "absolute",
        left: left,
        width: width,
        height: "10px",
        backgroundColor: "#42a5f5",
        borderRadius: "3px",
        zIndex: 1,
        pointerEvents: "none",
        bottom: "10px", // Adjusted to align with the bottom of the progress bar track
        opacity: 0.8,
    })
);

export const FrameProgressBar: React.FC<FrameProgressBarProps> = ({
    length,
    step,
    setStep,
    selectedFrames,
    setSelectedFrames,
    bookmarks,
    setBookmarks,
    connected,
    frameRate,
    setFrameRate,
    isFrameRendering,
}) => {
    const [inputValue, setInputValue] = useState<string>(String(step));
    const [fpsInputValue, setFpsInputValue] = useState<string>(String(frameRate));

    useEffect(() => {
        setInputValue(String(step));
    }, [step]);

    useEffect(() => {
        setFpsInputValue(String(frameRate));
    }, [frameRate]);

    const handleSliderChange = (event: Event, newValue: number | number[]) => {
        setStep(newValue as number);
    };

    const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        const value = event.target.value;
        setInputValue(value);
        const numValue = Number(value);
        if (!isNaN(numValue) && numValue >= 0 && numValue <= length) {
            setStep(numValue);
        }
    };

    const handleInputBlur = () => {
        const numValue = Number(inputValue);
        if (isNaN(numValue) || numValue < 0 || numValue > length) {
            setInputValue(String(step));
        }
    };

    const handleBookmarkClick = (event: React.MouseEvent, frameNumber: number) => {
        if (event.shiftKey) {
            const newBookmarks = { ...bookmarks };
            delete newBookmarks[frameNumber];
            setBookmarks(newBookmarks);
        } else {
            setStep(frameNumber);
        }
    };

    const handleToggleSelectedFrames = () => {
        setSelectedFrames({ ...selectedFrames, active: !selectedFrames.active });
    };

    const handleFpsInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        const value = event.target.value;
        setFpsInputValue(value);
        const numValue = Number(value);
        if (!isNaN(numValue) && numValue > 0) {
            setFrameRate(numValue);
        }
    };

    const handleFpsInputBlur = () => {
        const numValue = Number(fpsInputValue);
        if (isNaN(numValue) || numValue <= 0) {
            setFpsInputValue(String(frameRate));
        }
    };

    const renderBookmarks = () => {
        return Object.keys(bookmarks).map((key) => {
            const position = Number.parseInt(key);
            const leftPosition = `${(position / length) * 100}%`;
            return (
                <Tooltip key={`bookmark-${position}`} title={bookmarks[position]}>
                    <BookmarkIndicator left={leftPosition}>
                        <FaRegBookmark
                            onClick={(e) => handleBookmarkClick(e, position)}
                        />
                    </BookmarkIndicator>
                </Tooltip>
            );
        });
    };

    const renderSelectedFrames = () => {
        if (!selectedFrames || !selectedFrames.active || !selectedFrames.indices) {
            return null;
        }

        const framesToHighlight = [...selectedFrames.indices];
        framesToHighlight.sort((a, b) => a - b);

        const highlightBlocks: { start: number; end: number }[] = [];
        if (framesToHighlight.length > 0) {
            let currentBlock = { start: framesToHighlight[0], end: framesToHighlight[0] };
            for (let i = 1; i < framesToHighlight.length; i++) {
                if (framesToHighlight[i] === currentBlock.end + 1) {
                    currentBlock.end = framesToHighlight[i];
                } else {
                    highlightBlocks.push(currentBlock);
                    currentBlock = { start: framesToHighlight[i], end: framesToHighlight[i] };
                }
            }
            highlightBlocks.push(currentBlock);
        }

        return highlightBlocks.map((block, index) => {
            const leftPosition = `${(block.start / length) * 100}%`;
            const width = `${((block.end - block.start + 1) / length) * 100}%`;

            return (
                <SelectedFrameHighlight
                    key={`selected-block-${index}`}
                    left={leftPosition}
                    width={width}
                />
            );
        });
    };

    return (
        <ProgressBarContainer>
            {/* Frame Number Input */}
            <TextField
                variant="outlined"
                size="small"
                value={inputValue}
                onChange={handleInputChange}
                onBlur={handleInputBlur}
                type="number"
                inputProps={{
                    min: 0,
                    max: length,
                    step: 1,
                }}
                InputProps={{
                    endAdornment: (
                        <InputAdornment position="end">
                            <Typography variant="body2" sx={{ color: "#616161" }}>/ {length}</Typography>
                        </InputAdornment>
                    ),
                }}
                sx={{ width: 120, flexShrink: 0 }}
            />

            {/* Slider with Highlights and Bookmarks */}
            <SliderControlsWrapper>
                {renderSelectedFrames()}
                {renderBookmarks()}
                <Slider
                    value={step}
                    min={0}
                    max={length}
                    step={1}
                    onChange={handleSliderChange}
                    aria-labelledby="frame-slider"
                    sx={{
                        "& .MuiSlider-track": {
                            height: 8,
                            borderRadius: 4,
                            backgroundColor: selectedFrames.active ? 'transparent' : '#1976d2',
                        },
                        "& .MuiSlider-rail": {
                            height: 8,
                            borderRadius: 4,
                            backgroundColor: selectedFrames.active ? 'rgba(128, 128, 128, 0.3)' : '#e0e0e0',
                        },
                        "& .MuiSlider-thumb": {
                            width: 18,
                            height: 18,
                            // marginTop: '-50px',
                            backgroundColor: '#1976d2',
                            border: '2px solid white',
                            boxShadow: '0px 0px 5px rgba(0,0,0,0.2)',
                        },
                        width: '100%',
                        position: 'absolute',
                        bottom: 0,
                    }}
                />
            </SliderControlsWrapper>

            {/* FPS Custom Input */}
            <TextField
                variant="outlined"
                size="small"
                label="FPS"
                value={fpsInputValue}
                onChange={handleFpsInputChange}
                onBlur={handleFpsInputBlur}
                type="number"
                inputProps={{ min: 1 }}
                sx={{ width: 100, flexShrink: 0 }}
            />

            {/* Toggle for Selected Frames (Clickable Icon) */}
            <Tooltip title={selectedFrames.active ? "Hide Selections" : "Show Selections"}>
                <IconButton onClick={handleToggleSelectedFrames} sx={{ color: "#1976d2" }}>
                    {selectedFrames.active ? <VisibilityIcon /> : <VisibilityOffIcon />}
                </IconButton>
            </Tooltip>

            {/* Connection Status and Frame Loading Indicators */}
            <Box sx={{ display: 'flex', gap: 1.5, minWidth: 60, justifyContent: 'flex-end', flexShrink: 0 }}>
                {/* Connection Status Indicator */}
                <Tooltip title={connected ? "Connected to Server" : "Connecting..."}>
                    {connected ? (
                        <WifiIcon sx={{ color: "#4caf50", fontSize: 22 }} />
                    ) : (
                        <CircularProgress size={22} />
                    )}
                </Tooltip>

                {/* Frame Loading Indicator */}
                {isFrameRendering && (
                    <Tooltip title="Frame is loading...">
                        <CircularProgress size={22} />
                    </Tooltip>
                )}
            </Box>
        </ProgressBarContainer>
    );
};

export default FrameProgressBar;
