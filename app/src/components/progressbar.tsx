import React, { useState, useEffect } from "react";
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
import CloudDownloadIcon from '@mui/icons-material/CloudDownload'; // Icon for waiting for data

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

// Main container with improved alignment and padding
const ProgressBarContainer = styled(Box)(() => ({
    width: "100%",
    padding: "12px 16px",
    backgroundColor: "#f4f6f8",
    borderRadius: "8px",
    boxShadow: "0px 3px 6px rgba(0,0,0,0.16)",
    display: "flex",
    alignItems: "center", // Vertically aligns all direct children
    gap: "16px",
    flexWrap: 'wrap',
}));

// Simplified wrapper for the slider. No longer needs complex padding or height.
const SliderControlsWrapper = styled(Box)(() => ({
    position: "relative", // Remains relative for bookmarks and highlights
    flexGrow: 1,
    display: "flex",
    alignItems: "center",
    minWidth: 250,
}));

// Repositioned bookmark to be closer to the slider track
const BookmarkIndicator = styled("div")<{ left: string }>(({ left }) => ({
    position: "absolute",
    left: left,
    top: "50%", // Align with the vertical center of the wrapper
    zIndex: 2,
    // Shift horizontally to center on the mark, and vertically to sit just above the track
    transform: "translate(-50%, -110%)",
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

// Highlight style adjusted for the new, slimmer progress bar
const SelectedFrameHighlight = styled("div")<{ left: string; width: string }>(
    ({ left, width }) => ({
        position: "absolute",
        left: left,
        width: width,
        height: "6px", // Slightly thicker than the track for visibility
        backgroundColor: "#42a5f5",
        borderRadius: "3px",
        zIndex: 1,
        pointerEvents: "none",
        top: "50%", // Center vertically on the slider's centerline
        transform: "translateY(-50%)",
        opacity: 0.8,
    })
);

// Animation for the "waiting for data" icon
const waitingAnimation = {
    '@keyframes waiting': {
        '0%, 100%': { transform: 'translateY(0)' },
        '50%': { transform: 'translateY(-3px)' },
    },
    animation: 'waiting 1.5s ease-in-out infinite',
};


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

    const [debouncedConnected, setDebouncedConnected] = useState(connected);
    const [debouncedIsFrameRendering, setDebouncedIsFrameRendering] = useState(isFrameRendering);

    useEffect(() => {
        setInputValue(String(step));
    }, [step]);

    useEffect(() => {
        setFpsInputValue(String(frameRate));
    }, [frameRate]);

    useEffect(() => {
        const handler = setTimeout(() => setDebouncedConnected(connected), connected ? 0 : 100);
        return () => clearTimeout(handler);
    }, [connected]);

    useEffect(() => {
        const handler = setTimeout(() => setDebouncedIsFrameRendering(isFrameRendering), isFrameRendering ? 100 : 0);
        return () => clearTimeout(handler);
    }, [isFrameRendering]);

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
        if (!selectedFrames || !selectedFrames.active || !selectedFrames.indices) return null;

        const framesToHighlight = [...selectedFrames.indices].sort((a, b) => a - b);
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
            return <SelectedFrameHighlight key={`selected-block-${index}`} left={leftPosition} width={width} />;
        });
    };

    return (
        <ProgressBarContainer>
            <TextField
                variant="outlined"
                size="small"
                value={inputValue}
                onChange={handleInputChange}
                onBlur={handleInputBlur}
                type="number"
                inputProps={{ min: 0, max: length, step: 1 }}
                InputProps={{
                    endAdornment: <InputAdornment position="end">/ {length}</InputAdornment>,
                }}
                sx={{ width: 120, flexShrink: 0 }}
            />

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
                        // Removed absolute positioning for natural flex alignment
                        "& .MuiSlider-track": { height: 4, borderRadius: 2 },
                        "& .MuiSlider-rail": { height: 4, borderRadius: 2 },
                        "& .MuiSlider-thumb": {
                            width: 16,
                            height: 16,
                            backgroundColor: '#fff',
                            border: '2px solid currentColor',
                            '&:hover, &.Mui-focusVisible': {
                                boxShadow: '0 0 0 8px rgba(25, 118, 210, 0.16)',
                            },
                        },
                    }}
                />
            </SliderControlsWrapper>

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

            <Tooltip title={selectedFrames.active ? "Hide Selections" : "Show Selections"}>
                <IconButton onClick={handleToggleSelectedFrames} color="primary">
                    {selectedFrames.active ? <VisibilityIcon /> : <VisibilityOffIcon />}
                </IconButton>
            </Tooltip>
            
            {/* Unified status indicator logic */}
            <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', minWidth: 40, flexShrink: 0 }}>
                <Tooltip
                    title={
                        !debouncedConnected
                            ? "Connecting..."
                            : debouncedIsFrameRendering
                            ? "Waiting for frame data..."
                            : "Connected to Server"
                    }
                >
                    {/* Wrapper div prevents Tooltip warnings on conditionally rendered/disabled items */}
                    <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: 24, width: 24 }}>
                        {!debouncedConnected ? (
                            <CircularProgress size={24} />
                        ) : debouncedIsFrameRendering ? (
                            <CloudDownloadIcon sx={{ color: "#1976d2", ...waitingAnimation }} />
                        ) : (
                            <WifiIcon sx={{ color: "#4caf50" }} />
                        )}
                    </Box>
                </Tooltip>
            </Box>
        </ProgressBarContainer>
    );
};

export default FrameProgressBar;
