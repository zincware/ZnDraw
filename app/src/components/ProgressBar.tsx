import { Box, Slider, TextField, Typography, CircularProgress } from "@mui/material";
import { useState } from "react";
import WifiIcon from '@mui/icons-material/Wifi';
import CloudDownloadIcon from '@mui/icons-material/CloudDownload';


const FrameProgressBar = () => {
    const [currentFrame, setCurrentFrame] = useState(0);
    const [totalFrames, setTotalFrames] = useState(100);
    const [isEditing, setIsEditing] = useState(false);
    const [inputValue, setInputValue] = useState(currentFrame.toString());
    const [isLoading, setIsLoading] = useState(false);
    const [isConnected, setIsConnected] = useState(true);
    const [skipFrames, setSkipFrames] = useState(1);

    const waitingAnimation = {
        '@keyframes pulse': {
            '0%': { opacity: 1 },
            '50%': { opacity: 0.5 },
            '100%': { opacity: 1 },
        },
        animation: 'pulse 1.5s ease-in-out infinite',
    };

    const compactTextFieldSx = {
        '& .MuiInputBase-input': {
            textAlign: 'center',
            fontSize: '0.875rem',
            padding: '4px 8px',
        },
        '& .MuiInputLabel-root': {
            fontSize: '0.75rem',
        }
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
        if (!isNaN(newFrame) && newFrame >= 0 && newFrame <= totalFrames) {
            setCurrentFrame(newFrame);
        } else {
            setInputValue(currentFrame.toString());
        }
        setIsEditing(false);
    };

    const handleInputKeyPress = (event: React.KeyboardEvent) => {
        if (event.key === 'Enter') {
            handleInputSubmit();
        } else if (event.key === 'Escape') {
            setInputValue(currentFrame.toString());
            setIsEditing(false);
        }
    };

    const handleInputBlur = () => {
        handleInputSubmit();
    };

    const handleSkipFramesChange = (event: React.ChangeEvent<HTMLInputElement>) => {
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
            return (
                <WifiIcon sx={{ color: "#4caf50", fontSize: 22 }} />
            );
        } else {
            return (
                <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                    <CircularProgress size={20} />
                </Box>
            );
        }
    };

    return (
        <Box sx={{
            position: 'fixed',
            bottom: 0,
            left: 0,
            right: 0,
            padding: '5px 10px',
            backgroundColor: 'background.paper',
            border: '1px solid',
            borderColor: 'divider',
            boxShadow: '0 -2px 8px rgba(0,0,0,0.1)',
            zIndex: (theme) => theme.zIndex.drawer + 1,
            display: 'flex',
            alignItems: 'center',
            gap: 4,
        }}>
            <Box sx={{ minWidth: '100px', cursor: 'pointer' }} onClick={handleDisplayClick}>
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
                            width: '80px',
                            '& .MuiInputBase-input': {
                                textAlign: 'center',
                                fontSize: '0.875rem',
                            }
                        }}
                    />
                ) : (
                    <Typography
                        variant="body2"
                        sx={{
                            textAlign: 'center',
                            padding: '8px',
                            borderRadius: '4px',
                            '&:hover': {
                                backgroundColor: 'action.hover',
                            }
                        }}
                    >
                        {currentFrame} / {totalFrames}
                    </Typography>
                )}
            </Box>
            <Slider
                orientation="horizontal"
                value={currentFrame}
                onChange={(_e, val) => setCurrentFrame(val as number)}
                max={totalFrames}
                aria-label="Frame Progress"
                valueLabelDisplay="auto"
                sx={{
                    flexGrow: 1,
                    '& .MuiSlider-thumb': {
                        transition: 'none !important',
                    },
                    '& .MuiSlider-track': {
                        transition: 'none !important',
                    },
                    '& .MuiSlider-rail': {
                        transition: 'none !important',
                    }
                }}
            />
            <Box sx={{
                display: 'flex',
                alignItems: 'center',
                gap: 1,
                minWidth: '120px'
            }}>
                <TextField
                    variant="outlined"
                    size="small"
                    label="Skip"
                    value={skipFrames}
                    onChange={handleSkipFramesChange}
                    type="number"
                    slotProps={{
                        htmlInput: { min: 1 }
                    }}
                    sx={{ ...compactTextFieldSx, width: 80 }}
                />
                {renderConnectionStatus()}
            </Box>
        </Box>
    );
};

export default FrameProgressBar;