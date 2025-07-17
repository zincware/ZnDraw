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
import { styled, useTheme } from "@mui/material/styles";
import WifiIcon from "@mui/icons-material/Wifi";
import VisibilityIcon from "@mui/icons-material/Visibility";
import VisibilityOffIcon from "@mui/icons-material/VisibilityOff";
import CloudDownloadIcon from "@mui/icons-material/CloudDownload";
import { type IndicesState } from "./utils";

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

const SliderControlsWrapper = styled(Box)(() => ({
	position: "relative",
	flexGrow: 1,
	display: "flex",
	alignItems: "center",
	minWidth: 250,
	margin: "0 8px",
}));

// CORRECTED: Bookmark positioning logic
const BookmarkIndicator = styled("div")<{ left: string }>(({ left }) => ({
	position: "absolute",
	left: left,
	top: "50%", // Start from the vertical center of the parent
	// Horizontally center on the mark.
	// Vertically, shift it up by its own height (7px) plus half the rail height (2px)
	// to make the tip sit exactly on top of the slider rail.
	transform: "translate(-50%, -9px)",
	width: "12px",
	height: "7px", // The height of the triangle
	zIndex: 2,
	cursor: "pointer",
	"& svg": {
		display: "block",
		width: "100%",
		height: "100%",
		filter: "drop-shadow(0px 1px 1px rgba(0,0,0,0.2))",
	},
	"& svg path": {
		fill: "#ff9800",
		transition: "fill 0.2s ease-in-out",
	},
	"&:hover svg path": {
		fill: "#e65100", // Darken on hover
	},
}));

const compactTextFieldSx = {
	"& .MuiInputBase-root": { height: 30, fontSize: "0.875rem" },
	"& .MuiOutlinedInput-input": { padding: "4px 8px" },
	"& .MuiInputLabel-root": { top: "-6px" },
	"& .MuiInputLabel-shrink": { top: "0px" },
};

const waitingAnimation = {
	"@keyframes waiting": {
		"0%, 100%": { transform: "translateY(0)" },
		"50%": { transform: "translateY(-3px)" },
	},
	animation: "waiting 1.5s ease-in-out infinite",
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
	const theme = useTheme();
	const [inputValue, setInputValue] = useState<string>(String(step));
	const [fpsInputValue, setFpsInputValue] = useState<string>(String(frameRate));
	const [debouncedConnected, setDebouncedConnected] = useState(connected);
	const [debouncedIsFrameRendering, setDebouncedIsFrameRendering] =
		useState(isFrameRendering);

	useEffect(() => {
		setInputValue(String(step));
	}, [step]);
	useEffect(() => {
		setFpsInputValue(String(frameRate));
	}, [frameRate]);
	useEffect(() => {
		const handler = setTimeout(
			() => setDebouncedConnected(connected),
			connected ? 0 : 500,
		);
		return () => clearTimeout(handler);
	}, [connected]);
	useEffect(() => {
		const handler = setTimeout(
			() => setDebouncedIsFrameRendering(isFrameRendering),
			isFrameRendering ? 500 : 0,
		);
		return () => clearTimeout(handler);
	}, [isFrameRendering]);

	const handleSliderChange = (event: Event, newValue: number | number[]) => {
		setStep(newValue as number);
	};
	const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
		setInputValue(event.target.value);
		const numValue = Number(event.target.value);
		if (!isNaN(numValue) && numValue >= 0 && numValue <= length)
			setStep(numValue);
	};
	const handleInputBlur = () => {
		const numValue = Number(inputValue);
		if (isNaN(numValue) || numValue < 0 || numValue > length)
			setInputValue(String(step));
	};
	const handleBookmarkClick = (
		event: React.MouseEvent,
		frameNumber: number,
	) => {
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
		setFpsInputValue(event.target.value);
		const numValue = Number(event.target.value);
		if (!isNaN(numValue) && numValue > 0) setFrameRate(numValue);
	};
	const handleFpsInputBlur = () => {
		const numValue = Number(fpsInputValue);
		if (isNaN(numValue) || numValue <= 0) setFpsInputValue(String(frameRate));
	};

	const renderBookmarks = () =>
		Object.keys(bookmarks).map((key) => {
			const position = Number.parseInt(key);
			return (
				<Tooltip
					key={`bookmark-${position}`}
					title={bookmarks[position]}
					placement="top"
				>
					<BookmarkIndicator
						left={`${(position / (length - 1)) * 100}%`}
						onClick={(e) => handleBookmarkClick(e, position)}
					>
						<svg viewBox="0 0 12 7" xmlns="http://www.w3.org/2000/svg">
							<path d="M0 0 L12 0 L6 7 Z" />
						</svg>
					</BookmarkIndicator>
				</Tooltip>
			);
		});

	const renderSelectedFrames = () => {
		if (!selectedFrames.active || !selectedFrames.indices) return null;
		const frames = [...selectedFrames.indices].sort((a, b) => a - b);
		const blocks: { start: number; end: number }[] = [];
		if (frames.length > 0) {
			let currentBlock = { start: frames[0], end: frames[0] };
			for (let i = 1; i < frames.length; i++) {
				if (frames[i] === currentBlock.end + 1) {
					currentBlock.end = frames[i];
				} else {
					blocks.push(currentBlock);
					currentBlock = { start: frames[i], end: frames[i] };
				}
			}
			blocks.push(currentBlock);
		}
		return blocks.map((block, index) => (
			<Box
				key={`selected-block-${index}`}
				sx={{
					position: "absolute",
					left: `${(block.start / (length - 1)) * 100}%`,
					width: `${((block.end - block.start + 1) / (length - 1)) * 100}%`,
					height: "6px",
					backgroundColor: "primary.light",
					borderRadius: "3px",
					zIndex: 1,
					pointerEvents: "none",
					top: "50%",
					transform: "translateY(-50%)",
					opacity: 0.9,
				}}
			/>
		));
	};

	return (
		<Box
			sx={{
				width: "100%",
				padding: "4px 16px",
				backgroundColor: "background.paper",
				borderRadius: "8px",
				boxShadow: 1,
				display: "flex",
				alignItems: "center",
				gap: "12px",
				flexWrap: "wrap",
			}}
		>
			<TextField
				variant="outlined"
				size="small"
				value={inputValue}
				onChange={handleInputChange}
				onBlur={handleInputBlur}
				type="number"
				sx={{ ...compactTextFieldSx, width: 110 }}
				InputProps={{
					endAdornment: (
						<InputAdornment position="end">/{length - 1}</InputAdornment>
					),
				}}
			/>

			<SliderControlsWrapper>
				{renderSelectedFrames()}
				{renderBookmarks()}
				<Slider
					value={step}
					min={0}
					max={length - 1}
					step={1}
					onChange={handleSliderChange}
					sx={{
						color: selectedFrames.active ? "transparent" : "primary.main",
						"& .MuiSlider-rail": {
							height: 4,
							borderRadius: 2,
							backgroundColor: selectedFrames.active
								? theme.palette.action.disabled
								: theme.palette.action.hover,
							opacity: 1,
						},
						"& .MuiSlider-track": { height: 4, borderRadius: 2, transition: "none" },
						"& .MuiSlider-thumb": {
							width: 14,
							height: 14,
							backgroundColor: theme.palette.background.paper,
							border: "2px solid currentColor",
							color: selectedFrames.active
								? theme.palette.primary.light
								: "primary.main",
							transition: "none",
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
				sx={{ ...compactTextFieldSx, width: 80 }}
			/>

			<Tooltip
				title={selectedFrames.active ? "Hide Selections" : "Show Selections"}
			>
				<IconButton
					onClick={handleToggleSelectedFrames}
					size="small"
					color="primary"
				>
					{selectedFrames.active ? <VisibilityIcon /> : <VisibilityOffIcon />}
				</IconButton>
			</Tooltip>

			<Box
				sx={{
					display: "flex",
					alignItems: "center",
					justifyContent: "center",
					width: 24,
					height: 24,
				}}
			>
				<Tooltip
					title={
						!debouncedConnected
							? "Connecting..."
							: debouncedIsFrameRendering
								? "Waiting for data..."
								: "Connected"
					}
				>
					<Box sx={{ display: "flex" }}>
						{!debouncedConnected ? (
							<CircularProgress size={20} />
						) : debouncedIsFrameRendering ? (
							<CloudDownloadIcon
								sx={{ color: "#1976d2", fontSize: 22, ...waitingAnimation }}
							/>
						) : (
							<WifiIcon sx={{ color: "#4caf50", fontSize: 22 }} />
						)}
					</Box>
				</Tooltip>
			</Box>
		</Box>
	);
};

export default FrameProgressBar;
