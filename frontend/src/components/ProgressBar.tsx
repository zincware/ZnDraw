import CloudDownloadIcon from "@mui/icons-material/CloudDownload";
import FilterListIcon from "@mui/icons-material/FilterList";
import SyncIcon from "@mui/icons-material/Sync";
import SyncDisabledIcon from "@mui/icons-material/SyncDisabled";
import WifiIcon from "@mui/icons-material/Wifi";
import {
	Box,
	CircularProgress,
	IconButton,
	Slider,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import { useEffect, useRef, useState } from "react";
import { LAYOUT_CONSTANTS } from "../constants/layout";
import { useStepControl } from "../hooks/useStepControl";
import { socket } from "../socket";
import { useAppStore } from "../store";
import BookmarkLayer from "./BookmarkLayer";
import FrameSelectionInput from "./FrameSelectionInput";
import SelectionLayer from "./SelectionLayer";
import SelectionTrackOverlay from "./SelectionTrackOverlay";

const FrameProgressBar = () => {
	const [isEditing, setIsEditing] = useState(false);
	const [inputValue, setInputValue] = useState("0");
	const sliderContainerRef = useRef<HTMLDivElement>(null);
	const [sliderWidth, setSliderWidth] = useState(0);
	const [previewRange, setPreviewRange] = useState<[number, number] | null>(
		null,
	);

	const currentFrame = useAppStore((state) => state.currentFrame);
	const frameCount = useAppStore((state) => state.frameCount);
	const isConnected = useAppStore((state) => state.isConnected);
	const isLoading = useAppStore((state) => state.isLoading);
	const skipFrames = useAppStore((state) => state.skipFrames);
	const setSkipFrames = useAppStore((state) => state.setSkipFrames);
	const frame_selection = useAppStore((state) => state.frame_selection);
	const bookmarks = useAppStore((state) => state.bookmarks);
	const addBookmark = useAppStore((state) => state.addBookmark);
	const deleteBookmark = useAppStore((state) => state.deleteBookmark);
	const frameSelectionEnabled = useAppStore(
		(state) => state.frameSelectionEnabled,
	);
	const setFrameSelectionEnabled = useAppStore(
		(state) => state.setFrameSelectionEnabled,
	);
	const synchronizedMode = useAppStore((state) => state.synchronizedMode);
	const setSynchronizedMode = useAppStore((state) => state.setSynchronizedMode);
	const getIsFetching = useAppStore((state) => state.getIsFetching);
	const playing = useAppStore((state) => state.playing);
	const setPlaying = useAppStore((state) => state.setPlaying);
	const showSnackbar = useAppStore((state) => state.showSnackbar);

	// Long-fetch spinner: poll getIsFetching() (derived getter, not reactive),
	// show spinner after 500ms of continuous fetching, snackbar after 2s.
	const [showFetchSpinner, setShowFetchSpinner] = useState(false);
	const fetchStartRef = useRef<number | null>(null);
	const snackbarShownRef = useRef(false);

	// Reset tracking on frame change so timers restart
	useEffect(() => {
		fetchStartRef.current = null;
		snackbarShownRef.current = false;
		setShowFetchSpinner(false);
	}, [currentFrame]);

	useEffect(() => {
		const interval = setInterval(() => {
			const fetching = getIsFetching();
			if (fetching) {
				if (fetchStartRef.current === null) {
					fetchStartRef.current = Date.now();
				}
				const elapsed = Date.now() - fetchStartRef.current;
				if (elapsed >= 500) {
					setShowFetchSpinner(true);
				}
				if (elapsed >= 2000 && !snackbarShownRef.current) {
					snackbarShownRef.current = true;
					showSnackbar("Loading frame from source...", "info");
				}
			} else {
				if (fetchStartRef.current !== null) {
					fetchStartRef.current = null;
					snackbarShownRef.current = false;
					setShowFetchSpinner(false);
				}
			}
		}, 100);
		return () => clearInterval(interval);
	}, [getIsFetching, showSnackbar]);

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
		// Stop playback when scrubbing
		if (playing) {
			setPlaying(false);
		}
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
		if (
			frameSelectionEnabled &&
			(!frame_selection || frame_selection.length === 0)
		) {
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
		setInputValue((currentFrame + 1).toString());
	};

	const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
		setInputValue(event.target.value);
	};

	const handleInputSubmit = () => {
		const displayFrame = Number.parseInt(inputValue, 10);
		const newFrame = displayFrame - 1;
		if (!isNaN(newFrame) && newFrame >= 0 && newFrame <= frameCount - 1) {
			// Stop playback when manually setting frame
			if (playing) {
				setPlaying(false);
			}
			setStep(newFrame);
		} else {
			setInputValue((currentFrame + 1).toString());
		}
		setIsEditing(false);
	};

	const handleInputKeyPress = (event: React.KeyboardEvent) => {
		if (event.key === "Enter") {
			handleInputSubmit();
		} else if (event.key === "Escape") {
			setInputValue((currentFrame + 1).toString());
			setIsEditing(false);
		}
	};

	const handleInputBlur = () => {
		handleInputSubmit();
	};

	const handleSkipFramesChange = (
		event: React.ChangeEvent<HTMLInputElement>,
	) => {
		const value = Number.parseInt(event.target.value, 10);
		if (!isNaN(value) && value > 0) {
			setSkipFrames(value);
		}
	};

	// Bookmark click handler - stops playback and jumps to frame
	const handleBookmarkClick = (frame: number) => {
		if (playing) {
			setPlaying(false);
		}
		setStep(frame);
	};

	const renderConnectionStatus = () => {
		if (showFetchSpinner) {
			return (
				<Tooltip title="Loading frame from source...">
					<Box
						sx={{
							display: "flex",
							alignItems: "center",
							justifyContent: "center",
						}}
					>
						<CircularProgress size={20} />
					</Box>
				</Tooltip>
			);
		}
		if (isLoading) {
			return (
				<CloudDownloadIcon
					sx={{ color: "#1976d2", fontSize: 22, ...waitingAnimation }}
				/>
			);
		}
		if (isConnected) {
			return <WifiIcon sx={{ color: "#4caf50", fontSize: 22 }} />;
		}
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
						{currentFrame + 1} / {frameCount}
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
					onBookmarkClick={handleBookmarkClick}
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
					previewRange={previewRange}
				/>

				{/* Frame Selection Input - Shift+Click/Drag to select frame ranges */}
				<FrameSelectionInput
					frameCount={frameCount}
					containerWidth={sliderWidth}
					onPreviewChange={setPreviewRange}
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
					valueLabelFormat={(v) => v + 1}
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
								? "Disable frame selection filter (Escape to clear)"
								: "Enable frame selection filter (Escape to clear)"
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
							color: synchronizedMode ? "primary.main" : "text.secondary",
							...(getIsFetching() && synchronizedMode ? waitingAnimation : {}),
							"&:hover": {
								backgroundColor: "action.hover",
							},
						}}
					>
						{synchronizedMode ? (
							<SyncIcon fontSize="small" />
						) : (
							<SyncDisabledIcon fontSize="small" />
						)}
					</IconButton>
				</Tooltip>
				{renderConnectionStatus()}
			</Box>
		</Box>
	);
};

export default FrameProgressBar;
