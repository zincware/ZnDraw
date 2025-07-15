import { Box, Container, TextField, Tooltip } from "@mui/material";
import React, { useState, useEffect } from "react";
import "./progressbar.css";

import { FaRegBookmark } from "react-icons/fa";
import { IoMdCodeDownload } from "react-icons/io";
import { PiSelection, PiSelectionSlash } from "react-icons/pi";
import type { IndicesState } from "./utils";

interface JumpFrameProps {
	step: number;
	setStep: (step: number) => void;
	length: number;
}

const JumpFrame: React.FC<JumpFrameProps> = ({ step, setStep, length }) => {
	const handleBlur = (e: React.FocusEvent<HTMLInputElement>) => {
		if (e.target.value === "") {
			return;
		}
		const newStep = Number.parseInt(e.target.value, 10);
		if (newStep >= 0 && newStep < length) {
			setStep(newStep);
		} else {
			alert(`Invalid input. Please enter a number between 0 and ${length - 1}`);
		}
		e.target.value = "";
	};

	const handleKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
		if (e.key === "Enter") {
			e.preventDefault();
			e.currentTarget.blur();
		}
	};

	return (
		<TextField
			className="text-center user-select-none"
			placeholder={`${step === -1 ? length - 1 : step}/${length - 1}`}
			onBlur={handleBlur}
			onKeyDown={handleKeyDown}
			variant="outlined"
			size="small"
			style={{
				background: "transparent",
				zIndex: 1,
			}}
			InputProps={{
				style: {
					borderColor: "transparent",
					textAlign: "center",
				},
			}}
		/>
	);
};

interface ProgressBarProps {
	length: number;
	disabledFrames: number[];
	bookmarks: any;
	setBookmarks: any;
	step: number;
	setStep: (step: number) => void;
}

const ColoredTiles = ({
	length,
	disabledFrames,
	setStep,
	tickInterval,
}: {
	length: number;
	disabledFrames: number[];
	setStep: (step: number) => void;
	tickInterval: number;
}) => {
	const [disabledPositions, setdisabledPositions] = useState<number[]>([]);
	const [ticks, setTicks] = useState<number[]>([]);

	useEffect(() => {
		const disabledPositions = [...Array(length).keys()].filter((position) =>
			disabledFrames.includes(position),
		);
		setdisabledPositions(disabledPositions);
	}, [length, disabledFrames]);

	useEffect(() => {
		const ticks = [...Array(length).keys()].filter(
			(position) => position % tickInterval === 0,
		);
		setTicks(ticks);
	}, [length, tickInterval]);

	const onTileClick = (event: any) => {
		const rect = event.target.getBoundingClientRect();
		const x = event.clientX - rect.left;
		const position = Math.floor((x / rect.width) * length);
		setStep(position);
	};
	return (
		<>
			{ticks.map((position) => {
				const commonStyles = {
					left: `${(position / length) * 100}%`,
					width: `${100 / (length - 1)}%`,
					height: 25,
				};
				return (
					<div
						key={position}
						className={"position-absolute"}
						style={commonStyles}
					>
						<div className="progress-bar-tick-line bg-secondary" />
					</div>
				);
			})}
			<div
				className={`position-absolute ${disabledPositions.length === 0 ? "bg-body" : "bg-success"}`}
				style={{ width: "100%", height: 25 }}
				onClick={(e) => onTileClick(e)}
			/>

			{disabledPositions.map((position) => {
				const commonStyles = {
					left: `${(position / length) * 100}%`,
					width: `${100 / (length - 1)}%`,
					height: 25,
				};
				return (
					<div
						key={position}
						className={"position-absolute p-0 bg-body"}
						style={commonStyles}
					/>
				);
			})}
		</>
	);
};

interface FrameRateControlProps {
	frameRate: number;
	setFrameRate: (val: number) => void;
	isFrameRendering: boolean;
	connected: boolean;
}

const FrameRateControl: React.FC<FrameRateControlProps> = ({
	frameRate,
	setFrameRate,
	isFrameRendering,
	connected,
}) => {
	const [showLoadingIcon, setShowLoadingIcon] = useState(false);

	// Delay before showing ⧗ to avoid flickering
	useEffect(() => {
		let timeout: NodeJS.Timeout;
		if (isFrameRendering) {
			timeout = setTimeout(() => setShowLoadingIcon(true), 50);
		} else {
			setShowLoadingIcon(false);
			clearTimeout(timeout);
		}
		return () => clearTimeout(timeout);
	}, [isFrameRendering]);

	const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
		const val = Number.parseInt(e.target.value, 10);
		if (!isNaN(val)) {
			setFrameRate(Math.max(1, val)); // Ensure frame rate is at least 1
		}
	};

	return (
		<div
			className="d-flex justify-content-center align-items-center user-select-none"
			style={{ display: "flex", height: "100%" }}
		>
			<TextField
				type="number"
				value={frameRate}
				onChange={handleChange}
				InputProps={{
					inputProps: {
						min: 1,
						style: {
							textAlign: "center",
							fontSize: "10px",
						},
					},
				}}
				size="small"
				style={{
					width: "40px",
					height: "18px",
				}}
				title="Frame rate (1 = every frame, 2 = every 2nd frame, etc.)"
			/>

			{/* Frame loading icon */}
			<Tooltip title="Loading frame data...">
				<div
					style={{
						height: "100%",
						animation: showLoadingIcon ? "pulse 1s infinite" : "none",
						visibility: showLoadingIcon ? "visible" : "hidden",
						alignItems: "center",
						justifyContent: "center",
						display: "flex",
					}}
					title="Loading frame data..."
				>
					<IoMdCodeDownload />
				</div>
			</Tooltip>

			{/* Server connection spinner */}
			{!connected && (
				<Tooltip title="Not connected to server">
					<div
						style={{
							height: "25px",
							width: "25px",
						}}
					>
						<div
							className="spinner-border spinner-border-sm text-primary"
							role="status"
							style={{
								width: "15px",
								height: "15px",
							}}
						></div>
						<span className="visually-hidden">Loading...</span>
					</div>
				</Tooltip>
			)}
		</div>
	);
};

const Bookmarks = ({
	length,
	bookmarks,
	setBookmarks,
	setStep,
}: {
	length: number;
	bookmarks: { [key: number]: string };
	setBookmarks: (bookmarks: { [key: number]: string }) => void;
	setStep: (step: number) => void;
}) => {
	const handleBookmarkClick = (event: any, number: number) => {
		if (event.shiftKey) {
			const newBookmarks = { ...bookmarks };
			delete newBookmarks[number];
			setBookmarks(newBookmarks);
		} else {
			setStep(number);
		}
	};

	return (
		<>
			{Object.keys(bookmarks).map((key) => {
				const position = Number.parseInt(key);
				return (
					<Tooltip
						key={position}
						title={bookmarks[position]}
						slotProps={{
							tooltip: {
								style: { marginLeft: "0.375em" },
							},
						}}
					>
						<div
							className="position-absolute progress-bar-bookmark"
							style={{
								left: `${(position / length) * 100}%`,
								top: 7,
								marginLeft: "-0.375em",
							}}
						>
							<FaRegBookmark
								className="position-absolute"
								size={"0.75em"}
								onClick={(e) => handleBookmarkClick(e, position)}
							/>
						</div>
					</Tooltip>
				);
			})}
		</>
	);
};

const VLine = ({ length, step }: { length: number; step: number }) => {
	return (
		<div
			className="position-absolute p-0 progress-bar-v-line"
			style={{
				left: `${step === -1 ? ((length - 1) / length) * 100 : (step / length) * 100}%`,
				height: 25,
				width: "2px",
			}}
			// why do I need to overwrite the height and width here?
		/>
	);
};

const ProgressBar = ({
	length,
	disabledFrames,
	bookmarks,
	setBookmarks,
	step,
	setStep,
}: ProgressBarProps) => {
	const [tickInterval, setTickInterval] = useState<number>(1);

	useEffect(() => {
		setTickInterval(Math.floor(length / 100) + 1);
	}, [length]);

	return (
		<Box className="position-relative">
			<ColoredTiles
				length={length}
				disabledFrames={disabledFrames}
				setStep={setStep}
				tickInterval={tickInterval}
			/>
			<Bookmarks
				length={length}
				bookmarks={bookmarks}
				setBookmarks={setBookmarks}
				setStep={setStep}
			/>
			<VLine length={length} step={step} />
		</Box>
	);
};

interface FrameProgressBarProps {
	length: number;
	step: number;
	setStep: (step: number) => void;
	selectedFrames: IndicesState;
	bookmarks: any[]; // Replace with actual type if known
	setBookmarks: (bookmarks: any[]) => void; // Replace with actual type if known
	setSelectedFrames: (selectedFrames: IndicesState) => void;
	connected: boolean;
	frameRate: number;
	setFrameRate: (rate: number) => void;
	isFrameRendering: boolean;
}

const FrameProgressBar: React.FC<FrameProgressBarProps> = ({
	length,
	step,
	setStep,
	selectedFrames,
	bookmarks,
	setBookmarks,
	setSelectedFrames,
	connected,
	frameRate,
	setFrameRate,
	isFrameRendering,
}) => {
	const [linePosition, setLinePosition] = useState<number>(0);
	const [disabledFrames, setDisabledFrames] = useState<number[]>([]);
	const progressHandleParentRef = React.createRef<HTMLDivElement>();

	useEffect(() => {
		// disable frames are the frames that are not selected, if the selectedFrames is not empty
		if (selectedFrames.indices.size > 0 && selectedFrames.active) {
			const disabledFrames = [...Array(length).keys()].filter(
				(frame) => !selectedFrames.indices.has(frame),
			);
			setDisabledFrames(disabledFrames);
		} else {
			setDisabledFrames([]);
		}
	}, [selectedFrames, length]);

	const handleSelectionReset = () => {
		setSelectedFrames((prev) => ({
			indices: prev.indices,
			active: !prev.active,
		}));
	};

	useEffect(() => {
		// Calculate the linePosition based on the step, length, and window width
		setLinePosition(
			step === -1 ? ((length - 1) / length) * 100 : (step / length) * 100,
		);
	}, [step, length]);

	const handleMouseDown = (e) => {
		e.preventDefault();
		if (!progressHandleParentRef.current) {
			console.error(
				"progressHandleParentRef is not set - dragging will not work",
			);
			return;
		}
		// we need to store the parent rect in a variable because
		// once we move the mouse, the parent rect will change
		const parentRect = progressHandleParentRef.current.getBoundingClientRect();

		const handleMouseMove = (e) => {
			// compute the lineposition based on the mouse position and the length
			const x = e.clientX - parentRect.left;
			const position = Math.floor((x / parentRect.width) * length);
			setStep(Math.min(Math.max(position, 0), length - 1));
		};

		document.addEventListener("mousemove", handleMouseMove);

		document.addEventListener(
			"mouseup",
			() => {
				document.removeEventListener("mousemove", handleMouseMove);
			},
			{ once: true },
		);
	};

	return (
		<Container className="fixed-bottom px-0 py-0">
			<Box display="flex" width="100%">
				<Box width="16.67%">
					<Box>
						<Box
							className="d-flex bg-secondary justify-content-center align-items-center"
							style={{ height: 1 }}
						/>
					</Box>
					<Box>
						<Box
							className="d-flex bg-body justify-content-center align-items-center"
							style={{ height: 25 }}
						>
							<JumpFrame step={step} setStep={setStep} length={length} />
							<Tooltip title="toggle selection">
								<div>
									{" "}
									{selectedFrames.active ? (
										<PiSelection onClick={handleSelectionReset} />
									) : (
										<PiSelectionSlash onClick={handleSelectionReset} />
									)}
								</div>
							</Tooltip>
						</Box>
					</Box>
				</Box>
				<Box flexGrow={1}>
					<Box className="position-relative">
						<Box
							className="d-flex justify-content-center"
							ref={progressHandleParentRef}
						>
							<div
								className="handle"
								style={{ left: `${linePosition}%`, cursor: "pointer" }}
								onMouseDown={(e) => handleMouseDown(e)}
							>
								<div className="square" />
								<div className="triangle" />
							</div>
						</Box>
					</Box>
					<Box>
						<Box
							className="d-flex bg-secondary justify-content-center align-items-center"
							style={{ height: 1 }}
						/>
					</Box>

					<ProgressBar
						length={length}
						disabledFrames={disabledFrames}
						bookmarks={bookmarks}
						step={step}
						setStep={setStep}
						setBookmarks={setBookmarks}
					/>
				</Box>
				<Box width="8.33%">
					<Box>
						<Box
							className="d-flex bg-secondary justify-content-center align-items-center"
							style={{ height: 1 }}
						/>
					</Box>
					<Box>
						<Box
							className="d-flex bg-body justify-content-center align-items-center"
							style={{ height: 25 }}
						>
							<FrameRateControl
								frameRate={frameRate}
								setFrameRate={setFrameRate}
								isFrameRendering={isFrameRendering}
								connected={connected}
							/>
						</Box>
					</Box>
				</Box>
			</Box>
		</Container>
	);
};

export default FrameProgressBar;
