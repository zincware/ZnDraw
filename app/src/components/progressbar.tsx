import { Box, Tooltip } from "@mui/material";
import React, { useEffect, useState } from "react";
import { PiSelection, PiSelectionSlash } from "react-icons/pi";
import type { IndicesState } from "../utils";
import { FrameRateControl } from "./progress/FrameRateControl";
import { JumpFrame } from "./progress/JumpFrame";
import { ProgressBar } from "./progress/ProgressBar";
import "./progressbar.css";

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
		<Box
			sx={{
				position: "fixed",
				bottom: 0,
				left: 0,
				right: 0,
				bgcolor: "background.paper",
				display: "flex",
				width: "100%",
			}}
		>
			<Box sx={{ width: "16.67%" }}>
				<Box sx={{ height: 1, bgcolor: "divider" }} />
				<Box sx={{ height: 25, display: "flex", alignItems: "center", justifyContent: "center" }}>
					<JumpFrame step={step} setStep={setStep} length={length} />
					<Tooltip title="toggle selection">
						<Box sx={{ cursor: "pointer" }}>
							{" "}
							{selectedFrames.active ? (
								<PiSelection onClick={handleSelectionReset} />
							) : (
								<PiSelectionSlash onClick={handleSelectionReset} />
							)}
						</Box>
					</Tooltip>
				</Box>
			</Box>
			<Box sx={{ flexGrow: 1 }}>
				<Box sx={{ position: "relative" }}>
					<Box sx={{ display: "flex", justifyContent: "center" }} ref={progressHandleParentRef}>
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
				<Box sx={{ height: 1, bgcolor: "divider" }} />
				<ProgressBar
					length={length}
					disabledFrames={disabledFrames}
					bookmarks={bookmarks}
					step={step}
					setStep={setStep}
					setBookmarks={setBookmarks}
				/>
			</Box>
			<Box sx={{ width: "8.33%" }}>
				<Box sx={{ height: 1, bgcolor: "divider" }} />
				<Box sx={{ height: 25, display: "flex", alignItems: "center", justifyContent: "center" }}>
					<FrameRateControl
						frameRate={frameRate}
						setFrameRate={setFrameRate}
						isFrameRendering={isFrameRendering}
						connected={connected}
					/>
				</Box>
			</Box>
		</Box>
	);
};

export default FrameProgressBar;