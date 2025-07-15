import { Box } from "@mui/material";
import React, { useEffect, useState } from "react";
import { Bookmarks } from "./Bookmarks";
import { ColoredTiles } from "./ColoredTiles";
import { VLine } from "./VLine";

interface ProgressBarProps {
	length: number;
	disabledFrames: number[];
	bookmarks: any;
	setBookmarks: any;
	step: number;
	setStep: (step: number) => void;
}

export const ProgressBar = ({
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
			{/* <ColoredTiles
				length={length}
				disabledFrames={disabledFrames}
				setStep={setStep}
				tickInterval={tickInterval}
			/> */}
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