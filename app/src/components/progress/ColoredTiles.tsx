import React, { useEffect, useState } from "react";

export const ColoredTiles = ({
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
