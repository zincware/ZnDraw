import React from "react";

export const VLine = ({ length, step }: { length: number; step: number }) => {
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