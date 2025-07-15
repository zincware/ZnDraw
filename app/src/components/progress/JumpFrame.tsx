import { TextField } from "@mui/material";
import React from "react";

interface JumpFrameProps {
	step: number;
	setStep: (step: number) => void;
	length: number;
}

export const JumpFrame: React.FC<JumpFrameProps> = ({ step, setStep, length }) => {
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