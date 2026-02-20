import LocationOnIcon from "@mui/icons-material/LocationOn";
import MovieFilterIcon from "@mui/icons-material/MovieFilter";
import ViewInArIcon from "@mui/icons-material/ViewInAr";
import { Chip } from "@mui/material";
// FrameReference.tsx
import type React from "react";
import { useStepControl } from "../hooks/useStepControl";
import { useAppStore } from "../store";

interface FrameReferenceProps {
	frame: number;
}

export const FrameReference: React.FC<FrameReferenceProps> = ({ frame }) => {
	// Use individual selectors to prevent unnecessary re-renders
	const frameCount = useAppStore((state) => state.frameCount);
	const { setStep } = useStepControl();

	const handleClick = () => {
		if (frame >= 0 && frame < frameCount) {
			setStep(frame);
		}
	};

	return (
		<>
			{frame >= 0 && frame < frameCount ? (
				<Chip
					icon={<LocationOnIcon />}
					label={`${frame + 1}`}
					clickable
					color="primary"
					onClick={handleClick}
				/>
			) : (
				<span>@{frame + 1}</span>
			)}
		</>
	);
};
