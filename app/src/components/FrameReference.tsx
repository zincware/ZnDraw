// FrameReference.tsx
import React from "react";
import { Chip } from "@mui/material";
import MovieFilterIcon from "@mui/icons-material/MovieFilter";
import ViewInArIcon from "@mui/icons-material/ViewInAr";
import { useAppStore } from "../store";
import LocationOnIcon from "@mui/icons-material/LocationOn";
import { useStepControl } from "../hooks/useStepControl";

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
          label={`${frame}`}
          clickable
          color="primary"
          onClick={handleClick}
        />
      ) : (
        <span>@{frame}</span>
      )}
    </>
  );
};
