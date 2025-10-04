// FrameReference.tsx
import React from "react";
import { Chip } from "@mui/material";
import MovieFilterIcon from "@mui/icons-material/MovieFilter";
import ViewInArIcon from "@mui/icons-material/ViewInAr";
import { useAppStore } from "../store";
import LocationOnIcon from "@mui/icons-material/LocationOn";

interface FrameReferenceProps {
  frame: number;
}

export const FrameReference: React.FC<FrameReferenceProps> = ({ frame }) => {
  const { setCurrentFrame, frameCount } = useAppStore();

  const handleClick = () => {
    if (frame >= 0 && frame < frameCount) {
      setCurrentFrame(frame);
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
