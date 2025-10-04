import { Box } from "@mui/material";
import { useEffect } from "react";

interface SelectionLayerProps {
  selectedFrames: number[] | null;
  frameCount: number;
  currentFrame: number;
  containerWidth: number;
}

const SelectionLayer: React.FC<SelectionLayerProps> = ({
  selectedFrames,
  frameCount,
  currentFrame,
  containerWidth,
}) => {
  useEffect(() => {
    console.log("Selected frames updated:", selectedFrames);
    console.log("Container width:", containerWidth);
  }, [selectedFrames, containerWidth]);

  const getMarkerPosition = (frame: number): string => {
    if (frameCount <= 1) return "0%";

    // Calculate percentage position - matches MUI Slider calculation
    const percentage = (frame / (frameCount - 1)) * 100;

    return `${percentage}%`;
  };

  if (!selectedFrames || selectedFrames.length === 0 || containerWidth === 0) {
    console.log("Not rendering selection layer:", {
      selectedFrames,
      containerWidth,
    });
    return null;
  }

  console.log("Rendering", selectedFrames.length, "selection markers");

  return (
    <Box
      sx={{
        position: "absolute",
        top: "50%",
        left: 0,
        right: 0,
        transform: "translateY(-50%)",
        height: 8,
        pointerEvents: "none",
        zIndex: 0, // Below bookmarks and slider
      }}
    >
      {selectedFrames.map((frame) => {
        const isCurrentFrame = frame === currentFrame;
        const leftPercent = getMarkerPosition(frame);

        return (
          <Box
            key={frame}
            sx={{
              position: "absolute",
              left: leftPercent,
              transform: "translateX(-50%)",
              width: 2,
              height: "100%",
              backgroundColor: isCurrentFrame
                ? "secondary.main"
                : "primary.main",
              opacity: isCurrentFrame ? 0.5 : 0.25,
              pointerEvents: "none",
            }}
          />
        );
      })}
    </Box>
  );
};

export default SelectionLayer;
