import { Box } from "@mui/material";
import { useEffect } from "react";

interface SelectionTrackOverlayProps {
  selectedFrames: number[] | null;
  frameCount: number;
  containerWidth: number;
  enabled: boolean;
  currentFrame: number;
}

// Convert array of frames to contiguous ranges
const framesToRanges = (frames: number[]): Array<[number, number]> => {
  if (frames.length === 0) return [];

  const sorted = [...frames].sort((a, b) => a - b);
  const ranges: Array<[number, number]> = [];
  let start = sorted[0];
  let end = sorted[0];

  for (let i = 1; i < sorted.length; i++) {
    if (sorted[i] === end + 1) {
      end = sorted[i];
    } else {
      ranges.push([start, end]);
      start = sorted[i];
      end = sorted[i];
    }
  }
  ranges.push([start, end]);

  return ranges;
};

const SelectionTrackOverlay: React.FC<SelectionTrackOverlayProps> = ({
  selectedFrames,
  frameCount,
  containerWidth,
  enabled,
  currentFrame,
}) => {
  useEffect(() => {
    console.log("SelectionTrackOverlay:", { selectedFrames, enabled, containerWidth, currentFrame });
  }, [selectedFrames, enabled, containerWidth, currentFrame]);

  const getPercentPosition = (frame: number): number => {
    if (frameCount <= 1) return 0;
    return (frame / (frameCount - 1)) * 100;
  };

  if (containerWidth === 0) {
    return null;
  }

  // If disabled, show normal progress track (blue from 0 to current, grey from current to end)
  if (!enabled) {
    const segments = [
      { start: 0, end: currentFrame, selected: true },
      { start: currentFrame, end: frameCount - 1, selected: false },
    ];

    return (
      <Box
        sx={{
          position: 'absolute',
          top: 0,
          left: 0,
          right: 0,
          bottom: 0,
          pointerEvents: 'none',
          zIndex: 1,
        }}
      >
        {segments.map((segment, idx) => {
          const leftPercent = getPercentPosition(segment.start);
          const isLastSegment = idx === segments.length - 1;
          const visualEnd = isLastSegment ? frameCount - 1 : segments[idx + 1].start;
          const rightPercent = getPercentPosition(visualEnd);
          const widthPercent = rightPercent - leftPercent;

          if (widthPercent <= 0) return null;

          return (
            <Box
              key={idx}
              sx={{
                position: 'absolute',
                left: `${leftPercent}%`,
                width: `${widthPercent}%`,
                top: '50%',
                transform: 'translateY(-50%)',
                height: 4,
                backgroundColor: segment.selected ? 'primary.main' : 'action.disabledBackground',
                borderRadius: '2px',
              }}
            />
          );
        })}
      </Box>
    );
  }

  // Enabled mode: show segmented view based on selected frames
  if (!selectedFrames || selectedFrames.length === 0) {
    return null;
  }

  const ranges = framesToRanges(selectedFrames);

  // Generate segments for the entire timeline
  // Selected ranges = blue, gaps between = grey
  const segments: Array<{ start: number; end: number; selected: boolean }> = [];

  let currentPos = 0;
  ranges.forEach(([rangeStart, rangeEnd]) => {
    // Add grey segment before this range if there's a gap
    if (rangeStart > currentPos) {
      segments.push({
        start: currentPos,
        end: rangeStart - 1,
        selected: false,
      });
    }

    // Add blue segment for this range
    segments.push({
      start: rangeStart,
      end: rangeEnd,
      selected: true,
    });

    currentPos = rangeEnd + 1;
  });

  // Add grey segment after last range if needed
  if (currentPos < frameCount - 1) {
    segments.push({
      start: currentPos,
      end: frameCount - 1,
      selected: false,
    });
  }

  return (
    <Box
      sx={{
        position: 'absolute',
        top: 0,
        left: 0,
        right: 0,
        bottom: 0,
        pointerEvents: 'none',
        zIndex: 1, // Above the slider rail but below thumb
      }}
    >
      {segments.map((segment, idx) => {
        const leftPercent = getPercentPosition(segment.start);

        // Determine the visual end position for this segment
        // Extend to the start of the next segment, or to the end of the timeline
        const isLastSegment = idx === segments.length - 1;
        let visualEnd: number;

        if (isLastSegment) {
          visualEnd = frameCount - 1;
        } else {
          // Extend to just before the next segment starts
          visualEnd = segments[idx + 1].start;
        }

        const rightPercent = getPercentPosition(visualEnd);
        const widthPercent = rightPercent - leftPercent;

        return (
          <Box
            key={idx}
            sx={{
              position: 'absolute',
              left: `${leftPercent}%`,
              width: `${widthPercent}%`,
              top: '50%',
              transform: 'translateY(-50%)',
              height: 4,
              backgroundColor: segment.selected ? 'primary.main' : 'grey.400',
              borderRadius: '2px',
            }}
          />
        );
      })}
    </Box>
  );
};

export default SelectionTrackOverlay;
