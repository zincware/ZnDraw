import { Box, Tooltip } from "@mui/material";
import { useEffect } from "react";
import BookmarkIcon from '@mui/icons-material/Bookmark';

interface BookmarkLayerProps {
  bookmarks: Record<number, string> | null;
  frameCount: number;
  currentFrame: number;
  containerWidth: number;
  onBookmarkClick?: (frame: number, label: string) => void;
}

const BookmarkLayer: React.FC<BookmarkLayerProps> = ({
  bookmarks,
  frameCount,
  currentFrame,
  containerWidth,
  onBookmarkClick,
}) => {
  useEffect(() => {
    console.log("Bookmarks updated:", bookmarks);
    console.log("Bookmark container width:", containerWidth);
  }, [bookmarks, containerWidth]);

  const getMarkerPosition = (frame: number): string => {
    if (frameCount <= 1) return '0%';

    // Calculate percentage position - matches MUI Slider calculation
    const percentage = (frame / (frameCount - 1)) * 100;

    return `${percentage}%`;
  };

  if (!bookmarks || Object.keys(bookmarks).length === 0 || containerWidth === 0) {
    console.log("Not rendering bookmark layer:", { bookmarks, containerWidth });
    return null;
  }

  const bookmarkEntries = Object.entries(bookmarks).map(([frame, label]) => ({
    frame: parseInt(frame),
    label,
  }));

  console.log("Rendering", bookmarkEntries.length, "bookmarks");

  return (
    <Box
      sx={{
        position: 'absolute',
        top: '-20px',
        left: 0,
        right: 0,
        height: 16,
        pointerEvents: 'none', // Container doesn't block clicks
        zIndex: 2,
      }}
    >
      {bookmarkEntries.map(({ frame, label }) => {
        const isCurrentFrame = frame === currentFrame;
        const leftPercent = getMarkerPosition(frame);

        return (
          <Tooltip
            key={frame}
            title={
              <Box>
                <Box sx={{ fontWeight: 'bold', mb: 0.5 }}>{label}</Box>
                <Box sx={{ fontSize: '0.75rem', opacity: 0.8 }}>
                  Frame: {frame}
                </Box>
              </Box>
            }
            arrow
            placement="top"
            enterDelay={300}
            enterNextDelay={300}
          >
            <Box
              onClick={() => onBookmarkClick?.(frame, label)}
              sx={{
                position: 'absolute',
                left: leftPercent,
                transform: 'translateX(-50%)',
                cursor: 'pointer',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                transition: 'all 0.2s ease',
                pointerEvents: 'auto', // Only bookmark icon is clickable
                '&:hover': {
                  transform: 'translateX(-50%) scale(1.2)',
                },
              }}
            >
              <BookmarkIcon
                sx={{
                  fontSize: 20,
                  color: isCurrentFrame ? 'secondary.main' : 'primary.main',
                  fontWeight: 'normal',
                  '& path': {
                    fill: 'none',
                    stroke: 'currentColor',
                    strokeWidth: 1.5,
                  },
                  filter: isCurrentFrame
                    ? 'drop-shadow(0 0 4px rgba(156, 39, 176, 0.6))'
                    : 'drop-shadow(0 1px 2px rgba(0,0,0,0.2))',
                  transition: 'color 0.2s ease',
                }}
              />
            </Box>
          </Tooltip>
        );
      })}
    </Box>
  );
};

export default BookmarkLayer;
