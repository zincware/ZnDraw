import {
  Box,
  Tooltip,
  TextField,
  ClickAwayListener,
} from "@mui/material";
import { useEffect, useState, useMemo, memo, FC } from "react";
import BookmarkIcon from "@mui/icons-material/Bookmark";

interface SingleBookmarkProps {
  frame: number;
  label: string;
  position: string;
  isCurrentFrame: boolean;
  onBookmarkClick?: (frame: number, label: string) => void;
  onBookmarkDelete?: (frame: number) => void;
  onBookmarkEdit?: (frame: number, newLabel: string) => void;
}

const SingleBookmark: FC<SingleBookmarkProps> = memo(
  ({
    frame,
    label,
    position,
    isCurrentFrame,
    onBookmarkClick,
    onBookmarkDelete,
    onBookmarkEdit,
  }) => {
    const [isEditing, setIsEditing] = useState(false);
    const [editValue, setEditValue] = useState(label);

    const handleEditStart = () => {
      setIsEditing(true);
      setEditValue(label); // Reset edit value to current label on edit start
    };

    const handleEditSubmit = () => {
      const trimmedValue = editValue.trim();
      // UX Improvement: Only fire the callback if the value is not empty
      // and has actually changed from the original label.
      if (trimmedValue && trimmedValue !== label && onBookmarkEdit) {
        onBookmarkEdit(frame, trimmedValue);
      }
      setIsEditing(false);
    };

    const handleEditCancel = () => {
      setIsEditing(false);
    };

    if (isEditing) {
      return (
        <ClickAwayListener onClickAway={handleEditSubmit}>
          <Box
            sx={{
              position: "absolute",
              left: position,
              transform: "translateX(-50%)",
              top: "-8px",
              pointerEvents: "auto",
              zIndex: 10,
            }}
          >
            <TextField
              value={editValue}
              onChange={(e) => setEditValue(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") {
                  handleEditSubmit();
                } else if (e.key === "Escape") {
                  handleEditCancel();
                }
                e.stopPropagation();
              }}
              autoFocus
              size="small"
              sx={{
                width: "150px",
                "& .MuiInputBase-input": {
                  fontSize: "0.875rem",
                  padding: "4px 8px",
                },
                backgroundColor: "background.paper",
              }}
            />
          </Box>
        </ClickAwayListener>
      );
    }

    return (
      <Tooltip
        title={
          <Box>
            <Box sx={{ fontWeight: "bold", mb: 0.5 }}>{label}</Box>
            <Box sx={{ fontSize: "0.75rem", opacity: 0.8 }}>
              Frame: {frame}
            </Box>
            <Box sx={{ fontSize: "0.7rem", opacity: 0.6, mt: 0.5, fontStyle: "italic" }}>
              Double-click to edit, Shift+Click to delete
            </Box>
          </Box>
        }
        arrow
        placement="top"
        enterDelay={300}
        enterNextDelay={300}
      >
        <Box
          onClick={(e) => {
            if (e.shiftKey) {
              e.preventDefault();
              e.stopPropagation();
              onBookmarkDelete?.(frame);
            } else {
              onBookmarkClick?.(frame, label);
            }
          }}
          onDoubleClick={(e) => {
            e.preventDefault();
            e.stopPropagation();
            handleEditStart();
          }}
          sx={{
            position: "absolute",
            left: position,
            transform: "translateX(-50%)",
            cursor: "pointer",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            transition: "all 0.2s ease",
            pointerEvents: "auto",
            "&:hover": {
              transform: "translateX(-50%) scale(1.2)",
            },
          }}
        >
          <BookmarkIcon
            sx={{
              fontSize: 20,
              color: isCurrentFrame ? "secondary.main" : "primary.main",
              "& path": {
                fill: "none",
                stroke: "currentColor",
                strokeWidth: 1.5,
              },
              filter: isCurrentFrame
                ? "drop-shadow(0 0 4px rgba(156, 39, 176, 0.6))"
                : "drop-shadow(0 1px 2px rgba(0,0,0,0.2))",
              transition: "color 0.2s ease",
            }}
          />
        </Box>
      </Tooltip>
    );
  }
);

interface BookmarkLayerProps {
  bookmarks: Record<number, string> | null;
  frameCount: number;
  currentFrame: number;
  containerWidth: number;
  onBookmarkClick?: (frame: number, label: string) => void;
  onBookmarkDelete?: (frame: number) => void;
  onBookmarkEdit?: (frame: number, newLabel: string) => void;
}

const BookmarkLayer: FC<BookmarkLayerProps> = ({
  bookmarks,
  frameCount,
  currentFrame,
  containerWidth,
  onBookmarkClick,
  onBookmarkDelete,
  onBookmarkEdit,
}) => {
  // Performance Optimization: Memoize the conversion of the bookmarks object
  // to an array. This calculation will only re-run when the `bookmarks`
  // prop itself changes, not on every render.
  const bookmarkEntries = useMemo(() => {
    if (!bookmarks) return [];
    return Object.entries(bookmarks).map(([frame, label]) => ({
      frame: parseInt(frame),
      label,
    }));
  }, [bookmarks]);

  const getMarkerPosition = (frame: number): string => {
    if (frameCount <= 1) return "0%";
    const percentage = (frame / (frameCount - 1)) * 100;
    return `${percentage}%`;
  };

  if (!bookmarks || bookmarkEntries.length === 0 || containerWidth === 0) {
    return null;
  }

  return (
    <Box
      sx={{
        position: "absolute",
        top: "-20px",
        left: 0,
        right: 0,
        height: 16,
        pointerEvents: "none", // Container doesn't block clicks
        zIndex: 2,
      }}
    >
      {bookmarkEntries.map(({ frame, label }) => (
        <SingleBookmark
          key={frame}
          frame={frame}
          label={label}
          position={getMarkerPosition(frame)}
          isCurrentFrame={frame === currentFrame}
          onBookmarkClick={onBookmarkClick}
          onBookmarkDelete={onBookmarkDelete}
          onBookmarkEdit={onBookmarkEdit}
        />
      ))}
    </Box>
  );
};

export default memo(BookmarkLayer);
