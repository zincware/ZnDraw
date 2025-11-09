import { Box, Chip, Menu, MenuItem, Typography } from "@mui/material";
import BrushIcon from "@mui/icons-material/Brush";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import { useState, MouseEvent } from "react";
import { useAppStore } from "../../store";

/**
 * DrawingIndicator shows which curve is currently active for drawing
 * and allows quick switching between curves via dropdown
 */
export default function DrawingIndicator() {
  // Use individual selectors to prevent unnecessary re-renders
  const mode = useAppStore((state) => state.mode);
  const activeCurveForDrawing = useAppStore((state) => state.activeCurveForDrawing);
  const setActiveCurveForDrawing = useAppStore((state) => state.setActiveCurveForDrawing);
  const geometries = useAppStore((state) => state.geometries);
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);

  // Only show if drawing mode is active
  if (mode !== 'drawing' || !activeCurveForDrawing) {
    return null;
  }

  // Get all active curves
  const activeCurves = Object.entries(geometries)
    .filter(([_, g]) => g.type === "Curve" && g.data?.active !== false)
    .map(([key, _]) => key);

  const handleClick = (event: MouseEvent<HTMLDivElement>) => {
    // Only show dropdown if there are multiple curves
    if (activeCurves.length > 1) {
      setAnchorEl(event.currentTarget);
    }
  };

  const handleClose = () => {
    setAnchorEl(null);
  };

  const handleCurveSelect = (curveKey: string) => {
    setActiveCurveForDrawing(curveKey);
    handleClose();
  };

  const isMenuOpen = Boolean(anchorEl);
  const showDropdownIcon = activeCurves.length > 1;

  return (
    <>
      <Box
        sx={{
          position: "fixed",
          bottom: 80,
          right: 16,
          zIndex: (theme) => theme.zIndex.tooltip,
        }}
      >
        <Chip
          icon={<BrushIcon />}
          label={
            <Box sx={{ display: "flex", alignItems: "center", gap: 0.5 }}>
              <Typography variant="body2">Drawing on: {activeCurveForDrawing}</Typography>
              {showDropdownIcon && <ExpandMoreIcon fontSize="small" />}
            </Box>
          }
          color="primary"
          size="medium"
          onClick={handleClick}
          sx={{
            cursor: showDropdownIcon ? "pointer" : "default",
            "& .MuiChip-label": {
              px: 1,
            },
          }}
        />
      </Box>

      {/* Dropdown menu for switching curves */}
      {showDropdownIcon && (
        <Menu
          anchorEl={anchorEl}
          open={isMenuOpen}
          onClose={handleClose}
          anchorOrigin={{
            vertical: "top",
            horizontal: "right",
          }}
          transformOrigin={{
            vertical: "bottom",
            horizontal: "right",
          }}
        >
          {activeCurves.map((curveKey) => (
            <MenuItem
              key={curveKey}
              selected={curveKey === activeCurveForDrawing}
              onClick={() => handleCurveSelect(curveKey)}
            >
              {curveKey}
            </MenuItem>
          ))}
        </Menu>
      )}
    </>
  );
}
