import { Box, Typography, Paper, useTheme } from "@mui/material";
import { useAppStore } from "../../store";
import { Rnd } from "react-rnd";
import { usePropertyInspectorSettings } from "../../hooks/usePropertyInspectorSettings";
import { formatPropertyValue } from "../../utils/propertyFormatting";

/**
 * StaticInfoBox - Displays general scene information in a draggable window.
 * Default position is top-right corner.
 * Shows:
 * - Number of particles in the current frame
 * - Current frame number
 * - Total frame count
 * - Selection count (if any particles are selected)
 */
export default function StaticInfoBox() {
  const theme = useTheme();
  const {
    showInfoBoxes,
    particleCount,
    selection,
    playing,
    fps,
    frameLoadTime,
  } = useAppStore();

  // Performance optimization: only fetch global properties when boxes are visible
  const {
    enabledProperties,
    propertyValues,
    isEnabled,
  } = usePropertyInspectorSettings({
    category: "global",
    enabled: showInfoBoxes,
  });

  if (!showInfoBoxes) return null;

  const selectionCount = selection?.length || 0;

  return (
    <Rnd
      default={{
        x: window.innerWidth - 280,
        y: 84,
        width: 240,
        height: "auto",
      }}
      minWidth={200}
      minHeight={100}
      bounds=".drag-boundary-container"
      dragHandleClassName="info-drag-handle"
      style={{ zIndex: 1300 }}
      enableResizing={false}
    >
      <Paper
        elevation={3}
        sx={{
          width: "100%",
          height: "100%",
          display: "flex",
          flexDirection: "column",
          backgroundColor: theme.palette.mode === "dark" 
            ? "rgba(0, 0, 0, 0.8)" 
            : "rgba(255, 255, 255, 0.9)",
          backdropFilter: "blur(10px)",
          border: `1px solid ${theme.palette.mode === "dark" 
            ? "rgba(255, 255, 255, 0.1)" 
            : "rgba(0, 0, 0, 0.1)"}`,
        }}
      >
        <Box
          className="info-drag-handle"
          sx={{
            padding: 2,
            cursor: "move",
          }}
        >
          <Typography
            variant="h6"
            sx={{
              color: theme.palette.text.primary,
              fontWeight: "bold",
              marginBottom: 1,
              fontSize: "0.9rem",
            }}
          >
            Scene Info
          </Typography>
          <Box sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}>
            <Typography variant="body2" sx={{ color: theme.palette.text.secondary }}>
              Particles: <strong style={{ color: theme.palette.text.primary }}>{particleCount}</strong>
            </Typography>
            {selectionCount > 0 && (
              <Typography variant="body2" sx={{ color: theme.palette.warning.main }}>
                Selected: <strong>{selectionCount}</strong>
              </Typography>
            )}
            {playing && fps !== null && (
              <Typography variant="body2" sx={{ color: theme.palette.success.main }}>
                FPS: <strong>{fps.toFixed(1)}</strong>
              </Typography>
            )}
            {!playing && frameLoadTime !== null && (
              <Typography variant="body2" sx={{ color: theme.palette.info.main }}>
                Load: <strong>{frameLoadTime}ms</strong>
              </Typography>
            )}
            {isEnabled && propertyValues.map((query, index) => {
              const key = enabledProperties[index];

              // Skip if loading or error
              if (query.isLoading || query.isError || !query.data?.value) {
                return null;
              }

              const value = query.data.value;

              // Format value using utility function (no particleId for global props)
              const displayValue = formatPropertyValue(value);

              return (
                <Typography
                  key={key}
                  variant="caption"
                  sx={{
                    color: theme.palette.text.secondary,
                    fontFamily: "monospace",
                    fontSize: "0.75rem",
                  }}
                >
                  {key}: {displayValue}
                </Typography>
              );
            })}
          </Box>
        </Box>
      </Paper>
    </Rnd>
  );
}
