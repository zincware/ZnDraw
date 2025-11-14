import { useState, useEffect } from "react";
import { Box, Typography, Paper, useTheme } from "@mui/material";
import { useAppStore } from "../../store";
import { usePropertyInspectorSettings } from "../../hooks/usePropertyInspectorSettings";
import { formatPropertyValue } from "../../utils/propertyFormatting";

/**
 * HoverInfoBox - Displays contextual information at the mouse position.
 * Shows:
 * - Particle index when hovering over a particle
 * - Curve length when in drawing mode with an active curve
 */
export default function HoverInfoBox() {
  const theme = useTheme();
  const {
    showInfoBoxes,
    hoveredGeometryInstance,
    curveLength,
    mode,
  } = useAppStore();

  // Extract particle ID from hovered geometry instance
  const hoveredParticleId = hoveredGeometryInstance?.geometryKey === "particles"
    ? hoveredGeometryInstance.instanceId
    : null;

  const [mousePosition, setMousePosition] = useState({ x: 0, y: 0 });

  // Performance optimization: only fetch per-particle properties when boxes are visible and particle is hovered
  const shouldFetchProperties = showInfoBoxes && hoveredParticleId !== null;

  const {
    enabledProperties,
    propertyValues,
    isEnabled,
    categories,
  } = usePropertyInspectorSettings({
    category: "per-particle",
    enabled: shouldFetchProperties,
  });

  useEffect(() => {
    const handleMouseMove = (event: MouseEvent) => {
      setMousePosition({ x: event.clientX, y: event.clientY });
    };

    window.addEventListener("mousemove", handleMouseMove);
    return () => {
      window.removeEventListener("mousemove", handleMouseMove);
    };
  }, []);

  if (!showInfoBoxes) return null;

  // Determine what to show
  const showParticleInfo = hoveredParticleId != null; // Loose equality checks both null and undefined
  const showCurveInfo = mode === 'drawing' && curveLength > 0;
  const hasContent = showParticleInfo || showCurveInfo;

  if (!hasContent) return null;

  return (
    <Box
      sx={{
        position: "fixed",
        left: mousePosition.x + 20,
        top: mousePosition.y + 20,
        pointerEvents: "none",
        userSelect: "none",
        zIndex: 10000,
      }}
    >
      <Paper
        elevation={3}
        sx={{
          padding: 1.5,
          backgroundColor: theme.palette.mode === "dark" 
            ? "rgba(0, 0, 0, 0.85)" 
            : "rgba(255, 255, 255, 0.95)",
          backdropFilter: "blur(10px)",
          borderRadius: 1,
          minWidth: "120px",
          border: `1px solid ${theme.palette.mode === "dark" 
            ? "rgba(255, 255, 255, 0.15)" 
            : "rgba(0, 0, 0, 0.15)"}`,
        }}
      >
        <Box sx={{ display: "flex", flexDirection: "column", gap: 0.5 }}>
          {showParticleInfo && (
            <>
              <Typography
                variant="body2"
                sx={{ color: theme.palette.primary.main, fontWeight: "bold" }}
              >
                Particle: {hoveredParticleId}
              </Typography>
              {isEnabled && hoveredParticleId != null && enabledProperties.length > 0 && propertyValues.map((query, index) => {
                const key = enabledProperties[index];

                // Skip if loading or error
                if (query.isLoading || query.isError || !query.data?.value) {
                  return null;
                }

                const value = query.data.value;

                // Find property metadata for proper formatting
                const propInfo = categories?.perParticle.find(p => p.key === key);

                // Skip if property info not found (shouldn't happen but defensive)
                if (!propInfo) {
                  return null;
                }

                // Format value using utility function
                const displayValue = formatPropertyValue(value, hoveredParticleId, propInfo);

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
            </>
          )}
          {showCurveInfo && (
            <Typography
              variant="body2"
              sx={{ color: theme.palette.warning.main, fontWeight: "bold" }}
            >
              Length: {curveLength.toFixed(2)}
            </Typography>
          )}
        </Box>
      </Paper>
    </Box>
  );
}
