/**
 * PropertyInspectorRenderer - Custom JSONForms control for property inspection.
 *
 * Renders when schema has x-custom-type="property-inspector"
 * Displays property selector + table for per-particle and global properties.
 */

import { useMemo, useCallback } from "react";
import { withJsonFormsControlProps } from "@jsonforms/react";
import { rankWith, schemaMatches, and, uiTypeIs, ControlProps } from "@jsonforms/core";
import {
  Box,
  Typography,
  Divider,
  CircularProgress,
  Paper,
  Alert,
} from "@mui/material";
import { useAppStore } from "../../store";
import { useAvailableProperties } from "../../hooks/usePropertyInspector";
import PropertySelector from "./components/PropertySelector";

/**
 * PropertyInspectorRenderer Component
 *
 * ControlProps provided by withJsonFormsControlProps:
 * - data: string[] - Array of selected property keys
 * - handleChange: (path, value) => void - Update form data
 * - path: string - JSON path to this field
 * - label: string - Field label from schema
 * - schema: any - JSON schema for this field
 * - required: boolean - Whether field is required
 * - errors: string - Validation errors
 */
const PropertyInspectorRenderer = ({
  data,
  handleChange,
  path,
  label,
  errors,
}: ControlProps) => {
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const particleCount = useAppStore((state) => state.particleCount);
  const currentFrame = useAppStore((state) => state.currentFrame);

  // Fetch available properties with categorization
  const {
    data: categories,
    isLoading,
    isError,
  } = useAvailableProperties(roomId || undefined, currentFrame, particleCount);

  // data is the array of selected property keys (from enabled_properties field)
  const selectedKeys: string[] = useMemo(() => data || [], [data]);

  // Memoized handler to update JSONForms data
  const handlePropertyToggle = useCallback((key: string) => {
    const keys = new Set(selectedKeys);
    if (keys.has(key)) {
      keys.delete(key);
    } else {
      keys.add(key);
    }
    handleChange(path, Array.from(keys));
  }, [selectedKeys, handleChange, path]);

  const handleClearAll = useCallback(() => {
    handleChange(path, []);
  }, [handleChange, path]);

  if (!roomId) {
    return (
      <Paper variant="outlined" sx={{ mb: 2, p: 2 }}>
        <Alert severity="info">
          Property Inspector requires an active room connection
        </Alert>
      </Paper>
    );
  }

  if (!categories && isLoading) {
    return (
      <Paper variant="outlined" sx={{ mb: 2, p: 2 }}>
        <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", gap: 2 }}>
          <CircularProgress size={24} />
          <Typography>Loading properties...</Typography>
        </Box>
      </Paper>
    );
  }

  if (isError) {
    return (
      <Paper variant="outlined" sx={{ mb: 2, p: 2 }}>
        <Alert severity="error">Failed to load property metadata</Alert>
      </Paper>
    );
  }

  return (
    <Paper variant="outlined" sx={{ mb: 2, overflow: "hidden" }}>
      {/* Header */}
      <Box sx={{ p: 2, bgcolor: "action.hover" }}>
        <Typography variant="subtitle2" fontWeight="bold">
          {label || "Property Inspector"}
        </Typography>
        <Typography variant="caption" color="text.secondary">
          Select properties to display in InfoBoxes (press 'i' to toggle)
        </Typography>
      </Box>

      <Divider />

      {/* Property Selector - Collapsible multi-select */}
      <Box sx={{ maxHeight: 400, overflow: "auto" }}>
        <PropertySelector
          perParticleProps={categories?.perParticle || []}
          globalProps={categories?.global || []}
          selectedKeys={selectedKeys}
          onToggle={handlePropertyToggle}
          onClearAll={handleClearAll}
        />
      </Box>

      {errors && (
        <Alert severity="error" sx={{ m: 1 }}>
          {errors}
        </Alert>
      )}
    </Paper>
  );
};

/**
 * Tester function - determines when to use this renderer.
 * Priority 10 to override default array renderers.
 */
export const propertyInspectorTester = rankWith(
  10,
  and(
    schemaMatches((schema) => (schema as any)["x-custom-type"] === "property-inspector"),
    uiTypeIs("Control")
  )
);

// Export wrapped with JSONForms HOC
export default withJsonFormsControlProps(PropertyInspectorRenderer);
