/**
 * Hook to access PropertyInspector settings and fetch property values.
 * Integrates with the settings stored on the backend.
 *
 * Performance optimized:
 * - Category filtering to avoid fetching unused properties
 * - Conditional fetching based on enabled state
 * - Memoized property categorization
 */

import { useMemo } from "react";
import { useExtensionData } from "./useSchemas";
import { usePropertyValues, useAvailableProperties } from "./usePropertyInspector";
import { useAppStore } from "../store";

type PropertyCategory = "per-particle" | "global" | "all";

interface UsePropertyInspectorSettingsOptions {
  /** Category filter to only fetch relevant properties */
  category?: PropertyCategory;
  /** Whether to enable property fetching (e.g., based on showInfoBoxes) */
  enabled?: boolean;
}

/**
 * Hook to get the enabled properties from PropertyInspector settings.
 * Returns the list of enabled property keys and their values, optionally filtered by category.
 *
 * @param options - Configuration options for filtering and enabling
 * @returns Enabled properties, their values, and metadata
 */
export const usePropertyInspectorSettings = (
  options: UsePropertyInspectorSettingsOptions = {}
) => {
  const { category = "all", enabled = true } = options;
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  const userName = useAppStore((state) => state.userName);
  const currentFrame = useAppStore((state) => state.currentFrame);
  const particleCount = useAppStore((state) => state.particleCount);

  // Fetch property_inspector settings from backend
  const { data: settingsData } = useExtensionData(
    roomId || "",
    userName || "",
    "settings",
    "property_inspector"
  );

  // Extract enabled_properties array from settings
  const enabledProperties: string[] = settingsData?.enabled_properties || [];

  // Fetch property categorization metadata (only if enabled and we need filtering)
  const shouldFetchMetadata = enabled && category !== "all" && enabledProperties.length > 0;

  const { data: categories } = useAvailableProperties(
    roomId || undefined,
    currentFrame,
    particleCount,
    shouldFetchMetadata
  );

  // Memoize filtered properties based on category
  const filteredProperties = useMemo(() => {
    if (category === "all" || !categories) {
      return enabledProperties;
    }

    // Create lookup sets for efficient filtering
    const perParticleKeys = new Set(categories.perParticle.map(p => p.key));
    const globalKeys = new Set(categories.global.map(p => p.key));

    return enabledProperties.filter(key => {
      if (category === "per-particle") return perParticleKeys.has(key);
      if (category === "global") return globalKeys.has(key);
      return true;
    });
  }, [enabledProperties, categories, category]);

  // Fetch values for filtered properties (only if enabled)
  const propertyValues = usePropertyValues(
    roomId || undefined,
    currentFrame,
    filteredProperties,
    enabled && filteredProperties.length > 0
  );

  return {
    enabledProperties: filteredProperties,
    propertyValues,
    isEnabled: filteredProperties.length > 0,
    categories,
  };
};
