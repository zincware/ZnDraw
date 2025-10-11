/**
 * React Query hooks for Property Inspector functionality.
 * Handles fetching and categorizing frame properties.
 */

import { useQuery, useQueries, keepPreviousData } from "@tanstack/react-query";
import {
  getFrameMetadata,
  getFrames,
  categorizeProperties,
} from "../myapi/client";

/**
 * Hook to get available properties with categorization.
 * Fetches frame metadata and categorizes properties into per-particle and global.
 *
 * @param roomId - Room identifier
 * @param frameId - Frame index to fetch metadata for
 * @param particleCount - Number of particles in the frame
 * @param enabled - Whether to enable the query (default: true)
 * @returns React Query result with categorized properties
 */
export const useAvailableProperties = (
  roomId: string | undefined,
  frameId: number,
  particleCount: number,
  enabled: boolean = true
) => {
  return useQuery({
    queryKey: ["properties", "available", roomId, frameId],
    queryFn: async () => {
      if (!roomId) throw new Error("Room ID is required");
      const metadata = await getFrameMetadata(roomId, frameId);
      return categorizeProperties(metadata, particleCount);
    },
    enabled: enabled && !!roomId && particleCount > 0,
    staleTime: 5 * 60 * 1000, // 5 minutes - metadata rarely changes
    placeholderData: keepPreviousData, // Prevent flickering between frames
  });
};

/**
 * Hook to fetch selected property values.
 * Uses parallel queries for optimal performance.
 *
 * @param roomId - Room identifier
 * @param frameId - Frame index to fetch values for
 * @param propertyKeys - Array of property keys to fetch
 * @param enabled - Whether to enable the queries (default: true)
 * @returns Array of React Query results for each property
 */
export const usePropertyValues = (
  roomId: string | undefined,
  frameId: number,
  propertyKeys: string[],
  enabled: boolean = true
) => {
  return useQueries({
    queries: propertyKeys.map((key) => ({
      queryKey: ["property", "value", roomId, frameId, key],
      queryFn: async () => {
        if (!roomId) throw new Error("Room ID is required");
        const data = await getFrames(roomId, frameId, [key]);
        return {
          key,
          value: data?.[key],
        };
      },
      enabled: enabled && !!roomId && propertyKeys.length > 0,
      staleTime: 1000, // Cache for 1 second during rapid frame changes
      placeholderData: keepPreviousData, // Prevent flickering between frames
    })),
  });
};
