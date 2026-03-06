/**
 * React Query hooks for Property Inspector functionality.
 * Handles fetching and categorizing frame properties.
 */

import {
	keepPreviousData,
	useQueries,
	useQuery,
	useQueryClient,
} from "@tanstack/react-query";
import { categorizeProperties, getFrameMetadata } from "../myapi/client";
import { getFrameBatched } from "./useFrameBatch";

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
	enabled = true,
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
	enabled = true,
) => {
	const queryClient = useQueryClient();
	return useQueries({
		queries: propertyKeys.map((key) => ({
			queryKey: ["frame", roomId, frameId, key],
			queryFn: () => getFrameBatched(roomId!, frameId, key, queryClient),
			select: (data: Awaited<ReturnType<typeof getFrameBatched>>) => ({
				key,
				value: data?.[key],
			}),
			enabled: enabled && !!roomId && propertyKeys.length > 0,
			placeholderData: keepPreviousData,
		})),
	});
};
