import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import {
	listGeometries,
	getGeometry,
	createGeometry,
	deleteGeometry,
	GeometryListResponse,
	GeometryResponse,
} from "../myapi/client";

/**
 * Hook to fetch the list of geometry keys for a room
 */
export const useGeometriesList = (roomId: string | null) => {
	return useQuery<GeometryListResponse>({
		queryKey: ["geometries", roomId, "list"],
		queryFn: () => listGeometries(roomId!),
		enabled: !!roomId,
		staleTime: 30000, // 30 seconds
	});
};

/**
 * Hook to fetch a specific geometry by key
 */
export const useGeometry = (roomId: string | null, key: string | null) => {
	return useQuery<GeometryResponse>({
		queryKey: ["geometries", roomId, "detail", key],
		queryFn: () => getGeometry(roomId!, key!),
		enabled: !!roomId && !!key,
		staleTime: 30000,
	});
};

/**
 * Hook to create a new geometry
 */
export const useCreateGeometry = () => {
	const queryClient = useQueryClient();

	return useMutation({
		mutationFn: ({
			roomId,
			key,
			geometryType,
			geometryData,
		}: {
			roomId: string;
			key: string;
			geometryType: string;
			geometryData: Record<string, any>;
		}) => createGeometry(roomId, key, geometryType, geometryData),
		onSuccess: (_, variables) => {
			// Invalidate the geometries list to refetch
			queryClient.invalidateQueries({
				queryKey: ["geometries", variables.roomId, "list"],
			});
			// Also invalidate the specific geometry detail
			queryClient.invalidateQueries({
				queryKey: ["geometries", variables.roomId, "detail", variables.key],
			});
		},
	});
};

/**
 * Hook to delete a geometry
 */
export const useDeleteGeometry = () => {
	const queryClient = useQueryClient();

	return useMutation({
		mutationFn: ({ roomId, key }: { roomId: string; key: string }) =>
			deleteGeometry(roomId, key),
		onSuccess: (_, variables) => {
			// Remove the specific geometry from the cache
			queryClient.removeQueries({
				queryKey: ["geometries", variables.roomId, "detail", variables.key],
			});
			// Invalidate the list to refetch
			queryClient.invalidateQueries({
				queryKey: ["geometries", variables.roomId, "list"],
			});
		},
	});
};
