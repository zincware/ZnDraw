import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import {
	type GeometryListResponse,
	type GeometryResponse,
	createGeometry,
	deleteGeometry,
	getGeometry,
	listGeometries,
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
		onMutate: async (variables) => {
			await queryClient.cancelQueries({
				queryKey: ["geometries", variables.roomId, "list"],
			});
			const previous = queryClient.getQueryData<GeometryListResponse>([
				"geometries",
				variables.roomId,
				"list",
			]);
			if (previous) {
				queryClient.setQueryData<GeometryListResponse>(
					["geometries", variables.roomId, "list"],
					{
						...previous,
						items: {
							...previous.items,
							[variables.key]: {
								type: variables.geometryType,
								data: variables.geometryData,
								selection: [],
							},
						},
					},
				);
			}
			return { previous };
		},
		onError: (_err, variables, context) => {
			if (context?.previous) {
				queryClient.setQueryData(
					["geometries", variables.roomId, "list"],
					context.previous,
				);
			}
		},
		onSettled: (_data, _err, variables) => {
			queryClient.invalidateQueries({
				queryKey: ["geometries", variables.roomId, "list"],
			});
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
		onMutate: async (variables) => {
			await queryClient.cancelQueries({
				queryKey: ["geometries", variables.roomId, "list"],
			});
			const previous = queryClient.getQueryData<GeometryListResponse>([
				"geometries",
				variables.roomId,
				"list",
			]);
			if (previous) {
				const { [variables.key]: _, ...rest } = previous.items;
				queryClient.setQueryData<GeometryListResponse>(
					["geometries", variables.roomId, "list"],
					{ ...previous, items: rest },
				);
			}
			queryClient.removeQueries({
				queryKey: ["geometries", variables.roomId, "detail", variables.key],
			});
			return { previous };
		},
		onError: (_err, variables, context) => {
			if (context?.previous) {
				queryClient.setQueryData(
					["geometries", variables.roomId, "list"],
					context.previous,
				);
			}
		},
		onSettled: (_data, _err, variables) => {
			queryClient.invalidateQueries({
				queryKey: ["geometries", variables.roomId, "list"],
			});
		},
	});
};
