import { useMutation, useQuery, useQueryClient } from "@tanstack/react-query";
import {
	type FigureData,
	type FigureListResponse,
	createFigure,
	deleteFigure,
	getFigure,
	listFigures,
} from "../myapi/client";
import { useAppStore } from "../store";

// Define query keys for caching and invalidation
const figuresKeys = {
	all: (roomId: string) => ["figures", roomId] as const,
	list: (roomId: string) => [...figuresKeys.all(roomId), "list"] as const,
	detail: (roomId: string, key: string) =>
		[...figuresKeys.all(roomId), "detail", key] as const,
};

// --- Custom Hooks ---

/**
 * Fetches the list of all figure keys for the current room.
 */
export const useFigureList = () => {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);

	return useQuery({
		queryKey: figuresKeys.list(roomId!),
		queryFn: () => listFigures(roomId!),
		enabled: !!roomId, // Only run the query if a roomId is set
	});
};

/**
 * Fetches the data for a single figure.
 */
export const useFigure = (
	key: string | null,
	options?: { enabled?: boolean },
) => {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);

	return useQuery({
		queryKey: figuresKeys.detail(roomId!, key!),
		queryFn: () => getFigure(roomId!, key!),
		enabled:
			(options?.enabled !== undefined ? options.enabled : true) &&
			!!roomId &&
			!!key, // Only run if both roomId and key are set, and options.enabled is true
	});
};

/**
 * Provides a mutation function for creating or updating a figure.
 * Invalidates the figure list on success to trigger a refetch.
 */
export const useCreateFigure = () => {
	const queryClient = useQueryClient();
	const roomId = useAppStore((state) => state.roomId);

	return useMutation({
		mutationFn: (variables: { key: string; figure: FigureData }) =>
			createFigure(roomId!, variables.key, variables.figure),
		onMutate: async (variables) => {
			await queryClient.cancelQueries({
				queryKey: figuresKeys.list(roomId!),
			});
			const previous = queryClient.getQueryData<FigureListResponse>(
				figuresKeys.list(roomId!),
			);
			if (previous && !previous.items.includes(variables.key)) {
				queryClient.setQueryData<FigureListResponse>(
					figuresKeys.list(roomId!),
					{ ...previous, items: [...previous.items, variables.key] },
				);
			}
			return { previous };
		},
		onError: (_err, _variables, context) => {
			if (context?.previous) {
				queryClient.setQueryData(figuresKeys.list(roomId!), context.previous);
			}
		},
		onSettled: () => {
			queryClient.invalidateQueries({ queryKey: figuresKeys.list(roomId!) });
		},
	});
};

/**
 * Provides a mutation function for deleting a figure.
 * Invalidates the list and removes the specific figure from the cache.
 */
export const useDeleteFigure = () => {
	const queryClient = useQueryClient();
	const roomId = useAppStore((state) => state.roomId);

	return useMutation({
		mutationFn: (key: string) => deleteFigure(roomId!, key),
		onMutate: async (key) => {
			await queryClient.cancelQueries({
				queryKey: figuresKeys.list(roomId!),
			});
			const previous = queryClient.getQueryData<FigureListResponse>(
				figuresKeys.list(roomId!),
			);
			if (previous) {
				queryClient.setQueryData<FigureListResponse>(
					figuresKeys.list(roomId!),
					{ ...previous, items: previous.items.filter((k) => k !== key) },
				);
			}
			queryClient.removeQueries({
				queryKey: figuresKeys.detail(roomId!, key),
			});
			return { previous };
		},
		onError: (_err, _key, context) => {
			if (context?.previous) {
				queryClient.setQueryData(figuresKeys.list(roomId!), context.previous);
			}
		},
		onSettled: () => {
			queryClient.invalidateQueries({ queryKey: figuresKeys.list(roomId!) });
		},
	});
};
