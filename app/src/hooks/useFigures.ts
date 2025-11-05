import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import {
  listFigures,
  getFigure,
  createFigure,
  deleteFigure,
  FigureData,
} from "../myapi/client";
import { useAppStore } from "../store";
import { use } from "react";

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
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);

  return useMutation({
    mutationFn: (variables: { key: string; figure: FigureData }) =>
      createFigure(roomId!, variables.key, variables.figure),
    onSuccess: (_, variables) => {
      // When a figure is created, the list is now out of date. Invalidate it.
      queryClient.invalidateQueries({ queryKey: figuresKeys.list(roomId!) });
      // We can also pre-populate the cache for the new figure if desired
      // queryClient.setQueryData(figuresKeys.detail(roomId!, variables.key), { key: variables.key, figure: variables.figure });
    },
  });
};

/**
 * Provides a mutation function for deleting a figure.
 * Invalidates the list and removes the specific figure from the cache.
 */
export const useDeleteFigure = () => {
  const queryClient = useQueryClient();
  // Use individual selectors to prevent unnecessary re-renders
  const roomId = useAppStore((state) => state.roomId);
  //   const setSelectedFigureKey = useRoomStore((state) => state.setSelectedFigureKey); // TODO

  return useMutation({
    mutationFn: (key: string) => deleteFigure(roomId!, key),
    onSuccess: (_, key) => {
      // Invalidate the list to refetch it without the deleted key
      queryClient.invalidateQueries({ queryKey: figuresKeys.list(roomId!) });
      // Remove the stale data for the deleted figure from the cache immediately
      queryClient.removeQueries({ queryKey: figuresKeys.detail(roomId!, key) });

      // If the deleted figure was the selected one, unselect it
      //   if (useRoomStore.getState().selectedFigureKey === key) {
      //     setSelectedFigureKey(null);
      //   }
    },
  });
};
