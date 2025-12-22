/**
 * Settings hooks for managing user settings per room.
 *
 * Settings are always per-user, per-room (authenticated via JWT).
 * Unlike extensions, settings are never public and don't use workers/job queues.
 */

import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import {
	getSettings,
	updateSettings,
	type SettingsResponse,
} from "../myapi/client";
import { useAppStore } from "../store";

export type { SettingsResponse };

/**
 * Hook to fetch all settings (schema + data) for a room.
 * Returns both the JSON schema and current data for all settings categories.
 */
export const useSettings = (roomId: string) => {
	const userName = useAppStore((state) => state.userName);

	return useQuery({
		queryKey: ["settings", roomId, userName],
		queryFn: async () => {
			return await getSettings(roomId);
		},
		staleTime: Infinity, // Settings don't change often, rely on socket invalidation
		enabled: !!roomId && !!userName,
	});
};

/**
 * Hook to update settings categories for the current user.
 *
 * Accepts partial updates - only provided categories are sent.
 * Automatically updates the query cache on success.
 */
export const useUpdateSettings = () => {
	const queryClient = useQueryClient();
	const userName = useAppStore((state) => state.userName);

	return useMutation({
		mutationFn: async (variables: {
			roomId: string;
			data: Record<string, any>;
		}) => {
			const { roomId, data } = variables;
			return await updateSettings(roomId, data);
		},
		onSuccess: (_, variables) => {
			const { roomId, data: submittedData } = variables;

			// Optimistically update the cache with the submitted data
			const queryKey = ["settings", roomId, userName];
			queryClient.setQueryData(
				queryKey,
				(old: SettingsResponse | undefined) => {
					if (!old) return old;
					return {
						...old,
						data: {
							...old.data,
							...submittedData,
						},
					};
				},
			);
		},
		onError: (error) => {
			console.error("Error updating settings:", error);
		},
	});
};
