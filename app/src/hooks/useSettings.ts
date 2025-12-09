/**
 * Settings hooks for managing user settings per room.
 *
 * Settings are always per-user, per-room (authenticated via JWT).
 * Unlike extensions, settings are never public and don't use workers/job queues.
 */

import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import {
	getSettingsSchemas,
	getSetting,
	updateSetting,
	type SettingsSchema,
} from "../myapi/client";
import { useAppStore } from "../store";

export type { SettingsSchema };

/**
 * Hook to fetch settings schemas for a room.
 * Returns the JSON schemas for all settings categories.
 */
export const useSettingsSchemas = (roomId: string) => {
	return useQuery({
		queryKey: ["settingsSchemas", roomId],
		queryFn: async () => {
			return await getSettingsSchemas(roomId);
		},
		enabled: !!roomId,
	});
};

/**
 * Hook to fetch a specific settings category for the current user.
 *
 * Backend always returns valid data (defaults if not stored), so no need
 * for initialData on the frontend. This prevents race conditions where
 * frontend defaults could overwrite stored settings on page reload.
 *
 * @param roomId - Room identifier
 * @param category - Settings category (camera, studio_lighting, etc.)
 */
export const useSettingData = (roomId: string, category: string) => {
	const userName = useAppStore((state) => state.userName);

	return useQuery({
		queryKey: ["settingData", roomId, userName, category],
		queryFn: async () => {
			const result = await getSetting(roomId, category);
			return result.data;
		},
		staleTime: Infinity, // Settings don't change often, rely on socket invalidation
		enabled: !!roomId && !!userName && !!category,
	});
};

/**
 * Hook to update a settings category for the current user.
 *
 * Automatically updates the query cache on success.
 */
export const useUpdateSetting = () => {
	const queryClient = useQueryClient();
	const userName = useAppStore((state) => state.userName);

	return useMutation({
		mutationFn: async (variables: {
			roomId: string;
			category: string;
			data: any;
		}) => {
			const { roomId, category, data } = variables;
			return await updateSetting(roomId, category, data);
		},
		onSuccess: (_, variables) => {
			const { roomId, category, data: submittedData } = variables;

			// Optimistically update the cache with the submitted data
			const queryKey = ["settingData", roomId, userName, category];
			queryClient.setQueryData(queryKey, submittedData);
		},
		onError: (error) => {
			console.error("Error updating setting:", error);
		},
	});
};
