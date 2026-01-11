/**
 * Settings hooks for managing session settings per room.
 *
 * Settings are per-session, per-room (identified via session ID in URL path).
 * Each browser window/tab has its own settings.
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
	const sessionId = useAppStore((state) => state.sessionId);

	return useQuery({
		queryKey: ["settings", roomId, sessionId],
		queryFn: async () => {
			if (!sessionId) throw new Error("No session ID");
			return await getSettings(roomId, sessionId);
		},
		staleTime: Infinity, // Settings don't change often, rely on socket invalidation
		enabled: !!roomId && !!sessionId,
	});
};

/**
 * Hook to update settings categories for the current session.
 *
 * Accepts partial updates - only provided categories are sent.
 * Automatically updates the query cache on success.
 */
export const useUpdateSettings = () => {
	const queryClient = useQueryClient();
	const sessionId = useAppStore((state) => state.sessionId);

	return useMutation({
		mutationFn: async (variables: {
			roomId: string;
			data: Record<string, any>;
		}) => {
			if (!sessionId) throw new Error("No session ID");
			const { roomId, data } = variables;
			return await updateSettings(roomId, sessionId, data);
		},
		onSuccess: (_, variables) => {
			const { roomId, data: submittedData } = variables;

			// Optimistically update the cache with the submitted data
			const queryKey = ["settings", roomId, sessionId];
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
