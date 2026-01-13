import {
	useInfiniteQuery,
	useMutation,
	useQueryClient,
} from "@tanstack/react-query";
import type { ChatMessage, ChatMessagesResponse } from "../types/chat";
import { createChatMessage, editChatMessage } from "../myapi/client";

// Fetch chat messages with infinite scroll support
export const useChatMessages = (room: string, limit: number = 30) => {
	return useInfiniteQuery<ChatMessagesResponse>({
		queryKey: ["chat", room],
		queryFn: async ({ pageParam }) => {
			const params = new URLSearchParams();
			params.append("limit", limit.toString());
			if (pageParam) {
				params.append("before", pageParam.toString());
			}
			const response = await fetch(
				`/api/rooms/${room}/chat/messages?${params}`,
			);
			if (!response.ok) throw new Error("Failed to fetch messages");
			return response.json();
		},
		getNextPageParam: (lastPage) =>
			lastPage.metadata.hasMore ? lastPage.metadata.oldestTimestamp : undefined,
		initialPageParam: undefined,
		enabled: !!room,
		staleTime: 1000 * 60, // 1 minute
	});
};

// Send new message
export const useSendMessage = (room: string) => {
	return useMutation({
		mutationFn: async (content: string) => {
			const response = await createChatMessage(room, content);
			return response.message;
		},
		// Note: We don't invalidate queries here - socket events will update the cache
	});
};

// Edit message
export const useEditMessage = (room: string) => {
	return useMutation({
		mutationFn: async ({
			messageId,
			content,
		}: {
			messageId: string;
			content: string;
		}) => {
			const response = await editChatMessage(room, messageId, content);
			return response.message;
		},
		// Note: We don't invalidate queries here - socket events will update the cache
	});
};
