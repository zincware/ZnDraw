import {
	useInfiniteQuery,
	useMutation,
	useQueryClient,
} from "@tanstack/react-query";
import {
	createChatMessage,
	editChatMessage,
	getChatMessages,
} from "../myapi/client";
import type { ChatMessage, ChatMessagesResponse } from "../types/chat";

// Fetch chat messages with infinite scroll support
export const useChatMessages = (room: string, limit = 30) => {
	return useInfiniteQuery<ChatMessagesResponse>({
		queryKey: ["chat", room],
		queryFn: async ({ pageParam }) => {
			return getChatMessages(room, limit, pageParam as number | undefined);
		},
		getNextPageParam: (lastPage) =>
			lastPage.metadata.has_more
				? lastPage.metadata.oldest_timestamp
				: undefined,
		initialPageParam: undefined,
		enabled: !!room,
		staleTime: 1000 * 60, // 1 minute
	});
};

// Send new message
export const useSendMessage = (room: string) => {
	return useMutation({
		mutationFn: async (content: string) => {
			return await createChatMessage(room, content);
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
			messageId: number;
			content: string;
		}) => {
			return await editChatMessage(room, messageId, content);
		},
		// Note: We don't invalidate queries here - socket events will update the cache
	});
};
