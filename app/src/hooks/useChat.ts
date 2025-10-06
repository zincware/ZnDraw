import {
  useInfiniteQuery,
  useMutation,
  useQueryClient,
} from "@tanstack/react-query";
import { socket } from "../socket";
import type { ChatMessage, ChatMessagesResponse } from "../types/chat";

// Fetch chat messages with infinite scroll support
// Note: Using fetch directly here because socket.io handles chat updates
// and the response types match the types/chat.ts definitions
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
      return new Promise<ChatMessage>((resolve, reject) => {
        socket.emit("chat:message:create", { content }, (response: any) => {
          if (response.success) resolve(response.message);
          else reject(new Error(response.error));
        });
      });
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
      return new Promise<ChatMessage>((resolve, reject) => {
        socket.emit(
          "chat:message:edit",
          { messageId, content },
          (response: any) => {
            if (response.success) resolve(response.message);
            else reject(new Error(response.error));
          },
        );
      });
    },
    // Note: We don't invalidate queries here - socket events will update the cache
  });
};
