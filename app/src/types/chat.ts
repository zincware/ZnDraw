export interface ChatMessage {
  id: string;
  roomId: string;
  author: {
    id: string;
  };
  content: string;
  createdAt: number;
  updatedAt: number;
  isEdited: boolean;
}

export interface ChatMessagesResponse {
  messages: ChatMessage[];
  metadata: {
    hasMore: boolean;
    totalCount: number;
    oldestTimestamp: number | null;
    newestTimestamp: number | null;
  };
}
