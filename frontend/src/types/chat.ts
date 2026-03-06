export interface ChatMessage {
	id: number;
	room_id: string;
	user_id: string;
	content: string;
	created_at: string; // ISO datetime
	updated_at: string | null;
	email: string | null;
}

/** Socket event payload for `message_new`. Matches backend `MessageNew`. */
export interface MessageNewEvent {
	id: number;
	room_id: string;
	user_id: string;
	content: string;
	created_at: string;
	updated_at: string | null;
	email: string | null;
}

/** Socket event payload for `message_edited`. Matches backend `MessageEdited`. */
export interface MessageEditedEvent {
	id: number;
	room_id: string;
	content: string;
	updated_at: string;
}

export interface ChatMessagesResponse {
	items: ChatMessage[];
	metadata: {
		has_more: boolean;
		total_count: number;
		oldest_timestamp: number | null;
		newest_timestamp: number | null;
	};
}
