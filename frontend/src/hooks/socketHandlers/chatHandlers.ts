import { useAppStore } from "../../store";
import type { MessageEditedEvent, MessageNewEvent } from "../../types/chat";
import type { HandlerContext } from "./types";

/** Socket event payload for the `typing` event. */
export interface TypingEvent {
	user_id: string;
	email: string;
	is_typing: boolean;
}

interface ChatHandlers {
	onChatMessageNew: (data: MessageNewEvent) => void;
	onChatMessageUpdated: (data: MessageEditedEvent) => void;
	onTyping: (data: TypingEvent) => void;
}

export interface ChatHandlersResult {
	handlers: ChatHandlers;
	cleanup: () => void;
}

/**
 * Creates chat handler functions with cleanup for typing timeouts.
 *
 * The typingTimeouts Map is owned by this module (mutable state, not in
 * HandlerContext). The returned `cleanup` function clears all pending
 * timeouts and must be called in the useEffect cleanup.
 */
export function createChatHandlers(ctx: HandlerContext): ChatHandlersResult {
	const typingTimeouts = new Map<string, ReturnType<typeof setTimeout>>();

	function onChatMessageNew(data: MessageNewEvent) {
		ctx.queryClient.setQueryData(["chat", ctx.roomId], (oldData: any) => {
			if (!oldData) return oldData;
			const newPages = [...oldData.pages];
			const lastPageIndex = newPages.length - 1;

			if (lastPageIndex >= 0) {
				newPages[lastPageIndex] = {
					...newPages[lastPageIndex],
					items: [...newPages[lastPageIndex].items, data],
					metadata: {
						...newPages[lastPageIndex].metadata,
						total_count: newPages[lastPageIndex].metadata.total_count + 1,
					},
				};
			}
			return { ...oldData, pages: newPages };
		});

		// Increment unread count if chat panel is not active in any bar
		const state = useAppStore.getState();
		const chatIsActive =
			state.activeLeft === "chat" ||
			state.activeRight === "chat" ||
			state.activeBottom === "chat";
		if (!chatIsActive) {
			state.incrementChatUnread();
		}
	}

	function onChatMessageUpdated(data: MessageEditedEvent) {
		ctx.queryClient.setQueryData(["chat", ctx.roomId], (oldData: any) => {
			if (!oldData) return oldData;
			const newPages = oldData.pages.map((page: any) => ({
				...page,
				items: page.items.map((msg: any) =>
					msg.id === data.id ? { ...msg, ...data } : msg,
				),
			}));
			return { ...oldData, pages: newPages };
		});
	}

	function onTyping(data: TypingEvent) {
		const { addTypingUser, removeTypingUser } = useAppStore.getState();
		const email = data.email;

		// Clear previous timeout for this user
		const prev = typingTimeouts.get(email);
		if (prev) clearTimeout(prev);

		if (data.is_typing) {
			addTypingUser(email);
			// Auto-remove after 5s in case typing_stop is missed
			typingTimeouts.set(
				email,
				setTimeout(() => {
					removeTypingUser(email);
					typingTimeouts.delete(email);
				}, 5000),
			);
		} else {
			removeTypingUser(email);
			typingTimeouts.delete(email);
		}
	}

	function cleanup() {
		for (const timeout of typingTimeouts.values()) clearTimeout(timeout);
		typingTimeouts.clear();
	}

	return {
		handlers: { onChatMessageNew, onChatMessageUpdated, onTyping },
		cleanup,
	};
}
