import type { HandlerContext } from "./types";

/** Socket event payload for the `invalidate` event. */
export interface InvalidateEvent {
	roomId: string;
	userName: string;
	category: string;
	extension: string;
	sessionId?: string;
}

/** Socket event payload for the `schema_invalidate` event. */
export interface SchemaInvalidateEvent {
	category: string;
}

/**
 * Creates scene invalidation handler functions.
 *
 * Note: onSchemaInvalidate uses `ctx.appStoreRoomId` (not `ctx.roomId`)
 * to match the existing behavior for schema query keys.
 */
export function createSceneInvalidationHandlers(ctx: HandlerContext) {
	function onInvalidate(data: InvalidateEvent) {
		const { roomId, userName, category, extension, sessionId } = data;
		// Invalidate extension data queries (for modifiers, analysis, selections)
		ctx.queryClient.invalidateQueries({
			queryKey: ["extensionData", roomId, userName, category, extension],
		});
		// Invalidate settings queries (per-session)
		if (category === "settings" && sessionId) {
			ctx.queryClient.invalidateQueries({
				queryKey: ["settings", roomId, sessionId],
			});
		}
	}

	function onSchemaInvalidate(data: SchemaInvalidateEvent) {
		const { category } = data;
		ctx.queryClient.invalidateQueries({
			queryKey: ["schemas", ctx.appStoreRoomId, category],
		});
	}

	return { onInvalidate, onSchemaInvalidate };
}
