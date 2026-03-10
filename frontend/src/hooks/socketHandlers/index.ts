/**
 * Barrel re-exports for socket handler factories and types.
 *
 * Each handler factory accepts a HandlerContext and returns an object of
 * event handler functions. useSocketManager calls these factories inside
 * its useEffect, registers the handlers on socket.on(), and deregisters
 * them in the cleanup function.
 */

// Types
export type { HandlerContext } from "./types";

// Utility
export { createInvalidateHandler } from "./utils";

// Handler factories
export { createConnectionHandlers } from "./connectionHandlers";
export { createFrameHandlers } from "./frameHandlers";
export { createGeometryHandlers } from "./geometryHandlers";
export { createChatHandlers } from "./chatHandlers";
export { createSceneInvalidationHandlers } from "./sceneInvalidationHandlers";
export { createFigureHandlers } from "./figureHandlers";
export { createRoomHandlers } from "./roomHandlers";
