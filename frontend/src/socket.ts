import { io } from "socket.io-client";
import { acquireToken, type AuthResult } from "./utils/auth";
import { getShareToken } from "./utils/shareToken";

export const socket = io(undefined, { autoConnect: false });

/**
 * (Re)connect the socket with a valid JWT.
 *
 * Always acquires a fresh token via acquireToken(), which validates
 * the token in localStorage or creates a new guest token.
 * If the current URL path contains a room id and a share token is
 * known for that room, it is included in the auth payload so the
 * server can grant access accordingly.
 */
export async function connectWithAuth(): Promise<AuthResult> {
	const result = await acquireToken();
	if (socket.connected) socket.disconnect();
	const roomIdMatch = /^\/rooms\/([^/]+)/.exec(window.location.pathname);
	const roomId = roomIdMatch?.[1];
	const shareToken = roomId ? getShareToken(roomId) : undefined;
	socket.auth = { token: result.token, share_token: shareToken };
	socket.connect();
	return result;
}
