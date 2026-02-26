import { io } from "socket.io-client";
import { acquireToken, type AuthResult } from "./utils/auth";

export const socket = io(undefined, { autoConnect: false });

/**
 * (Re)connect the socket with a valid JWT.
 *
 * Always acquires a fresh token via acquireToken(), which validates
 * the token in localStorage or creates a new guest token.
 */
export async function connectWithAuth(): Promise<AuthResult> {
	const result = await acquireToken();
	if (socket.connected) socket.disconnect();
	socket.auth = { token: result.token };
	socket.connect();
	return result;
}
