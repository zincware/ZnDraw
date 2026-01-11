import { io } from "socket.io-client";
import { getToken } from "./utils/auth";

const URL =
	process.env.NODE_ENV === "production" ? undefined : "http://localhost:5000";

export const socket = io(URL, {
	autoConnect: false,
	auth: (cb) => {
		// Send JWT token in auth on connect/reconnect
		// SessionId is created in room:join handler, not at connect time
		const token = getToken();
		cb({ token });
	},
});
