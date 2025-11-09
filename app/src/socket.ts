import { io } from "socket.io-client";
import { getToken } from "./utils/auth";
import { useAppStore } from "./store"


const URL = process.env.NODE_ENV === 'production' ? undefined : 'http://localhost:5000';

export const socket = io(URL, {
  autoConnect: false,
  auth: (cb) => {
    // Send JWT token in auth on connect/reconnect
    // This callback is called every time the socket connects
    const token = getToken();
    const sessionId = useAppStore.getState().sessionId;
    cb({ token , sessionId });
  }
});
