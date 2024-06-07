import { io } from "socket.io-client";
import { useEffect, useRef } from "react";

import * as THREE from "three";

// TODO: get the server url from the environment
const serverUrl = window.location.origin;
// const serverUrl = "http://localhost:1235"; // for local development
export const socket = io(serverUrl);

export const sendStep = (step: number, fromSockets: any) => {
  const timeoutRef = useRef(null);

  useEffect(() => {
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current);
    }
    if (fromSockets.current) {
      fromSockets.current = false;
    } else {
      timeoutRef.current = setTimeout(() => {
        socket.emit("room:step:set", step);
      }, 100);
    }

    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
    };
  }, [step]);
};

export const sendCamera = (data: { position: number[]; target: number[] }) => {
  const timeoutRef = useRef(null);

  useEffect(() => {
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current);
    }

    timeoutRef.current = setTimeout(() => {
      socket.emit("room:camera:set", data);
    }, 100);

    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
    };
  }, [data]);
};

export const sendSelection = (selection: Set<number>, fromSockets: any) => {
  // TODO: only sync if no change has happend for 100 milliseconds
  useEffect(() => {
    if (fromSockets.current) {
      fromSockets.current = false;
    } else {
      socket.emit("room:selection:set", { 0: Array.from(selection) });
    }
  }, [selection]);
};

export const sendBookmarks = (
  bookmarks: { number: string },
  fromSockets: any,
) => {
  useEffect(() => {
    if (fromSockets.current) {
      fromSockets.current = false;
    } else {
      socket.emit("room:bookmarks:set", bookmarks);
    }
  }, [bookmarks]);
};

export const sendPoints = (points: THREE.Vector3[], fromSockets: any) => {
  const timeoutRef = useRef(null);

  useEffect(() => {
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current);
    }

    if (fromSockets.current) {
      fromSockets.current = false;
    } else {
      timeoutRef.current = setTimeout(() => {
        const message = { 0: points.map((point) => point.toArray()) };
        socket.emit("room:points:set", message);
      }, 100);
    }
  }, [points]);
};
