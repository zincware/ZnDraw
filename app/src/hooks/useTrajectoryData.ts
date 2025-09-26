import { useCallback } from 'react';
import { useAppStore } from '../store';
import { decode } from '@msgpack/msgpack';

export const useTrajectoryData = () => {
  const { roomId, setFrameData, setLoadingFrame } = useAppStore();

  const loadFrame = useCallback(async (frameIndex: number, keys: string[]) => {
      fetch(`/api/frames/${roomId}`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ indices: [frameIndex], keys: keys }),
      })
        .then((res) => {
          if (res.ok) {
            return res.arrayBuffer();
          }
          console.error("Received error response from /api/frames:", res);
          throw new Error(`Request failed with status ${res.status}`);
        })
        .then((arrayBuffer) => {
          const decodedData = decode(arrayBuffer);
          console.log("Received data from /api/frames:", decodedData);
          setFrameData(frameIndex, decodedData as any);
          setLoadingFrame(false);
        })
        .catch((error) => {
          console.error("Error fetching from /api/frames:", error);
          setLoadingFrame(false);
        });
    
  }, [roomId, setFrameData, setLoadingFrame]);
  
  return { loadFrame };
};