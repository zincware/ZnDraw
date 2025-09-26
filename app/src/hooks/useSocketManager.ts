import { useEffect } from 'react';
import { socket } from '../socket';
import { useAppStore } from '../store';

export const useSocketManager = (room: string) => {
  const { setConnected, setFrameCount } = useAppStore();

  useEffect(() => {
    function onConnect() {
      setConnected(true, room);
      socket.emit('join_room', { room });
      socket.emit('frames:count', {}, (response: any) => {
        if (!response || !response.success) {
          console.error(`Failed to get frame count: ${response ? response.error : 'No response'}`);
          return;
        }
        setFrameCount(response.count);
      });

    }    
    socket.on('connect', onConnect);

    return () => {
      socket.off('connect', onConnect);
    };
  }, [room, setConnected, setFrameCount]);
};