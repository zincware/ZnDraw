import { useEffect } from 'react';
import { socket } from '../socket';
import { useAppStore } from '../store';
import { useFormStore } from '../formStore';

export const useSocketManager = (room: string) => {
  const { setConnected, setFrameCount, isConnected } = useAppStore();
  const { setFormConfigs } = useFormStore();

  useEffect(() => {
    if (!room) return;
    if (isConnected) return;
    socket.connect();
  }, [room]);

  useEffect(() => {
    function onConnect() {
      console.log('Socket connected and joining room:', room);
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
    function onSchema(data: any) {
      console.log('Received schema:', data);
      setFormConfigs(data);
    }
    function onDisconnect() {
      console.log('Socket disconnected');
      setConnected(false, room);
    }

    socket.on('disconnect', onDisconnect);
    socket.on('connect', onConnect);
    socket.on('schema', onSchema);

    return () => {
      socket.off('connect', onConnect);
      socket.off('schema', onSchema);
      socket.off('disconnect', onDisconnect);
    };
  }, [room, setConnected, setFrameCount, setFormConfigs]);
};