import { useEffect } from 'react';
import { socket } from '../socket';
import { useAppStore } from '../store';
import { useParams } from 'react-router-dom';
import { useQueryClient } from '@tanstack/react-query';

export const useSocketManager = () => {
  const { roomId: room, userId } = useParams<{ roomId: string, userId: string }>();
  const { setConnected, setFrameCount, isConnected, setCurrentFrame } = useAppStore();
  const queryClient = useQueryClient();

  useEffect(() => {
    if (!room) return;
    if (isConnected) return;
    socket.connect();
  }, [room]);

  useEffect(() => {
    function onConnect() {
      console.log('Socket connected and joining room:', room);
      setConnected(true, room, userId);
      socket.emit('join_room', { room, userId });
    }
    function onDisconnect() {
      console.log('Socket disconnected');
      setConnected(false, room, userId);
    }
    function onLenUpdate(data: any) {
      if (data && typeof data.count === 'number') {
        setFrameCount(data.count);
      } else {
        console.error('Invalid len_frames data:', data);
      }
    }

    function onFrameUpdate(data: any) {
      const { frame } = data;
      setCurrentFrame(frame);
    }

    function onInvalidate(data: any) {
      const { roomId, userId, category, extension } = data;
      queryClient.invalidateQueries({
          queryKey: ['extensionData', roomId, userId, category, extension],
        });
      console.log(`Invalidated extension data for user ${userId}, category ${category}, extension ${extension} in room ${roomId}`);
    }

    function onSchemaInvalidate(data: any) {
      const { roomId, category } = data;
      queryClient.invalidateQueries({
          queryKey: ['schemas', roomId, category],
        });
      console.log(`Invalidated schemas for category ${category} in room ${roomId}`);
    }

    function onQueueUpdate(data: any) {
      const { roomId, category, extension, queueLength, idleWorkers, progressingWorkers } = data;
      // Update the cached schema data with new queue length and worker counts
      queryClient.setQueryData(['schemas', roomId, category], (oldData: any) => {
        if (!oldData || !oldData[extension]) return oldData;
        return {
          ...oldData,
          [extension]: {
            ...oldData[extension],
            queueLength,
            idleWorkers,
            progressingWorkers,
          },
        };
      });
      console.log(`Queue updated for ${category}/${extension}: ${queueLength} queued, ${idleWorkers} idle, ${progressingWorkers} progressing`);
    }

    socket.on('disconnect', onDisconnect);
    socket.on('connect', onConnect);
    socket.on('len_frames', onLenUpdate);
    socket.on('frame_update', onFrameUpdate);
    socket.on('invalidate', onInvalidate);
    socket.on('invalidate:schema', onSchemaInvalidate);
    socket.on('queue:update', onQueueUpdate);

    return () => {
      socket.off('connect', onConnect);
      socket.off('disconnect', onDisconnect);
      socket.off('len_frames', onLenUpdate);
      socket.off('frame_update', onFrameUpdate);
      socket.off('invalidate', onInvalidate);
      socket.off('invalidate:schema', onSchemaInvalidate);
      socket.off('queue:update', onQueueUpdate);
    };
  }, [room, setConnected, setFrameCount, userId, isConnected, setCurrentFrame, queryClient]);
};
