import { useEffect } from 'react';
import { socket } from '../socket';
import { useAppStore } from '../store';
import { useParams } from 'react-router-dom';
import { useQueryClient } from '@tanstack/react-query';

export const useSocketManager = () => {
  const { roomId: room, userId } = useParams<{ roomId: string, userId: string }>();
  const { setConnected, setFrameCount, isConnected, setCurrentFrame, setFrameSelection, setSelection, setBookmarks } = useAppStore();
  const queryClient = useQueryClient();

  useEffect(() => {
    if (!room) return;
    if (isConnected) return;
    socket.connect();
  }, [room, isConnected]);

  useEffect(() => {
    async function onConnect() {
      console.log('Socket connected and joining room:', room);
      setConnected(true, room || '', userId || '');
      socket.emit('join_room', { room, userId });

      // Post to join endpoint to get initial room data
      try {
        const response = await fetch(`/api/rooms/${room}/join`, {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({}),
        });

        if (!response.ok) {
          console.error(`Failed to join room: ${response.status} ${response.statusText}`);
          return;
        }

        const data = await response.json();

        // Update Zustand store with room data
        if (typeof data.frameCount === 'number') {
          setFrameCount(data.frameCount);
        }
        if (data.selection !== undefined) {
          setSelection(data.selection);
        }
        if (data.frame_selection !== undefined) {
          setFrameSelection(data.frame_selection);
        }
        if (data.step !== undefined) {
          setCurrentFrame(data.step);
        }
        if (data.bookmarks !== undefined) {
          setBookmarks(data.bookmarks);
        }

        console.log('Room joined successfully:', data);
      } catch (error) {
        console.error('Error joining room:', error);
      }
    }
    function onDisconnect() {
      console.log('Socket disconnected');
      setConnected(false, room || '', userId || '');
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

      // 1. Update the cached schema data with new queue length and worker counts
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

      // 2. Invalidate jobs list to refetch all jobs (covers created/started/completed/failed/deleted)
      queryClient.invalidateQueries({ queryKey: ['jobs', room] });

      console.log(`Queue updated for ${category}/${extension}: ${queueLength} queued, ${idleWorkers} idle, ${progressingWorkers} progressing`);
    }

    function onSelectionUpdate(data: any) {
      console.log(data);
      setSelection(data["indices"] || null);
    }

    function onFrameSelectionUpdate(data: any) {
      console.log(data);
      setFrameSelection(data["indices"] || null);
    }

    function onBookmarksUpdate(data: any) {
      console.log('Bookmarks update:', data);
      setBookmarks(data["bookmarks"] || null);
    }

    function onFramesInvalidate(data: any) {
      const { roomId, operation, affectedIndex, affectedFrom } = data;
      console.log('Invalidating frame cache:', { roomId, operation, affectedIndex, affectedFrom });

      queryClient.invalidateQueries({
        predicate: (query) => {
          const [type, qRoomId, frameIndex] = query.queryKey;

          // Only invalidate frame queries for this room
          if (type !== 'frame' || qRoomId !== roomId) return false;

          // frameIndex might be undefined for non-frame queries
          if (typeof frameIndex !== 'number') return false;

          if (operation === 'replace') {
            // Only invalidate the specific replaced frame
            return frameIndex === affectedIndex;
          } else if (operation === 'delete' || operation === 'insert') {
            // Invalidate all frames from the affected position onward
            return frameIndex >= affectedFrom;
          }

          return false;
        }
      });
    }

    function onChatMessageNew(data: any) {
      queryClient.setQueryData(['chat', room], (oldData: any) => {
        if (!oldData) return oldData;
        const newPages = [...oldData.pages];
        const lastPageIndex = newPages.length - 1;
        if (lastPageIndex >= 0) {
          newPages[lastPageIndex] = {
            ...newPages[lastPageIndex],
            messages: [...newPages[lastPageIndex].messages, data],
            metadata: {
              ...newPages[lastPageIndex].metadata,
              totalCount: newPages[lastPageIndex].metadata.totalCount + 1,
            }
          };
        }
        return { ...oldData, pages: newPages };
      });
    }

    function onChatMessageUpdated(data: any) {
      queryClient.setQueryData(['chat', room], (oldData: any) => {
        if (!oldData) return oldData;
        const newPages = oldData.pages.map((page: any) => ({
          ...page,
          messages: page.messages.map((msg: any) =>
            msg.id === data.id ? data : msg
          ),
        }));
        return { ...oldData, pages: newPages };
      });
    }

    socket.on('disconnect', onDisconnect);
    socket.on('connect', onConnect);
    socket.on('len_frames', onLenUpdate);
    socket.on('frame_update', onFrameUpdate);
    socket.on('invalidate', onInvalidate);
    socket.on('invalidate:schema', onSchemaInvalidate);
    socket.on('queue:update', onQueueUpdate);
    socket.on('selection:update', onSelectionUpdate);
    socket.on('frame_selection:update', onFrameSelectionUpdate);
    socket.on('bookmarks:update', onBookmarksUpdate);
    socket.on('frames:invalidate', onFramesInvalidate);
    socket.on('chat:message:new', onChatMessageNew);
    socket.on('chat:message:updated', onChatMessageUpdated);

    return () => {
      socket.off('connect', onConnect);
      socket.off('disconnect', onDisconnect);
      socket.off('len_frames', onLenUpdate);
      socket.off('frame_update', onFrameUpdate);
      socket.off('invalidate', onInvalidate);
      socket.off('invalidate:schema', onSchemaInvalidate);
      socket.off('queue:update', onQueueUpdate);
      socket.off('selection:update', onSelectionUpdate);
      socket.off('frame_selection:update', onFrameSelectionUpdate);
      socket.off('bookmarks:update', onBookmarksUpdate);
      socket.off('frames:invalidate', onFramesInvalidate);
      socket.off('chat:message:new', onChatMessageNew);
      socket.off('chat:message:updated', onChatMessageUpdated);
    };
  }, [room, setConnected, setFrameCount, userId, isConnected, setCurrentFrame, queryClient, setBookmarks]);
};
