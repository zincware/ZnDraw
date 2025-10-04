import { useEffect } from 'react';
import { socket } from '../socket';
import { useAppStore } from '../store';
import { useQueryClient } from '@tanstack/react-query';

export const useSocketManager = () => {
  const { setConnected, setFrameCount, isConnected, setCurrentFrame, setFrameSelection, setSelection, setBookmarks, roomId, userId, joinToken, setGeometries } = useAppStore();
  const queryClient = useQueryClient();

  useEffect(() => {
    if (!joinToken) {
      console.warn('No join token available, cannot connect socket.');
      return;
    }

    async function onConnect() {
      console.log('Socket connected and joining room:',   roomId, userId);
      setConnected(true);

    }
    function onDisconnect() {
      console.log('Socket disconnected');
      setConnected(false);
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
      queryClient.invalidateQueries({ queryKey: ['jobs', roomId] });

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
          // Query keys are expected to be in the format: ['frame', roomId, frameIndex]
          const [type, qRoomId, frameIndex] = query.queryKey;

          // Only invalidate frame queries for this room
          if (type !== 'frame' || qRoomId !== roomId) return false;

          // Ensure we are only dealing with queries for specific frames
          if (typeof frameIndex !== 'number') return false;

          if (operation === 'replace') {
            // Handles single-item replacement ONLY.
            return frameIndex === affectedIndex;
          } else if (
            operation === 'delete' ||
            operation === 'insert' ||
            operation === 'bulk_replace'
          ) {
            // Handles any operation that shifts indices or modifies a range.
            // Invalidate all frames from the affected position onward.
            return frameIndex >= affectedFrom;
          }

          // Default to not invalidating if the operation is unknown.
          return false;
        }
      });
    }

    function onChatMessageNew(data: any) {
      queryClient.setQueryData(['chat', roomId], (oldData: any) => {
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
      queryClient.setQueryData(['chat', roomId], (oldData: any) => {
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

    function onGeometriesInvalidate(data: any) {
      // get new geometries from data
      fetch('/api/rooms/' + roomId + '/geometries')
        .then(response => {
          if (!response.ok) {
            throw new Error(`Failed to fetch geometries: ${response.status} ${response.statusText}`);
          }
          return response.json();
        })
        .then(geometries => {
          console.log('Received updated geometries:', geometries);
          setGeometries(geometries);
        })
        .catch(error => {
          console.error('Error fetching geometries:', error);
        });
    }

    function onConnectError(err: any) {
      console.error('Socket connection error:', err);
    }

    socket.on('disconnect', onDisconnect);
    socket.on('connect', onConnect);
    socket.on('connect_error', onConnectError);
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
    socket.on('invalidate:geometry', onGeometriesInvalidate);

    socket.auth = { token: joinToken };
    socket.connect();

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
      socket.off('invalidate:geometry', onGeometriesInvalidate);
    };
  }, [joinToken, roomId, userId, setConnected, setFrameCount, userId, isConnected, setCurrentFrame, queryClient, setBookmarks, setSelection, setFrameSelection, setGeometries]);
};
