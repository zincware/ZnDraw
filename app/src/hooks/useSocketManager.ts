import { useEffect } from 'react';
import { socket } from '../socket';
import { useAppStore } from '../store';
import { useParams } from 'react-router-dom';
import { useQueryClient } from '@tanstack/react-query';


export const useSocketManager = () => {
  const { roomId: room, userId } = useParams<{ roomId: string, userId: string }>();
  const { setConnected, setFrameCount, isConnected, setPresenter, setPresenterSid, setCurrentFrame } = useAppStore();
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

    function onPresenterTokenGranted() {
      console.log("I am now the presenter!");
      setPresenter(true);
      // Start a timer to renew the token
      // startTokenRenewal();
    }
    function onPresenterTokenDenied() {
      console.log("Someone else is presenting.");
      setPresenter(false);
    }
    function onPresenterUpdate(data: any) {
      // data = { presenterSid: 'some-session-id' } or { presenterSid: null }
      console.log('Presenter update:', data);
      setPresenterSid(data.presenterSid);
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

    socket.on('disconnect', onDisconnect);
    socket.on('connect', onConnect);
    socket.on('len_frames', onLenUpdate);

    socket.on('presenter_token_granted', onPresenterTokenGranted);
    socket.on('presenter_token_denied', onPresenterTokenDenied);
    socket.on('presenter_update', onPresenterUpdate);
    socket.on('frame_update', onFrameUpdate);
    socket.on('invalidate', onInvalidate);

    return () => {
      socket.off('connect', onConnect);
      socket.off('disconnect', onDisconnect);
      socket.off('len_frames', onLenUpdate);

      socket.off('presenter_token_granted', onPresenterTokenGranted);
      socket.off('presenter_token_denied', onPresenterTokenDenied);
      socket.off('presenter_update', onPresenterUpdate);
      socket.off('frame_update', onFrameUpdate);
      socket.off('invalidate', onInvalidate);
    };
  }, [room, setConnected, setFrameCount, userId, isConnected, setPresenter, setPresenterSid, setCurrentFrame, queryClient]);
};
