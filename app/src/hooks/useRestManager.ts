import { useEffect, useMemo, useCallback, useRef } from 'react';
import { useAppStore } from '../store';
import { useParams, useSearchParams } from 'react-router-dom';
import { set, throttle } from 'lodash';
import { useQueryClient } from '@tanstack/react-query';

export const useRestJoinManager = () => {
    const { setClientId, setRoomId, setUserId, setCurrentFrame, setFrameCount, setSelection, setFrameSelection, setBookmarks, setJoinToken } = useAppStore();
    const { roomId: room, userId } = useParams<{ roomId: string, userId: string }>();
    const [searchParams] = useSearchParams();
    const queryClient = useQueryClient();

    const abortControllerRef = useRef<AbortController | null>(null);

    const joinRoom = useCallback(async () => {
        if (!room || !userId) {
            return;
        }

        // 2. Create a new AbortController for this specific request.
        const controller = new AbortController();
        abortControllerRef.current = controller; // Store it in the ref

        console.log('Joining room via REST:', room, userId);

        // Get template from query parameters
        const template = searchParams.get('template');

        // Build request body
        const requestBody: { userId: string; template?: string } = { userId };
        if (template) {
            requestBody.template = template;
        }

        try {
            const response = await fetch(`/api/rooms/${room}/join`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(requestBody),
                // 3. Pass the controller's signal to the fetch options.
                signal: controller.signal,
            });

            if (!response.ok) {
                console.error(`Failed to join room: ${response.status} ${response.statusText}`);
                return;
            }

            const data = await response.json();
            console.log('Join response data:', data);

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
            if (data.clientId) {
                setClientId(data.clientId);
            }
            if (data.joinToken) {
                setJoinToken(data.joinToken);
            }
            setRoomId(room);
            setUserId(userId);

            // Store settings in query cache for each category
            if (data.settings) {
                for (const [categoryName, categoryData] of Object.entries(data.settings)) {
                    const queryKey = ['extensionData', room, userId, 'settings', categoryName];
                    queryClient.setQueryData(queryKey, categoryData);
                }
            }

        } catch (error) {
            // 4. Check if the error was due to the request being aborted.
            if (error.name === 'AbortError') {
                console.log('Fetch aborted on unmount or re-run.');
            } else {
                console.error('Error joining room:', error);
            }
        } finally {
            // Clean up the ref if the controller for this fetch is the current one
            if (abortControllerRef.current === controller) {
                abortControllerRef.current = null;
            }
        }

    }, [room, userId, searchParams, queryClient, setClientId, setRoomId, setUserId, setCurrentFrame, setFrameCount, setSelection, setFrameSelection, setBookmarks, setJoinToken]);

    const throttledJoin = useMemo(() =>
        throttle(joinRoom, 10000, { leading: true, trailing: false }),
        [joinRoom]
    );

    useEffect(() => {
        throttledJoin();

        return () => {
            abortControllerRef.current?.abort(); // Abort the fetch
            throttledJoin.cancel(); // Cancel any pending throttled execution
        };
    }, [room, userId, throttledJoin]);

    return {};
};