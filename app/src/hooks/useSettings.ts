// hooks/useSettings.ts
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';

// This function now needs the userId to fetch the correct settings
const fetchSettingsData = async (roomId: string, userId: string) => {
  // You would adjust your API endpoint to accept a userId,
  // perhaps in the header or as a query param. For now, we assume the backend knows.
  const response = await fetch(`/api/rooms/${roomId}/settings`);
  if (!response.ok) throw new Error('Failed to fetch settings');
  return response.json();
};

const saveSettingsData = async (roomId: string, userId: string, data: any) => {
    // Similarly, the POST request identifies the user
    const response = await fetch(`/api/rooms/${roomId}/settings`, { /* ... */ });
    if (!response.ok) throw new Error('Failed to save settings');
    return response.json();
};

export const useSettingsData = (roomId: string, userId: string, { enabled = true } = {}) => {
  return useQuery({
    queryKey: ['settings', 'data', roomId, userId],
    queryFn: () => fetchSettingsData(roomId, userId),
    enabled: !!roomId && !!userId && enabled,
  });
};

export const useSaveSettings = (roomId: string, userId: string) => {
  const queryClient = useQueryClient();
  return useMutation({
    mutationFn: (data: any) => saveSettingsData(roomId, userId, data),
    onSuccess: () => {
      queryClient.invalidateQueries({ queryKey: ['settings', 'data', roomId, userId] });
    },
  });
};