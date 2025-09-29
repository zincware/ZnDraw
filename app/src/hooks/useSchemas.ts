import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';

const fetchSchemas = async (room: string, option: string) => {
  const response = await fetch(`/api/rooms/${room}/schema/${option}`);
  if (!response.ok) {
    throw new Error(`Failed to fetch schema for ${option}`);
  }
  return response.json();
};

export const useSchemas = (room: string, option: string) => {
  return useQuery({
    queryKey: ['schemas', room, option],
    queryFn: () => fetchSchemas(room, option),
    staleTime: Infinity, // These schemas are static, they never go stale
    enabled: !!room && !!option, // Only run the query if room and option are available
  });
};

const fetchSchemaData = async (room: string, user: string, option: string, action: string) => {
  const response = await fetch(`/api/rooms/${room}/actions-data?userId=${user}&room=${room}&action=${action}&method=${option}`, {
    method: 'GET',
    headers: {
      'Content-Type': 'application/json',
    },
  });
  if (!response.ok) {
    throw new Error(`Failed to fetch schema data for ${option}`);
  }
  const result = await response.json();
  console.log(`Fetched schema ${action}/${option} data:`, result);
  return result["data"];
};

export const useSchemaData = (room: string, user: string, action: string, option: string) => {
  return useQuery({
    queryKey: ['schemaData', room, user, action, option], // The key is critical
    queryFn: () => fetchSchemaData(room, user, option, action),
    // Remove staleTime: Infinity if you expect this data to be updated externally
    staleTime: 1000 * 60, // e.g., stale after 1 minute
    enabled: !!room && !!user && !!option && !!action,
  });
};

// NEW: Mutation for submitting the form data
const submitAction = async (variables: {
  roomId: string;
  userId: string;
  action: string;
  method: string;
  data: any;
}) => {
  const { roomId, ...payload } = variables;
  const response = await fetch(`/api/rooms/${roomId}/actions`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload),
  });
  if (!response.ok) {
    throw new Error('Network response was not ok');
  }
  return response.json();
};

export const useSubmitAction = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: submitAction,
    // On success, MANUALLY update the query's cached data
    onSuccess: (serverData, variables) => {
      console.log('Action submitted successfully:', serverData);
      const { roomId, userId, action, method, data: submittedData } = variables;

      // Define the exact query key for the data we want to update
      const queryKey = ['schemaData', roomId, userId, action, method];

      // Replace the cached data with the data the user just submitted.
      // This is an instant, local update with no network request.
      queryClient.setQueryData(queryKey, submittedData);
    },
    onError: (error) => {
      console.error('Error submitting action:', error);
    },
  });
};
