import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';

const fetchSchemas = async (room: string, category: string) => {
  const response = await fetch(`/api/rooms/${room}/schema/${category}`);
  if (!response.ok) {
    throw new Error(`Failed to fetch schema for ${category}`);
  }
  return response.json();
};

export const useSchemas = (room: string, category: string) => {
  return useQuery({
    queryKey: ['schemas', room, category],
    queryFn: () => fetchSchemas(room, category),
    staleTime: Infinity, // These schemas are static, they never go stale
    enabled: !!room && !!category, // Only run the query if room and category are available
  });
};

const fetchExtensionData = async (room: string, user: string, category: string, extension: string) => {
  const response = await fetch(`/api/rooms/${room}/extension-data/${category}/${extension}?userId=${user}`, {
    method: 'GET',
    headers: {
      'Content-Type': 'application/json',
    },
  });
  if (!response.ok) {
    throw new Error(`Failed to fetch extension data for ${category}/${extension}`);
  }
  const result = await response.json();
  console.log(`Fetched extension ${category}/${extension} data:`, result);
  return result["data"];
};

export const useExtensionData = (room: string, user: string, category: string, extension: string) => {
  return useQuery({
    queryKey: ['extensionData', room, user, category, extension], // The key is critical
    queryFn: () => fetchExtensionData(room, user, category, extension),
    // Remove staleTime: Infinity if you expect this data to be updated externally
    staleTime: 1000 * 60, // e.g., stale after 1 minute
    enabled: !!room && !!user && !!category && !!extension,
  });
};

// NEW: Mutation for submitting the extension data
const submitExtension = async (variables: {
  roomId: string;
  userId: string;
  category: string;
  extension: string;
  data: any;
}) => {
  const { roomId, userId, category, extension, data } = variables;
  const response = await fetch(`/api/rooms/${roomId}/extensions/${category}/${extension}?userId=${userId}`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  });
  if (!response.ok) {
    throw new Error('Network response was not ok');
  }
  return response.json();
};

export const useSubmitExtension = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: submitExtension,
    // On success, MANUALLY update the query's cached data
    onSuccess: (serverData, variables) => {
      console.log('Extension submitted successfully:', serverData);
      const { roomId, userId, category, extension, data: submittedData } = variables;

      // Define the exact query key for the data we want to update
      const queryKey = ['extensionData', roomId, userId, category, extension];

      // Replace the cached data with the data the user just submitted.
      // This is an instant, local update with no network request.
      queryClient.setQueryData(queryKey, submittedData);
    },
    onError: (error) => {
      console.error('Error submitting extension:', error);
    },
  });
};
