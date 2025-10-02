import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query';

export interface ExtensionMetadata {
  schema: any;
  provider: "celery" | number;  // "celery" for server-side, or number of workers for client-side
  queueLength: number;
  idleWorkers: number;
  progressingWorkers: number;
}

export interface SchemasResponse {
  [extensionName: string]: ExtensionMetadata;
}

const fetchSchemas = async (room: string, category: string): Promise<SchemasResponse> => {
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
    enabled: !!room && !!category, // Only run the query if room and category are available
  });
};

const fetchExtensionData = async (room: string, user: string, category: string, extension: string) => {
  const response = await fetch(`/api/rooms/${room}/extensions/${category}/${extension}/data?userId=${user}`, {
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
  const response = await fetch(`/api/rooms/${roomId}/extensions/${category}/${extension}/submit`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ userId, data }),
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

// Job interfaces
export interface Job {
  id: string;
  room: string;
  category: string;
  extension: string;
  data: any;
  user_id: string;
  status: 'queued' | 'running' | 'completed' | 'failed';
  provider: string;
  created_at: string;
  started_at: string;
  completed_at: string;
  worker_id: string;
  error: string;
  result: any;
}

// Fetch active jobs for a room
const fetchJobs = async (room: string): Promise<Job[]> => {
  const response = await fetch(`/api/rooms/${room}/jobs`);
  if (!response.ok) {
    throw new Error('Failed to fetch jobs');
  }
  const result = await response.json();
  return result.jobs;
};

export const useJobs = (room: string) => {
  return useQuery({
    queryKey: ['jobs', room],
    queryFn: () => fetchJobs(room),
    enabled: !!room,
    refetchInterval: 5000, // Refetch every 5 seconds
  });
};

// Fetch a specific job
const fetchJob = async (room: string, jobId: string): Promise<Job> => {
  const response = await fetch(`/api/rooms/${room}/jobs/${jobId}`);
  if (!response.ok) {
    throw new Error('Failed to fetch job');
  }
  const result = await response.json();
  return result.job;
};

export const useJob = (room: string, jobId: string) => {
  return useQuery({
    queryKey: ['job', room, jobId],
    queryFn: () => fetchJob(room, jobId),
    enabled: !!room && !!jobId,
  });
};
