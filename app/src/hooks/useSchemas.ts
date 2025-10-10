import { getFrameMetadata, getSchemas, getExtensionData, submitExtension as submitExtensionApi, listJobs, getJob } from "../myapi/client";
import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";

export interface ExtensionMetadata {
  schema: any;
  provider: "celery" | number; // "celery" for server-side, or number of workers for client-side
  queueLength: number;
  idleWorkers: number;
  progressingWorkers: number;
}

export interface SchemasResponse {
  [extensionName: string]: ExtensionMetadata;
}

export const useSchemas = (room: string, category: string) => {
  return useQuery({
    queryKey: ["schemas", room, category],
    queryFn: () => getSchemas(room, category),
    enabled: !!room && !!category, // Only run the query if room and category are available
  });
};

export const useExtensionData = (
  room: string,
  user: string,
  category: string,
  extension: string,
) => {
  const initialData = {
    "studio_lighting": {
      "key_light_intensity": 1.2,
      "fill_light_intensity": 0.3,
      "rim_light_intensity": 1.5,
      "background_color": "#333840",
      "contact_shadow": true
    },
    "camera": {
      "camera": "PerspectiveCamera"
    }
  }

  return useQuery({
    queryKey: ["extensionData", room, user, category, extension],
    queryFn: async () => {
      const result = await getExtensionData(room, user, category, extension);
      console.log(`Fetched extension ${category}/${extension} data:`, result);
      return result.data;
    },
    // Remove staleTime: Infinity if you expect this data to be updated externally
    initialData: initialData[extension] || undefined,
    staleTime: Infinity,
    enabled: !!room && !!user && !!category && !!extension,
  });
};

export const useSubmitExtension = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: async (variables: {
      roomId: string;
      userId: string;
      category: string;
      extension: string;
      data: any;
    }) => {
      return submitExtensionApi(variables);
    },
    // On success, MANUALLY update the query's cached data
    onSuccess: (serverData, variables) => {
      console.log("Extension submitted successfully:", serverData);
      const {
        roomId,
        userId,
        category,
        extension,
        data: submittedData,
      } = variables;

      // Define the exact query key for the data we want to update
      const queryKey = ["extensionData", roomId, userId, category, extension];

      // Replace the cached data with the data the user just submitted.
      // This is an instant, local update with no network request.
      queryClient.setQueryData(queryKey, submittedData);
    },
    onError: (error) => {
      console.error("Error submitting extension:", error);
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
  status: "queued" | "running" | "completed" | "failed";
  provider: string;
  created_at: string;
  started_at: string;
  completed_at: string;
  worker_id: string;
  error: string;
  result: any;
}

export const useJobs = (room: string) => {
  return useQuery({
    queryKey: ["jobs", room],
    queryFn: async () => {
      const result = await listJobs(room);
      return result.jobs;
    },
    enabled: !!room,
    refetchInterval: 5000, // Refetch every 5 seconds
  });
};

export const useJob = (room: string, jobId: string) => {
  return useQuery({
    queryKey: ["job", room, jobId],
    queryFn: async () => {
      const result = await getJob(room, jobId);
      return result.job;
    },
    enabled: !!room && !!jobId,
  });
};


export const useFrameMetadata = (roomId: string, frameId: number = 0) => {
  return useQuery({
    queryKey: ['metadata', roomId, frameId], // TODO: need to invalidate!
    queryFn: () => getFrameMetadata(roomId, frameId),
    enabled: !!roomId, // Only run the query if a roomId is available
    staleTime: 1000 * 60 * 5, // Metadata for frame 0 is unlikely to change often
  });
};