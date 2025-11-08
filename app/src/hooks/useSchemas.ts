import { getFrameKeys, getFrameMetadata, getSchemas, getExtensionData, submitExtension as submitExtensionApi, listJobs, getJob } from "../myapi/client";
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
    queryFn: async () => {
      // Performance monitoring for schema fetch
      const metricName = `schemas-fetch-${room}-${category}`;
      performance.mark(`${metricName}-start`);

      try {
        const result = await getSchemas(room, category);

        performance.mark(`${metricName}-end`);
        performance.measure(metricName, `${metricName}-start`, `${metricName}-end`);

        const measure = performance.getEntriesByName(metricName)[0];
        if (import.meta.env.DEV && measure) {
          console.log(`[Performance] Schemas fetch (${category}) took ${Math.round(measure.duration)}ms`);
        }

        // Clean up marks to avoid memory leaks
        performance.clearMarks(`${metricName}-start`);
        performance.clearMarks(`${metricName}-end`);
        performance.clearMeasures(metricName);

        return result;
      } catch (error) {
        // Clean up marks even on error
        performance.clearMarks(`${metricName}-start`);
        performance.clearMarks(`${metricName}-end`);
        throw error;
      }
    },
    enabled: !!room && !!category, // Only run the query if room and category are available
  });
};

export const useExtensionData = (
  room: string,
  user: string,
  category: string,
  extension: string,
) => {
  const initialData: Record<string, any> = {
    "studio_lighting": {
      "ambient_light": 0.35,
      "key_light": 0.7,
      "fill_light": 0.4,
      "rim_light": 0.5,
      "hemisphere_light": 0.3,
      "background_color": "default",
    },
    "camera": {
      "camera": "PerspectiveCamera",
      "show_crosshair": false,
      "far_plane": 300,
      "near_plane": 1
    }
  }

  return useQuery({
    queryKey: ["extensionData", room, user, category, extension],
    queryFn: async () => {
      const result = await getExtensionData(room, category, extension);
      if (import.meta.env.DEV) {
        console.log(`Fetched extension ${category}/${extension} data:`, result);
      }
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
      userName: string;
      category: string;
      extension: string;
      data: any;
    }) => {
      return submitExtensionApi(variables);
    },
    // On success, MANUALLY update the query's cached data
    onSuccess: (serverData, variables) => {
      if (import.meta.env.DEV) {
        console.log("Extension submitted successfully:", serverData);
      }
      const {
        roomId,
        userName,
        category,
        extension,
        data: submittedData,
      } = variables;

      // Define the exact query key for the data we want to update
      const queryKey = ["extensionData", roomId, userName, category, extension];

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
  userName: string;
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


export const useFrameKeys = (
  roomId: string,
  frameId: number = 0,
  enabled: boolean = true
) => {
  return useQuery({
    queryKey: ['frame-keys', roomId, frameId],
    queryFn: async () => {
      return await getFrameKeys(roomId, frameId);
    },
    enabled: !!roomId && enabled,
    staleTime: 1000 * 60 * 5, // Keys rarely change
  });
};

export const useFrameMetadata = (
  roomId: string,
  frameId: number = 0,
  enabled: boolean = true
) => {
  return useQuery({
    queryKey: ['metadata', roomId, frameId], // TODO: need to invalidate!
    queryFn: async () => {
      // Performance monitoring for metadata fetch
      const metricName = `metadata-fetch-${roomId}-${frameId}`;
      performance.mark(`${metricName}-start`);

      try {
        const result = await getFrameMetadata(roomId, frameId);

        performance.mark(`${metricName}-end`);
        performance.measure(metricName, `${metricName}-start`, `${metricName}-end`);

        const measure = performance.getEntriesByName(metricName)[0];
        if (import.meta.env.DEV && measure) {
          console.log(`[Performance] Metadata fetch took ${Math.round(measure.duration)}ms`);
        }

        // Clean up marks to avoid memory leaks
        performance.clearMarks(`${metricName}-start`);
        performance.clearMarks(`${metricName}-end`);
        performance.clearMeasures(metricName);

        return result;
      } catch (error) {
        // Clean up marks even on error
        performance.clearMarks(`${metricName}-start`);
        performance.clearMarks(`${metricName}-end`);
        throw error;
      }
    },
    enabled: !!roomId && enabled, // Only run the query if roomId is available AND explicitly enabled
    staleTime: 1000 * 60 * 5, // Metadata for frame 0 is unlikely to change often
  });
};