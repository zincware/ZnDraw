import { useState, useEffect, useCallback } from "react";
import { getFrameKeys, getFrameMetadata, getSchemas, getExtensionData, submitExtension as submitExtensionApi, listJobs, getJob, type Job } from "../myapi/client";
import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import { useAppStore } from "../store";
import { type JobStatus, type ExtensionStats } from "../types/jobs";
import { socket } from "../socket";

export type { Job };

export interface ExtensionMetadata {
  schema: any;
  provider: "celery" | number; // "celery" for server-side, or number of workers for client-side
  public: boolean; // Whether this is a global/public extension
  name: string;
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
  category: string,
  extension: string,
) => {
  // Get userName from app store instead of props (data is per-user, but API uses JWT)
  const userName = useAppStore((state) => state.userName);

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
    queryKey: ["extensionData", room, userName, category, extension],
    queryFn: async () => {
      const result = await getExtensionData(room, category, extension);
      return result.data;
    },
    // Remove staleTime: Infinity if you expect this data to be updated externally
    initialData: initialData[extension] || undefined,
    staleTime: Infinity,
    enabled: !!room && !!userName && !!category && !!extension,
  });
};

export const useSubmitExtension = () => {
  const queryClient = useQueryClient();
  // Get userName from app store for cache invalidation (API uses JWT)
  const userName = useAppStore((state) => state.userName);

  return useMutation({
    mutationFn: async (variables: {
      roomId: string;
      category: string;
      extension: string;
      data: any;
      isPublic?: boolean; // Whether this is a global/public extension
    }) => {
      return submitExtensionApi(variables);
    },
    // On success, MANUALLY update the query's cached data
    onSuccess: (serverData, variables) => {
      const {
        roomId,
        category,
        extension,
        data: submittedData,
      } = variables;

      // Define the exact query key for the data we want to update
      // userName is from app store (extension data is per-user)
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

export const useJobs = (room: string) => {
  const [jobs, setJobs] = useState<any[]>([]);
  const [isLoading, setIsLoading] = useState(true); // Start with true for initial load
  const [hasLoaded, setHasLoaded] = useState(false); // Track if we've loaded once

  const refetch = useCallback(async (showLoading = false) => {
    if (!room) return;

    // Only show loading skeleton on initial load
    if (showLoading) {
      setIsLoading(true);
    }

    try {
      const result = await listJobs(room);
      const jobsList = result.jobs || [];
      setJobs(jobsList);
      setHasLoaded(true);
    } catch (error) {
      console.error("[useJobs] Error fetching jobs:", error);
      setJobs([]); // Set empty array on error
    } finally {
      if (showLoading) {
        setIsLoading(false);
      }
    }
  }, [room]);

  // Fetch jobs on mount (show loading skeleton)
  useEffect(() => {
    if (room) {
      refetch(true); // Show loading on initial mount
    }
  }, [room, refetch]);

  // Register socket handler for job state changes
  useEffect(() => {
    if (!room) return;

    const handleJobStateChanged = (data: any) => {
      const { room: jobRoom, jobId, status, metadata } = data;

      // Only update if this is for our room
      if (jobRoom !== room) return;

      // Optimistic update: immediately update the job status in local state
      setJobs((prevJobs) => {
        const jobIndex = prevJobs.findIndex((job) => job.id === jobId);
        if (jobIndex !== -1) {
          // Update existing job
          const updatedJobs = [...prevJobs];
          updatedJobs[jobIndex] = {
            ...updatedJobs[jobIndex],
            status,
            ...metadata,
          };
          return updatedJobs;
        }
        // Job not found locally, will be fetched
        return prevJobs;
      });

      // Refetch in background to sync with backend (no loading spinner)
      refetch(false);
    };

    // Register socket listener
    socket.on("job:state_changed", handleJobStateChanged);

    // Cleanup on unmount
    return () => {
      socket.off("job:state_changed", handleJobStateChanged);
    };
  }, [room, refetch]);

  return {
    data: jobs,
    isLoading: isLoading && !hasLoaded, // Only show loading if we haven't loaded yet
    refetch: () => refetch(true), // Manual refetch shows loading
  };
};

export const useJob = (room: string, jobId: string) => {
  return useQuery({
    queryKey: ["job", room, jobId],
    queryFn: async () => {
      try {
        const job = await getJob(room, jobId);
        return job;
      } catch (error) {
        console.error("[useJob] Error fetching job:", error);
        throw error; // Re-throw for individual job errors
      }
    },
    enabled: !!room && !!jobId,
  });
};

export const useExtensionStats = (
  room: string,
  scope: "user" | "global",
  category: string,
  extension: string
) => {
  return useQuery<ExtensionStats>({
    queryKey: ["extensionStats", room, scope, category, extension],
    queryFn: async () => {
      // TODO: Implement getExtensionStats in client.ts
      // For now, return default values
      return {
        idleWorkers: 0,
        busyWorkers: 0,
        pendingJobs: 0,
      };
    },
    enabled: !!room && !!category && !!extension,
    refetchInterval: 10000, // Slower polling for stats (10s)
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