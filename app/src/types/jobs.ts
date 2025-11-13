/**
 * Job status types matching backend JobStatus enum
 */
export type JobStatus = "pending" | "assigned" | "processing" | "completed" | "failed";

/**
 * Job state change event data received via socket
 */
export interface JobStateChangeEvent {
  roomId: string;
  jobId: string;
  status: JobStatus;
  metadata: {
    assignedAt?: string;
    startedAt?: string;
    completedAt?: string;
    error?: string;
    workerId?: string;
  };
}

/**
 * Duration spent in each job state (in milliseconds)
 */
export interface JobStateDuration {
  pending?: number;
  assigned?: number;
  processing?: number;
  total?: number;
}

/**
 * Extension statistics showing worker availability and job counts
 */
export interface ExtensionStats {
  idleWorkers: number;
  busyWorkers: number;
  pendingJobs: number;
}
