import { type Job } from "../hooks/useSchemas";
import { type JobStateDuration, type JobStatus } from "../types/jobs";

/**
 * Calculate duration spent in each job state
 * @param job Job object with timestamps
 * @returns JobStateDuration object with durations in milliseconds
 */
export function calculateJobDurations(job: Job): JobStateDuration {
  const durations: JobStateDuration = {};

  const createdAt = new Date(job.created_at).getTime();
  const assignedAt = job.assigned_at ? new Date(job.assigned_at).getTime() : null;
  const startedAt = job.started_at ? new Date(job.started_at).getTime() : null;
  const completedAt = job.completed_at ? new Date(job.completed_at).getTime() : null;
  const now = Date.now();

  // Calculate pending duration (created -> assigned)
  if (assignedAt) {
    durations.pending = assignedAt - createdAt;
  } else if (job.status === "pending") {
    durations.pending = now - createdAt;
  }

  // Calculate assigned duration (assigned -> started)
  if (assignedAt) {
    if (startedAt) {
      durations.assigned = startedAt - assignedAt;
    } else if (job.status === "assigned") {
      durations.assigned = now - assignedAt;
    }
  }

  // Calculate processing duration (started -> completed)
  if (startedAt) {
    if (completedAt) {
      durations.processing = completedAt - startedAt;
    } else if (job.status === "processing") {
      durations.processing = now - startedAt;
    }
  }

  // Calculate total duration
  if (completedAt) {
    durations.total = completedAt - createdAt;
  } else {
    durations.total = now - createdAt;
  }

  return durations;
}

/**
 * Format duration in a human-readable format
 * @param milliseconds Duration in milliseconds
 * @returns Formatted string (e.g., "2m 30s", "45s", "1.5s")
 */
export function formatDuration(milliseconds: number | undefined): string {
  if (milliseconds === undefined || milliseconds === null) {
    return "—";
  }

  const seconds = Math.floor(milliseconds / 1000);
  const minutes = Math.floor(seconds / 60);
  const hours = Math.floor(minutes / 60);

  if (hours > 0) {
    const remainingMinutes = minutes % 60;
    return `${hours}h ${remainingMinutes}m`;
  } else if (minutes > 0) {
    const remainingSeconds = seconds % 60;
    return `${minutes}m ${remainingSeconds}s`;
  } else if (seconds > 0) {
    return `${seconds}s`;
  } else {
    return `${Math.round(milliseconds)}ms`;
  }
}

/**
 * Get a color for a job status
 * @param status Job status
 * @returns MUI color string
 */
export function getJobStatusColor(status: JobStatus): "default" | "warning" | "info" | "success" | "error" {
  switch (status) {
    case "pending":
      return "default";
    case "assigned":
      return "warning";
    case "processing":
      return "info";
    case "completed":
      return "success";
    case "failed":
      return "error";
    default:
      return "default";
  }
}

/**
 * Get a human-readable label for a job status
 * @param status Job status
 * @returns Label string
 */
export function getJobStatusLabel(status: JobStatus): string {
  switch (status) {
    case "pending":
      return "Pending";
    case "assigned":
      return "Assigned";
    case "processing":
      return "Processing";
    case "completed":
      return "Completed";
    case "failed":
      return "Failed";
    default:
      return status;
  }
}

/**
 * Calculate progress percentage for a job (0-100)
 * @param job Job object
 * @returns Progress percentage (0-100)
 */
export function calculateJobProgress(job: Job): number {
  switch (job.status) {
    case "pending":
      return 0;
    case "assigned":
      return 25;
    case "processing":
      return 50;
    case "completed":
      return 100;
    case "failed":
      return 100;
    default:
      return 0;
  }
}
