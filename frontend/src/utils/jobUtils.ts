import type { Task } from "../hooks/useSchemas";
import type { TaskStateDuration, TaskStatus } from "../types/jobs";

/**
 * Calculate duration spent in each task state
 * @param task Task object with timestamps
 * @returns TaskStateDuration object with durations in milliseconds
 */
export function calculateTaskDurations(task: Task): TaskStateDuration {
	const durations: TaskStateDuration = {};

	const createdAt = new Date(task.created_at).getTime();
	const startedAt = task.started_at
		? new Date(task.started_at).getTime()
		: null;
	const completedAt = task.completed_at
		? new Date(task.completed_at).getTime()
		: null;
	const now = Date.now();

	// Calculate pending + claimed duration (created -> started)
	// In the new model, we don't have a separate claimed_at timestamp,
	// so waiting = time from created to started
	if (startedAt) {
		durations.waiting = startedAt - createdAt;
	} else if (task.status === "pending" || task.status === "claimed") {
		durations.waiting = now - createdAt;
	}

	// Calculate running duration (started -> completed)
	if (startedAt) {
		if (completedAt) {
			durations.running = completedAt - startedAt;
		} else if (task.status === "running") {
			durations.running = now - startedAt;
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
	}
	if (minutes > 0) {
		const remainingSeconds = seconds % 60;
		return `${minutes}m ${remainingSeconds}s`;
	}
	if (seconds > 0) {
		return `${seconds}s`;
	}
	return `${Math.round(milliseconds)}ms`;
}

/**
 * Get a color for a task status
 * @param status Task status
 * @returns MUI color string
 */
export function getTaskStatusColor(
	status: TaskStatus,
): "default" | "warning" | "info" | "success" | "error" {
	switch (status) {
		case "pending":
			return "default";
		case "claimed":
			return "warning";
		case "running":
			return "info";
		case "completed":
			return "success";
		case "failed":
			return "error";
		case "cancelled":
			return "default";
		default:
			return "default";
	}
}

/**
 * Get a human-readable label for a task status
 * @param status Task status
 * @returns Label string
 */
export function getTaskStatusLabel(status: TaskStatus): string {
	switch (status) {
		case "pending":
			return "Pending";
		case "claimed":
			return "Claimed";
		case "running":
			return "Running";
		case "completed":
			return "Completed";
		case "failed":
			return "Failed";
		case "cancelled":
			return "Cancelled";
		default:
			return status;
	}
}

/**
 * Calculate progress percentage for a task (0-100)
 * @param task Task object
 * @returns Progress percentage (0-100)
 */
export function calculateTaskProgress(task: Task): number {
	switch (task.status) {
		case "pending":
			return 0;
		case "claimed":
			return 25;
		case "running":
			return 50;
		case "completed":
			return 100;
		case "failed":
			return 100;
		case "cancelled":
			return 100;
		default:
			return 0;
	}
}
