/**
 * Task status types matching backend TaskStatus enum values.
 */
export type TaskStatus =
	| "pending"
	| "claimed"
	| "running"
	| "completed"
	| "failed"
	| "cancelled";

/**
 * Task status event data received via socket.
 * Broadcast when task state changes (task_status_event).
 */
export interface TaskStatusEvent {
	id: string;
	name: string;
	room_id: string;
	status: TaskStatus;
	created_at: string;
	started_at: string | null;
	completed_at: string | null;
	queue_position: number | null;
	worker_id: string | null;
	error: string | null;
}

/**
 * Duration spent in each task state (in milliseconds)
 */
export interface TaskStateDuration {
	pending?: number;
	claimed?: number;
	running?: number;
	total?: number;
	/** Combined pending + claimed time (waiting for work to start) */
	waiting?: number;
}

/**
 * Extension statistics showing worker availability and job counts
 */
export interface ExtensionStats {
	idleWorkers: number;
	busyWorkers: number;
	pendingJobs: number;
}
