import {
	keepPreviousData,
	useQuery,
	useQueryClient,
} from "@tanstack/react-query";
import { useEffect } from "react";
import { type Task, getFrameMetadata, listTasksForJob } from "../myapi/client";
import { socket } from "../socket";
import type { TaskStatusEvent } from "../types/jobs";

export type { Task };

export const useTasks = (
	roomId: string,
	jobName: string | null,
	options?: {
		status?: import("../types/jobs").TaskStatus;
		limit?: number;
		offset?: number;
	},
) => {
	const queryClient = useQueryClient();
	const status = options?.status;
	const limit = options?.limit ?? 10;
	const offset = options?.offset ?? 0;

	const query = useQuery({
		queryKey: ["tasks", roomId, jobName, { status, limit, offset }],
		queryFn: () =>
			listTasksForJob(roomId, jobName as string, {
				status,
				limit,
				offset,
			}),
		enabled: !!roomId && !!jobName,
		staleTime: 1000 * 30,
		placeholderData: keepPreviousData,
	});

	// Invalidate on task_status_event for this job
	useEffect(() => {
		if (!roomId || !jobName) return;

		const handleTaskStatusChanged = (data: TaskStatusEvent) => {
			if (data.room_id !== roomId) return;
			if (data.name !== jobName) return;
			queryClient.invalidateQueries({
				queryKey: ["tasks", roomId, jobName],
			});
		};

		socket.on("task_status_event", handleTaskStatusChanged);
		return () => {
			socket.off("task_status_event", handleTaskStatusChanged);
		};
	}, [roomId, jobName, queryClient]);

	return query;
};

export const useFrameMetadata = (
	roomId: string,
	frameId = 0,
	enabled = true,
) => {
	return useQuery({
		queryKey: ["metadata", roomId, frameId],
		queryFn: () => getFrameMetadata(roomId, frameId),
		enabled: !!roomId && enabled,
		staleTime: 1000 * 60 * 5,
	});
};
