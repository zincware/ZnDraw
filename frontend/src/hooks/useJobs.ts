import { useQuery, useQueryClient } from "@tanstack/react-query";
import { useEffect, useMemo } from "react";
import { getJobSchema, listJobs } from "../myapi/client";
import { socket } from "../socket";

type DisplayScope = "public" | "internal" | "private";

export interface ParsedJobName {
	scope: string;
	displayScope: DisplayScope;
	category: string;
	name: string;
}

/** Parse job full_name into scope, displayScope, category, and name components. */
export const parseJobName = (fullName: string): ParsedJobName | null => {
	const parts = fullName.split(":");
	if (parts.length !== 3) return null;

	const scope = parts[0];
	let displayScope: DisplayScope;
	if (scope === "@global") {
		displayScope = "public";
	} else if (scope === "@internal") {
		displayScope = "internal";
	} else {
		displayScope = "private";
	}

	return {
		scope,
		displayScope,
		category: parts[1],
		name: parts[2],
	};
};

/**
 * Extract default values from a JSON schema's properties.
 * Replaces the old job_defaults field.
 */
export const extractDefaults = (
	schema: Record<string, unknown>,
): Record<string, unknown> => {
	const properties = schema.properties as
		| Record<string, Record<string, unknown>>
		| undefined;
	if (!properties) return {};

	const defaults: Record<string, unknown> = {};
	for (const [key, prop] of Object.entries(properties)) {
		if (prop.default !== undefined) {
			defaults[key] = prop.default;
		}
	}
	return defaults;
};

/** Fetch all jobs for a room, with socket invalidation on jobs_invalidate. */
export const useJobs = (roomId: string) => {
	const queryClient = useQueryClient();

	const query = useQuery({
		queryKey: ["jobs", roomId],
		queryFn: () => listJobs(roomId),
		enabled: !!roomId,
		staleTime: 1000 * 60 * 5, // 5 minutes
	});

	// Invalidate on jobs_invalidate socket event
	useEffect(() => {
		if (!roomId) return;

		const handleJobsInvalidate = () => {
			queryClient.invalidateQueries({ queryKey: ["jobs", roomId] });
			queryClient.invalidateQueries({ queryKey: ["jobSchema", roomId] });
			queryClient.invalidateQueries({ queryKey: ["tasks", roomId] });
		};

		socket.on("jobs_invalidate", handleJobsInvalidate);
		return () => {
			socket.off("jobs_invalidate", handleJobsInvalidate);
		};
	}, [roomId, queryClient]);

	return query;
};

/** Filter jobs by category (modifiers, analysis, selections, etc.). */
export const useJobsByCategory = (roomId: string, category: string) => {
	const { data: jobsData, isLoading, isError } = useJobs(roomId);

	const filteredJobs = useMemo(() => {
		if (!jobsData) return [];

		return jobsData.items.filter((job) => {
			const parsed = parseJobName(job.full_name);
			return parsed?.category === category;
		});
	}, [jobsData, category]);

	return { data: filteredJobs, isLoading, isError };
};

/** Fetch schema and defaults for a specific job. */
export const useJobSchema = (roomId: string, jobName: string | null) => {
	return useQuery({
		queryKey: ["jobSchema", roomId, jobName],
		queryFn: async () => {
			if (!jobName) throw new Error("jobName is required");
			return getJobSchema(roomId, jobName);
		},
		enabled: !!roomId && !!jobName,
		staleTime: 1000 * 60 * 5, // 5 minutes
	});
};
