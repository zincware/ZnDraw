import AssignmentIcon from "@mui/icons-material/Assignment";
import CancelIcon from "@mui/icons-material/Cancel";
import CheckCircleIcon from "@mui/icons-material/CheckCircle";
import ErrorIcon from "@mui/icons-material/Error";
import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import PendingIcon from "@mui/icons-material/HourglassEmpty";
import PlayArrowIcon from "@mui/icons-material/PlayArrow";
import {
	Box,
	Chip,
	Collapse,
	IconButton,
	Paper,
	Stack,
	Typography,
} from "@mui/material";
import { useState } from "react";
import type { Task } from "../hooks/useSchemas";
import {
	calculateTaskDurations,
	formatDuration,
	getTaskStatusColor,
	getTaskStatusLabel,
} from "../utils/jobUtils";

interface TaskListItemProps {
	task: Task;
	showExtensionName?: boolean;
}

/**
 * Displays individual task as a slim one-liner with status, time, and ID
 */
export const TaskListItem = ({
	task,
	showExtensionName = false,
}: TaskListItemProps) => {
	const [expanded, setExpanded] = useState(false);
	const durations = calculateTaskDurations(task);
	const statusColor = getTaskStatusColor(task.status);
	const statusLabel = getTaskStatusLabel(task.status);

	const getStatusIcon = () => {
		switch (task.status) {
			case "pending":
				return <PendingIcon fontSize="small" />;
			case "claimed":
				return <AssignmentIcon fontSize="small" />;
			case "running":
				return <PlayArrowIcon fontSize="small" />;
			case "completed":
				return <CheckCircleIcon fontSize="small" />;
			case "failed":
				return <ErrorIcon fontSize="small" />;
			case "cancelled":
				return <CancelIcon fontSize="small" />;
			default:
				return <PendingIcon fontSize="small" />;
		}
	};

	// Get time display: runtime for finished tasks, time since created for others
	const getTimeDisplay = () => {
		if (
			task.status === "completed" ||
			task.status === "failed" ||
			task.status === "cancelled"
		) {
			return durations.total ? formatDuration(durations.total) : "-";
		}

		// For pending/claimed/running: show time since created
		const createdAt = new Date(task.created_at);
		const now = new Date();
		const elapsed = Math.floor((now.getTime() - createdAt.getTime()) / 1000);
		return formatDuration(elapsed * 1000); // formatDuration expects milliseconds
	};

	// Get queue position display for pending tasks
	const getQueuePositionDisplay = () => {
		if (
			task.status === "pending" &&
			task.queue_position !== undefined &&
			task.queue_position !== null
		) {
			return `#${task.queue_position + 1}`; // Display 1-indexed position
		}
		return null;
	};

	const queuePosition = getQueuePositionDisplay();

	return (
		<Paper
			elevation={1}
			sx={{
				p: 1.5,
				cursor: "pointer",
				"&:hover": {
					bgcolor: "action.hover",
				},
			}}
			onClick={() => setExpanded(!expanded)}
		>
			<Stack direction="row" alignItems="center" spacing={1.5}>
				<Chip
					icon={getStatusIcon()}
					label={statusLabel}
					color={statusColor}
					size="small"
					sx={{ minWidth: 100 }}
				/>

				{queuePosition && (
					<Chip
						label={queuePosition}
						size="small"
						variant="outlined"
						sx={{ minWidth: 40 }}
					/>
				)}

				<Typography variant="body2" fontWeight="medium">
					{getTimeDisplay()}
				</Typography>

				<Typography
					variant="body2"
					color="text.secondary"
					sx={{ fontFamily: "monospace", fontSize: "0.85rem" }}
				>
					{task.id.substring(0, 8)}
				</Typography>

				{/* Show worker ID for claimed/running tasks */}
				{task.worker_id &&
					(task.status === "claimed" || task.status === "running") && (
						<Chip
							label={`worker:${task.worker_id.substring(0, 6)}`}
							size="small"
							variant="outlined"
							sx={{ fontFamily: "monospace", fontSize: "0.75rem" }}
						/>
					)}

				<Box sx={{ flex: 1 }} />

				<IconButton
					size="small"
					onClick={(e) => {
						e.stopPropagation();
						setExpanded(!expanded);
					}}
				>
					{expanded ? (
						<ExpandLessIcon fontSize="small" />
					) : (
						<ExpandMoreIcon fontSize="small" />
					)}
				</IconButton>
			</Stack>

			{/* Expanded details */}
			<Collapse in={expanded}>
				<Box sx={{ mt: 2, pt: 2, borderTop: 1, borderColor: "divider" }}>
					<Stack spacing={1}>
						<Box>
							<Typography variant="caption" color="text.secondary">
								Task ID:
							</Typography>
							<Typography variant="body2" fontFamily="monospace">
								{task.id}
							</Typography>
						</Box>

						{showExtensionName && (
							<Box>
								<Typography variant="caption" color="text.secondary">
									Name:
								</Typography>
								<Typography variant="body2">{task.job_name}</Typography>
							</Box>
						)}

						<Box>
							<Typography variant="caption" color="text.secondary">
								Created:
							</Typography>
							<Typography variant="body2">
								{new Date(task.created_at).toLocaleString()}
							</Typography>
						</Box>

						{task.started_at && (
							<Box>
								<Typography variant="caption" color="text.secondary">
									Started:
								</Typography>
								<Typography variant="body2">
									{new Date(task.started_at).toLocaleString()}
								</Typography>
							</Box>
						)}

						{task.completed_at && (
							<Box>
								<Typography variant="caption" color="text.secondary">
									Completed:
								</Typography>
								<Typography variant="body2">
									{new Date(task.completed_at).toLocaleString()}
								</Typography>
							</Box>
						)}

						{task.worker_id && (
							<Box>
								<Typography variant="caption" color="text.secondary">
									Worker:
								</Typography>
								<Typography variant="body2" fontFamily="monospace">
									{task.worker_id.substring(0, 12)}...
								</Typography>
							</Box>
						)}

						{task.error && (
							<Box>
								<Typography variant="caption" color="error">
									Error:
								</Typography>
								<Typography variant="body2" color="error">
									{task.error}
								</Typography>
							</Box>
						)}

						{/* Timing information */}
						{(durations.waiting || durations.running) && (
							<Box sx={{ mt: 1 }}>
								<Typography
									variant="caption"
									color="text.secondary"
									display="block"
									sx={{ mb: 0.5 }}
								>
									Timing:
								</Typography>
								<Stack direction="row" spacing={2}>
									{durations.waiting && (
										<Typography variant="body2">
											<strong>Wait:</strong> {formatDuration(durations.waiting)}
										</Typography>
									)}
									{durations.running && (
										<Typography variant="body2">
											<strong>Run:</strong> {formatDuration(durations.running)}
										</Typography>
									)}
									{durations.total && (
										<Typography variant="body2">
											<strong>Total:</strong> {formatDuration(durations.total)}
										</Typography>
									)}
								</Stack>
							</Box>
						)}
					</Stack>
				</Box>
			</Collapse>
		</Paper>
	);
};
