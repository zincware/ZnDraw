import ExpandLessIcon from "@mui/icons-material/ExpandLess";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import HistoryIcon from "@mui/icons-material/History";
import RefreshIcon from "@mui/icons-material/Refresh";
import {
	Box,
	Button,
	Collapse,
	Divider,
	IconButton,
	Skeleton,
	Stack,
	TablePagination,
	Typography,
} from "@mui/material";
import { useEffect, useState } from "react";
import { useTasks } from "../hooks/useSchemas";
import { TaskListItem } from "./JobListItem";

interface TaskHistoryPanelProps {
	roomId: string;
	jobName: string | null;
}

/**
 * Shows recent task history for the current extension with server-side pagination.
 */
export const TaskHistoryPanel = ({
	roomId,
	jobName,
}: TaskHistoryPanelProps) => {
	const [expanded, setExpanded] = useState(true);
	const [page, setPage] = useState(0);
	const [rowsPerPage, setRowsPerPage] = useState(5);

	const {
		data: tasksResponse,
		isLoading,
		refetch,
	} = useTasks(roomId, jobName, {
		limit: rowsPerPage,
		offset: page * rowsPerPage,
	});

	const tasks = tasksResponse?.items ?? [];
	const totalTasks = tasksResponse?.total ?? 0;

	// Reset page when job changes
	// biome-ignore lint/correctness/useExhaustiveDependencies: reset page on jobName change
	useEffect(() => {
		setPage(0);
	}, [jobName]);

	const handleChangePage = (_event: unknown, newPage: number) => {
		setPage(newPage);
	};

	const handleChangeRowsPerPage = (
		event: React.ChangeEvent<HTMLInputElement>,
	) => {
		setRowsPerPage(Number.parseInt(event.target.value, 10));
		setPage(0);
	};

	if (!jobName) return null;

	return (
		<Box sx={{ mt: 2 }}>
			<Divider sx={{ mb: 2 }} />

			<Stack
				direction="row"
				alignItems="center"
				justifyContent="space-between"
				sx={{ mb: 2 }}
			>
				<Stack direction="row" alignItems="center" spacing={1}>
					<IconButton
						size="small"
						onClick={() => setExpanded(!expanded)}
						aria-label={
							expanded ? "collapse recent tasks" : "expand recent tasks"
						}
					>
						{expanded ? (
							<ExpandLessIcon fontSize="small" />
						) : (
							<ExpandMoreIcon fontSize="small" />
						)}
					</IconButton>
					<HistoryIcon fontSize="small" color="action" />
					<Typography variant="subtitle2" color="text.secondary">
						Recent Tasks
					</Typography>
				</Stack>

				<Button
					size="small"
					startIcon={<RefreshIcon />}
					onClick={() => refetch()}
					disabled={isLoading}
				>
					Refresh
				</Button>
			</Stack>

			<Collapse in={expanded}>
				{isLoading && !tasks.length ? (
					<Stack spacing={1}>
						<Skeleton
							variant="rectangular"
							height={48}
							sx={{ borderRadius: 1 }}
						/>
						<Skeleton
							variant="rectangular"
							height={48}
							sx={{ borderRadius: 1 }}
						/>
						<Skeleton
							variant="rectangular"
							height={48}
							sx={{ borderRadius: 1 }}
						/>
					</Stack>
				) : tasks.length === 0 ? (
					<Box sx={{ textAlign: "center", py: 2 }}>
						<Typography variant="body2" color="text.secondary">
							No tasks yet. Run an extension to see task history.
						</Typography>
					</Box>
				) : (
					<>
						<Stack spacing={1}>
							{tasks.map((task) => (
								<TaskListItem
									key={task.id}
									task={task}
									showExtensionName={false}
								/>
							))}
						</Stack>
						<TablePagination
							component="div"
							count={totalTasks}
							page={page}
							onPageChange={handleChangePage}
							rowsPerPage={rowsPerPage}
							onRowsPerPageChange={handleChangeRowsPerPage}
							rowsPerPageOptions={[5, 10, 25]}
						/>
					</>
				)}
			</Collapse>
		</Box>
	);
};
