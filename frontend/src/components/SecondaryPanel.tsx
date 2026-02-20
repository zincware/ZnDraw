import { materialCells } from "@jsonforms/material-renderers";
import { JsonForms } from "@jsonforms/react";
import SaveIcon from "@mui/icons-material/Save";
import {
	Box,
	Button,
	Chip,
	Divider,
	Fade,
	FormControl,
	InputLabel,
	MenuItem,
	Select,
	Tooltip,
	Typography,
} from "@mui/material";
import type { SelectChangeEvent } from "@mui/material";
import { memo, useCallback, useEffect, useMemo, useState } from "react";
import { useFormStore } from "../formStore";
import {
	extractDefaults,
	parseJobName,
	useJobSchema,
	useJobsByCategory,
} from "../hooks/useJobs";
import { useTasks } from "../hooks/useSchemas";
import { submitTask } from "../myapi/client";
import { selectIsRoomReadOnly, useAppStore } from "../store";
import { loadFormData, saveFormData } from "../utils/formStorage";
import { customRenderers, injectDynamicEnums } from "../utils/jsonforms";
import { TaskHistoryPanel } from "./JobHistoryPanel";
import { JobStatusChips } from "./JobStatusChips";
import { FormSkeleton, PanelSkeleton } from "./shared/LoadingSkeletons";

interface SecondaryPanelProps {
	panelTitle: string;
}

const SecondaryPanel = ({ panelTitle }: SecondaryPanelProps) => {
	const roomId = useAppStore((state) => state.roomId);
	const userName = useAppStore((state) => state.user?.email ?? null);
	const geometries = useAppStore((state) => state.geometries);
	const roomReadOnly = useAppStore(selectIsRoomReadOnly);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const [localFormData, setLocalFormData] = useState<Record<string, unknown>>(
		{},
	);
	const [isSubmitting, setIsSubmitting] = useState(false);

	const { selectedExtensions, setSelectedExtension } = useFormStore();
	const selectedJobName = selectedExtensions[panelTitle] || null;

	const {
		data: jobs,
		isLoading: isLoadingJobs,
		isError: isJobsError,
	} = useJobsByCategory(roomId ?? "", panelTitle);

	const {
		data: schemaResponse,
		isLoading: isLoadingSchema,
		isError: isSchemaError,
	} = useJobSchema(roomId ?? "", selectedJobName);

	const { data: pendingResponse } = useTasks(roomId ?? "", selectedJobName, {
		status: "pending",
		limit: 0,
	});
	const pendingTaskCount = pendingResponse?.total ?? 0;

	const selectedJob = useMemo(() => {
		if (!jobs || !selectedJobName) return null;
		return jobs.find((j) => j.full_name === selectedJobName);
	}, [jobs, selectedJobName]);

	const currentSchema = schemaResponse?.schema;

	useEffect(() => {
		if (!schemaResponse || !selectedJobName || !userName || !roomId) return;

		const storedData = loadFormData(roomId, userName, selectedJobName);
		const mergedData = {
			...extractDefaults(schemaResponse.schema),
			...storedData,
		};
		setLocalFormData(mergedData);
	}, [schemaResponse, roomId, userName, selectedJobName]);

	const handleFormChange = useCallback(
		({ data }: { data: Record<string, unknown> }) => {
			const safeData = data ?? {};
			setLocalFormData(safeData);

			if (selectedJobName && userName && roomId) {
				saveFormData(roomId, userName, selectedJobName, safeData);
			}
		},
		[roomId, userName, selectedJobName],
	);

	const handleSelectionChange = (event: SelectChangeEvent<string>) => {
		setSelectedExtension(panelTitle, event.target.value || null);
	};

	const handleSubmit = async () => {
		if (!selectedJobName || !roomId || !userName) return;

		setIsSubmitting(true);
		try {
			await submitTask(roomId, selectedJobName, localFormData);
			showSnackbar("Task submitted successfully", "success");
		} catch (error) {
			console.error("Failed to submit task:", error);
			showSnackbar("Failed to submit task", "error");
		} finally {
			setIsSubmitting(false);
		}
	};

	// dynamic-atom-props enums are fetched by DynamicEnumRenderer itself;
	// injectDynamicEnums only handles dynamic-geometries here.
	const dynamicSchema = useMemo(() => {
		if (!currentSchema) return null;
		return injectDynamicEnums(currentSchema as Record<string, unknown>, geometries);
	}, [currentSchema, geometries]);

	if (!roomId || !userName) {
		return <Typography sx={{ p: 2 }}>Joining room...</Typography>;
	}

	if (isJobsError || isSchemaError) {
		return (
			<Typography color="error" sx={{ p: 2 }}>
				Failed to load configuration.
			</Typography>
		);
	}

	return (
		<Box sx={{ display: "flex", flexDirection: "column", height: "100%" }}>
			<Typography variant="h6" sx={{ p: 2, pb: 1, flexShrink: 0 }}>
				{panelTitle}
			</Typography>
			<Divider />

			<Box sx={{ p: 2, pb: 12, flexGrow: 1, overflowY: "auto" }}>
				{isLoadingJobs ? (
					<PanelSkeleton />
				) : (
					<>
						<FormControl fullWidth sx={{ mb: 2 }}>
							<InputLabel id="panel-select-label">
								{panelTitle} Method
							</InputLabel>
							<Select
								labelId="panel-select-label"
								value={selectedJobName || ""}
								label={`${panelTitle} Method`}
								onChange={handleSelectionChange}
							>
								{jobs?.map((job) => {
									const parsed = parseJobName(job.full_name);
									return (
										<MenuItem key={job.full_name} value={job.full_name}>
											<Box
												sx={{ display: "flex", alignItems: "center", gap: 1 }}
											>
												<Typography>{parsed?.name}</Typography>
												<Chip
													label={parsed?.displayScope}
													size="small"
													color={
														parsed?.displayScope === "internal"
															? "primary"
															: parsed?.displayScope === "public"
																? "success"
																: "default"
													}
												/>
											</Box>
										</MenuItem>
									);
								})}
							</Select>
						</FormControl>

						{selectedJob && (
							<JobStatusChips
								jobName={selectedJob.full_name}
								workerCount={selectedJob.workers.length}
								pendingTaskCount={pendingTaskCount}
							/>
						)}

						{selectedJobName && (
							<>
								{panelTitle !== "settings" && (
									<Tooltip title={roomReadOnly ? "Room is locked" : ""}>
										<span>
											<Button
												variant="contained"
												startIcon={<SaveIcon />}
												onClick={handleSubmit}
												disabled={
													isSubmitting || isLoadingSchema || roomReadOnly
												}
												fullWidth
												color="primary"
												sx={{ mb: 2 }}
											>
												{isSubmitting ? "Running..." : "Run Extension"}
											</Button>
										</span>
									</Tooltip>
								)}

								{isLoadingSchema ? (
									<FormSkeleton />
								) : dynamicSchema ? (
									<Fade in={!isLoadingSchema} timeout={200}>
										<Box>
											<JsonForms
												key={selectedJobName}
												schema={dynamicSchema}
												data={localFormData}
												renderers={customRenderers}
												cells={materialCells}
												onChange={handleFormChange}
											/>
										</Box>
									</Fade>
								) : null}

								{panelTitle !== "settings" && (
									<TaskHistoryPanel roomId={roomId} jobName={selectedJobName} />
								)}
							</>
						)}
					</>
				)}
			</Box>
		</Box>
	);
};

export default memo(SecondaryPanel);
