import { materialCells } from "@jsonforms/material-renderers";
import { JsonForms } from "@jsonforms/react";
import ClearIcon from "@mui/icons-material/Clear";
import DeleteIcon from "@mui/icons-material/Delete";
import EditIcon from "@mui/icons-material/Edit";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import RadioButtonCheckedIcon from "@mui/icons-material/RadioButtonChecked";
import RadioButtonUncheckedIcon from "@mui/icons-material/RadioButtonUnchecked";
import SaveIcon from "@mui/icons-material/Save";
import {
	Accordion,
	AccordionDetails,
	AccordionSummary,
	Box,
	Button,
	Chip,
	Dialog,
	DialogActions,
	DialogContent,
	DialogContentText,
	DialogTitle,
	Divider,
	Fade,
	FormControl,
	IconButton,
	InputLabel,
	MenuItem,
	Select,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import type { SelectChangeEvent } from "@mui/material";
import {
	DataGrid,
	type GridColDef,
	type GridRenderCellParams,
	type GridRowParams,
} from "@mui/x-data-grid";
import { useCallback, useEffect, useMemo, useState } from "react";
import { useFormStore } from "../formStore";
import {
	extractDefaults,
	parseJobName,
	useJobSchema,
	useJobsByCategory,
} from "../hooks/useJobs";
import {
	createUpdateSelectionGroup,
	deleteSelectionGroup,
	updateSelection,
	submitTask,
} from "../myapi/client";
import {
	selectActiveSelectionGroup,
	selectIsRoomReadOnly,
	useAppStore,
} from "../store";
import { loadFormData, saveFormData } from "../utils/formStorage";
import { customRenderers, injectDynamicEnums } from "../utils/jsonforms";
import { TaskHistoryPanel } from "./JobHistoryPanel";
import { JobStatusChips } from "./JobStatusChips";
import {
	FormSkeleton,
	SelectionToolsSkeleton,
} from "./shared/LoadingSkeletons";

interface SelectionGroupRow {
	id: string;
	name: string;
	isCurrent: boolean;
	isActive: boolean;
	selections: Record<string, number>; // geometry -> count
	rawData: Record<string, number[]>; // geometry -> indices array
}

export default function SelectionsPanel() {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);
	const userName = useAppStore((state) => state.user?.email ?? null);
	const selections = useAppStore((state) => state.selections);
	const selectionGroups = useAppStore((state) => state.selectionGroups);
	const activeSelectionGroup = useAppStore(selectActiveSelectionGroup);
	const updateSelectionForGeometry = useAppStore(
		(state) => state.updateSelectionForGeometry,
	);
	const geometries = useAppStore((state) => state.geometries);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const roomReadOnly = useAppStore(selectIsRoomReadOnly);

	// Selection Groups state
	const [currentGroupName, setCurrentGroupName] = useState("");
	const [editingGroupName, setEditingGroupName] = useState<string | null>(null);
	const [editNameValue, setEditNameValue] = useState("");
	const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
	const [groupToDelete, setGroupToDelete] = useState<string | null>(null);
	const [isLoading, setIsLoading] = useState(false);

	// Selection Tools state
	const [localFormData, setLocalFormData] = useState<Record<string, unknown>>(
		{},
	);
	const [isSubmitting, setIsSubmitting] = useState(false);
	const panelTitle = "selections";

	const { selectedExtensions, setSelectedExtension } = useFormStore();
	const selectedJobName = selectedExtensions[panelTitle] || null;

	// Fetch jobs for this category
	const {
		data: jobs,
		isLoading: isLoadingJobs,
		isError: isJobsError,
	} = useJobsByCategory(roomId!, panelTitle);

	// Fetch schema and defaults for selected job
	const {
		data: schemaResponse,
		isLoading: isLoadingSchema,
		isError: isSchemaError,
	} = useJobSchema(roomId!, selectedJobName);

	// Get the selected job summary for worker count
	const selectedJob = useMemo(() => {
		if (!jobs || !selectedJobName) return null;
		return jobs.find((j) => j.full_name === selectedJobName);
	}, [jobs, selectedJobName]);

	const currentSchema = schemaResponse?.schema;

	// Initialize form data from defaults + localStorage
	useEffect(() => {
		if (!schemaResponse || !selectedJobName || !userName || !roomId) return;

		const storedData = loadFormData(roomId, userName, selectedJobName);
		const mergedData = {
			...extractDefaults(schemaResponse.schema),
			...storedData,
		};
		setLocalFormData(mergedData);
	}, [schemaResponse, roomId, userName, selectedJobName]);

	// Handle form changes - update local state and persist to localStorage
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

	const dynamicSchema = useMemo(() => {
		if (!currentSchema) return null;
		return injectDynamicEnums(
			currentSchema as Record<string, unknown>,
			geometries,
		);
	}, [currentSchema, geometries]);

	// Create form options from jobs list
	const formOptions = useMemo(() => {
		if (!jobs || !Array.isArray(jobs)) return [];

		return jobs.map((job) => {
			const parsed = parseJobName(job.full_name);
			return {
				jobName: job.full_name,
				displayName: parsed?.name || job.full_name,
				scope: parsed?.displayScope || "private",
			};
		});
	}, [jobs]);

	// Selection Groups logic
	const allGeometryKeys = useMemo(() => {
		const keys = new Set<string>();

		Object.keys(selections).forEach((key) => {
			if (selections[key] && selections[key].length > 0) {
				keys.add(key);
			}
		});

		Object.values(selectionGroups).forEach((group) => {
			Object.keys(group).forEach((key) => {
				if (group[key] && group[key].length > 0) {
					keys.add(key);
				}
			});
		});

		return Array.from(keys).sort();
	}, [selections, selectionGroups]);

	const rows = useMemo((): SelectionGroupRow[] => {
		const result: SelectionGroupRow[] = [];

		const currentSelectionCounts: Record<string, number> = {};
		allGeometryKeys.forEach((key) => {
			currentSelectionCounts[key] = selections[key]?.length || 0;
		});

		result.push({
			id: "current",
			name: "",
			isCurrent: true,
			isActive: false,
			selections: currentSelectionCounts,
			rawData: selections,
		});

		const groupNames = Object.keys(selectionGroups).sort();
		groupNames.forEach((groupName) => {
			const groupData = selectionGroups[groupName];
			const groupCounts: Record<string, number> = {};

			allGeometryKeys.forEach((key) => {
				groupCounts[key] = groupData[key]?.length || 0;
			});

			result.push({
				id: groupName,
				name: groupName,
				isCurrent: false,
				isActive: activeSelectionGroup === groupName,
				selections: groupCounts,
				rawData: groupData,
			});
		});

		return result;
	}, [selections, selectionGroups, activeSelectionGroup, allGeometryKeys]);

	const hasCurrentSelection = useMemo(
		() =>
			Object.values(selections).some(
				(indices) => indices && indices.length > 0,
			),
		[selections],
	);

	const handleSaveCurrentSelection = useCallback(async () => {
		if (!roomId || !currentGroupName.trim()) return;

		if (selectionGroups[currentGroupName.trim()]) {
			showSnackbar(
				`Group "${currentGroupName.trim()}" already exists`,
				"error",
			);
			return;
		}

		setIsLoading(true);
		try {
			await createUpdateSelectionGroup(
				roomId,
				currentGroupName.trim(),
				selections,
			);
			setCurrentGroupName("");
			showSnackbar(`Group "${currentGroupName.trim()}" saved`, "success");
		} catch (err) {
			console.error("Failed to save selection group:", err);
			showSnackbar("Failed to save selection group", "error");
		} finally {
			setIsLoading(false);
		}
	}, [roomId, currentGroupName, selections, selectionGroups, showSnackbar]);

	const handleSaveNameEdit = useCallback(
		async (oldName: string) => {
			if (
				!roomId ||
				!editNameValue.trim() ||
				editNameValue.trim() === oldName
			) {
				setEditingGroupName(null);
				return;
			}

			if (
				selectionGroups[editNameValue.trim()] &&
				editNameValue.trim() !== oldName
			) {
				showSnackbar(`Group "${editNameValue.trim()}" already exists`, "error");
				setEditingGroupName(null);
				return;
			}

			setIsLoading(true);
			try {
				await deleteSelectionGroup(roomId, oldName);
				await createUpdateSelectionGroup(
					roomId,
					editNameValue.trim(),
					selectionGroups[oldName],
				);
				setEditingGroupName(null);
				showSnackbar(`Group renamed to "${editNameValue.trim()}"`, "success");
			} catch (err) {
				console.error("Failed to rename group:", err);
				showSnackbar("Failed to rename group", "error");
				setEditingGroupName(null);
			} finally {
				setIsLoading(false);
			}
		},
		[roomId, editNameValue, selectionGroups, showSnackbar],
	);

	const handleRowClick = useCallback(
		async (params: GridRowParams<SelectionGroupRow>) => {
			const row = params.row;

			if (row.isCurrent || editingGroupName) return;

			if (!roomId) return;

			setIsLoading(true);
			try {
				// Apply each geometry's selection via REST API
				await Promise.all(
					Object.entries(row.rawData).map(([geoKey, indices]) =>
						updateSelection(roomId, geoKey, indices),
					),
				);
				// Success - no notification needed
			} catch (err) {
				console.error("Failed to load selection group:", err);
				showSnackbar(`Failed to load group "${row.name}"`, "error");
			} finally {
				setIsLoading(false);
			}
		},
		[roomId, editingGroupName, showSnackbar],
	);

	const handleClearCurrent = useCallback(() => {
		if (!hasCurrentSelection) return;

		allGeometryKeys.forEach((key) => {
			if (selections[key] && selections[key].length > 0) {
				updateSelectionForGeometry(key, []);
			}
		});
	}, [
		hasCurrentSelection,
		allGeometryKeys,
		selections,
		updateSelectionForGeometry,
	]);

	const handleDeleteConfirm = useCallback(async () => {
		if (!roomId || !groupToDelete) return;

		setIsLoading(true);
		setDeleteDialogOpen(false);

		try {
			await deleteSelectionGroup(roomId, groupToDelete);
			showSnackbar(`Group "${groupToDelete}" deleted`, "success");
		} catch (err) {
			console.error("Failed to delete group:", err);
			showSnackbar(`Failed to delete group "${groupToDelete}"`, "error");
		} finally {
			setIsLoading(false);
			setGroupToDelete(null);
		}
	}, [roomId, groupToDelete, showSnackbar]);

	const columns = useMemo((): GridColDef[] => {
		const cols: GridColDef[] = [
			{
				field: "name",
				headerName: "Name",
				width: 200,
				flex: 1,
				renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
					const row = params.row;

					if (row.isCurrent) {
						return (
							<Box
								sx={{
									display: "flex",
									alignItems: "center",
									gap: 1,
									width: "100%",
									height: "100%",
								}}
							>
								<TextField
									size="small"
									fullWidth
									placeholder="Type name to save..."
									value={currentGroupName}
									onChange={(e) => setCurrentGroupName(e.target.value)}
									onKeyDown={(e) => {
										if (e.key === "Enter" && currentGroupName.trim()) {
											handleSaveCurrentSelection();
										}
									}}
									onBlur={() => {
										if (currentGroupName.trim()) {
											handleSaveCurrentSelection();
										}
									}}
									disabled={isLoading || !hasCurrentSelection}
									sx={{ "& .MuiInputBase-root": { fontSize: "0.875rem" } }}
								/>
							</Box>
						);
					}

					const isEditing = editingGroupName === row.id;

					return (
						<Box
							sx={{
								display: "flex",
								alignItems: "center",
								gap: 1,
								width: "100%",
								height: "100%",
							}}
						>
							{row.isActive ? (
								<RadioButtonCheckedIcon
									sx={{ color: "success.main" }}
									fontSize="small"
								/>
							) : (
								<RadioButtonUncheckedIcon
									sx={{ color: "action.disabled" }}
									fontSize="small"
								/>
							)}
							{isEditing ? (
								<TextField
									size="small"
									fullWidth
									autoFocus
									value={editNameValue}
									onChange={(e) => setEditNameValue(e.target.value)}
									onKeyDown={(e) => {
										if (e.key === "Enter") {
											handleSaveNameEdit(row.id);
										} else if (e.key === "Escape") {
											setEditingGroupName(null);
										}
									}}
									onBlur={() => handleSaveNameEdit(row.id)}
									sx={{ "& .MuiInputBase-root": { fontSize: "0.875rem" } }}
								/>
							) : (
								<Typography
									sx={{
										fontWeight: row.isActive ? "bold" : "normal",
										fontSize: "0.875rem",
										flex: 1,
									}}
									onDoubleClick={() => {
										setEditingGroupName(row.id);
										setEditNameValue(row.name);
									}}
								>
									{row.name}
								</Typography>
							)}
						</Box>
					);
				},
			},
		];

		allGeometryKeys.forEach((geometryKey) => {
			cols.push({
				field: geometryKey,
				headerName: geometryKey,
				width: 100,
				align: "center",
				headerAlign: "center",
				renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
					const row = params.row;
					const count = row.selections[geometryKey] || 0;

					if (count === 0) {
						return (
							<Box
								sx={{
									display: "flex",
									alignItems: "center",
									justifyContent: "center",
									height: "100%",
								}}
							>
								<Typography
									sx={{ color: "text.disabled", fontSize: "0.875rem" }}
								>
									-
								</Typography>
							</Box>
						);
					}

					let chipColor: "primary" | "success" | "default" = "default";
					if (row.isCurrent) {
						chipColor = "primary";
					} else if (row.isActive) {
						chipColor = "success";
					}

					return (
						<Box
							sx={{
								display: "flex",
								alignItems: "center",
								justifyContent: "center",
								height: "100%",
							}}
						>
							<Chip label={count} size="small" color={chipColor} />
						</Box>
					);
				},
			});
		});

		cols.push({
			field: "actions",
			headerName: "Actions",
			width: 100,
			sortable: false,
			renderCell: (params: GridRenderCellParams<SelectionGroupRow>) => {
				const row = params.row;

				if (row.isCurrent) {
					return (
						<Box
							sx={{
								display: "flex",
								alignItems: "center",
								justifyContent: "center",
								height: "100%",
							}}
						>
							<Tooltip title="Clear all current selections">
								<span>
									<IconButton
										size="small"
										onClick={(e) => {
											e.stopPropagation();
											handleClearCurrent();
										}}
										disabled={!hasCurrentSelection || isLoading}
									>
										<ClearIcon fontSize="small" />
									</IconButton>
								</span>
							</Tooltip>
						</Box>
					);
				}

				return (
					<Box
						sx={{
							display: "flex",
							alignItems: "center",
							justifyContent: "center",
							gap: 0.5,
							height: "100%",
						}}
					>
						<Tooltip title="Edit name">
							<IconButton
								size="small"
								onClick={(e) => {
									e.stopPropagation();
									setEditingGroupName(row.id);
									setEditNameValue(row.name);
								}}
								disabled={isLoading}
							>
								<EditIcon fontSize="small" />
							</IconButton>
						</Tooltip>
						<Tooltip title="Delete group">
							<IconButton
								size="small"
								onClick={(e) => {
									e.stopPropagation();
									setGroupToDelete(row.id);
									setDeleteDialogOpen(true);
								}}
								disabled={isLoading}
							>
								<DeleteIcon fontSize="small" />
							</IconButton>
						</Tooltip>
					</Box>
				);
			},
		});

		return cols;
	}, [
		allGeometryKeys,
		currentGroupName,
		editingGroupName,
		editNameValue,
		isLoading,
		hasCurrentSelection,
		handleSaveCurrentSelection,
		handleSaveNameEdit,
		handleClearCurrent,
	]);

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
		<Box
			sx={{
				height: "100%",
				width: "100%",
				display: "flex",
				flexDirection: "column",
			}}
		>
			<Typography variant="h6" sx={{ p: 2, pb: 1 }}>
				Selections
			</Typography>
			<Divider />

			<Box
				sx={{
					flexGrow: 1,
					overflowY: "auto",
					overscrollBehavior: "contain",
					p: 2,
				}}
			>
				{isLoadingJobs ? (
					<SelectionToolsSkeleton />
				) : (
					<>
						{/* Selection Tools Section */}
						<Accordion defaultExpanded>
							<AccordionSummary expandIcon={<ExpandMoreIcon />}>
								<Typography variant="subtitle1" fontWeight="medium">
									Selection Tools
								</Typography>
							</AccordionSummary>
							<AccordionDetails>
								<FormControl fullWidth sx={{ mb: 2 }}>
									<InputLabel id="panel-select-label">
										Selection Method
									</InputLabel>
									<Select
										labelId="panel-select-label"
										value={selectedJobName || ""}
										label="Selection Method"
										onChange={handleSelectionChange}
									>
										{formOptions.map((option) => (
											<MenuItem key={option.jobName} value={option.jobName}>
												<Box
													sx={{ display: "flex", alignItems: "center", gap: 1 }}
												>
													<Typography>{option.displayName}</Typography>
													<Chip
														label={option.scope}
														size="small"
														color={
															option.scope === "internal"
																? "primary"
																: option.scope === "public"
																	? "success"
																	: "default"
														}
													/>
												</Box>
											</MenuItem>
										))}
									</Select>
								</FormControl>

								{selectedJob && (
									<JobStatusChips
										jobName={selectedJob.full_name}
										workerCount={selectedJob.workers.length}
										pendingTaskCount={0}
									/>
								)}

								{dynamicSchema && (
									<>
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

										{isLoadingSchema ? (
											<FormSkeleton />
										) : (
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
										)}

										{/* Task History Panel */}
										<TaskHistoryPanel
											roomId={roomId}
											jobName={selectedJobName}
										/>
									</>
								)}
							</AccordionDetails>
						</Accordion>

						{/* Selection Groups Section */}
						<Accordion defaultExpanded sx={{ mt: 2 }}>
							<AccordionSummary expandIcon={<ExpandMoreIcon />}>
								<Typography variant="subtitle1" fontWeight="medium">
									Selection Groups
								</Typography>
							</AccordionSummary>
							<AccordionDetails sx={{ p: 0 }}>
								<Box sx={{ height: 400, width: "100%" }}>
									<DataGrid
										rows={rows}
										columns={columns}
										disableRowSelectionOnClick
										disableColumnMenu
										hideFooter
										onRowClick={handleRowClick}
										getRowClassName={(params) => {
											if (params.row.isCurrent) return "current-selection-row";
											if (params.row.isActive) return "active-group-row";
											return "";
										}}
										sx={{
											"& .MuiDataGrid-columnHeaders": {
												minHeight: "48px !important",
												maxHeight: "48px !important",
											},
											"& .MuiDataGrid-columnHeader": {
												padding: "8px",
											},
											"& .current-selection-row": {
												bgcolor: "rgba(33, 150, 243, 0.08)",
												borderTop: "2px solid",
												borderTopColor: "primary.main",
												"&:hover": {
													bgcolor: "rgba(33, 150, 243, 0.12)",
												},
											},
											"& .active-group-row": {
												bgcolor: "rgba(76, 175, 80, 0.08)",
												"&:hover": {
													bgcolor: "rgba(76, 175, 80, 0.12)",
												},
											},
											"& .MuiDataGrid-row": {
												cursor: "pointer",
											},
											"& .MuiDataGrid-row.current-selection-row": {
												cursor: "default",
											},
										}}
									/>
								</Box>
							</AccordionDetails>
						</Accordion>

						{Object.keys(selectionGroups).length === 0 && (
							<Typography
								color="text.secondary"
								sx={{ p: 2, textAlign: "center" }}
							>
								Select geometry instances and type a name above to save a group.
							</Typography>
						)}

						{/* Delete Confirmation Dialog */}
						<Dialog
							open={deleteDialogOpen}
							onClose={() => setDeleteDialogOpen(false)}
						>
							<DialogTitle>Delete Selection Group</DialogTitle>
							<DialogContent>
								<DialogContentText>
									Are you sure you want to delete the group "{groupToDelete}"?
									This action cannot be undone.
								</DialogContentText>
							</DialogContent>
							<DialogActions>
								<Button onClick={() => setDeleteDialogOpen(false)}>
									Cancel
								</Button>
								<Button
									onClick={handleDeleteConfirm}
									color="error"
									variant="contained"
								>
									Delete
								</Button>
							</DialogActions>
						</Dialog>
					</>
				)}
			</Box>
		</Box>
	);
}
