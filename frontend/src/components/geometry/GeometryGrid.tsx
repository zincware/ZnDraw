import DeleteIcon from "@mui/icons-material/Delete";
import EditIcon from "@mui/icons-material/Edit";
import LockIcon from "@mui/icons-material/Lock";
import RadioButtonCheckedIcon from "@mui/icons-material/RadioButtonChecked";
import RadioButtonUncheckedIcon from "@mui/icons-material/RadioButtonUnchecked";
import SearchIcon from "@mui/icons-material/Search";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import {
	Box,
	IconButton,
	InputAdornment,
	Switch,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import {
	DataGrid,
	GridActionsCellItem,
	type GridColDef,
	type GridRenderCellParams,
	type GridRowParams,
} from "@mui/x-data-grid";
import { useMemo, useState } from "react";
import { useDefaultCamera } from "../../hooks/useDefaultCamera";
import {
	useCreateGeometry,
	useDeleteGeometry,
} from "../../hooks/useGeometries";
import { useAppStore } from "../../store";
import { useGeometryStore } from "../../stores/geometryStore";
import DeleteConfirmDialog from "./DeleteConfirmDialog";

interface GeometryGridProps {
	geometries: Array<{ key: string; type: string }>;
}

const GeometryGrid = ({ geometries }: GeometryGridProps) => {
	// Use individual selectors to prevent unnecessary re-renders
	const roomId = useAppStore((state) => state.roomId);
	const userId = useAppStore((state) => state.user?.id ?? null);
	const geometriesData = useAppStore((state) => state.geometries);
	const superuserLock = useAppStore((state) => state.superuserLock);
	const userLock = useAppStore((state) => state.userLock);
	const isSuperuser = useAppStore((state) => state.user?.is_superuser ?? false);
	const currentUserEmail = useAppStore((state) => state.user?.email ?? null);
	const activeCurveForDrawing = useAppStore(
		(state) => state.activeCurveForDrawing,
	);
	const setActiveCurveForDrawing = useAppStore(
		(state) => state.setActiveCurveForDrawing,
	);
	const attachedCameraKey = useAppStore((state) => state.attachedCameraKey);
	const attachToCamera = useAppStore((state) => state.attachToCamera);
	const { setMode, setSelectedKey, searchFilter, setSearchFilter } =
		useGeometryStore();
	const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
	const [geometryToDelete, setGeometryToDelete] = useState<string | null>(null);

	const { mutate: deleteGeometry } = useDeleteGeometry();
	const { mutate: updateGeometry } = useCreateGeometry();
	const { defaultCamera, setDefaultCamera, isSettingDefault } =
		useDefaultCamera();

	/** Check if the current user can edit a geometry. */
	const canEdit = (key: string): boolean => {
		if (isSuperuser) return true;
		const geom = geometriesData[key];
		if (!geom) return false;
		const owner = geom.data?.owner as string | null | undefined;
		// Owned by another user — always blocked
		if (owner && owner !== userId) return false;
		// Edit lock by someone else — always blocked
		if (userLock && userLock !== currentUserEmail) return false;
		// Superuser lock — owners can still edit their own
		if (superuserLock && owner !== userId) return false;
		return true;
	};

	const handleEdit = (key: string) => {
		setSelectedKey(key);
		setMode("edit");
	};

	const handleDeleteClick = (key: string) => {
		setGeometryToDelete(key);
		setDeleteDialogOpen(true);
	};

	const handleDeleteConfirm = () => {
		if (geometryToDelete && roomId) {
			deleteGeometry({ roomId, key: geometryToDelete });
		}
		setDeleteDialogOpen(false);
		setGeometryToDelete(null);
	};

	const handleDeleteCancel = () => {
		setDeleteDialogOpen(false);
		setGeometryToDelete(null);
	};

	const handleActiveToggle = (key: string, currentActive: boolean) => {
		const geometry = geometriesData[key];
		if (!geometry || !roomId) return;

		// Send full data to preserve all fields (including owner)
		updateGeometry({
			roomId,
			key,
			geometryType: geometry.type,
			geometryData: {
				...geometry.data,
				active: !currentActive,
			},
		});
	};

	const handleSelectionToggle = (key: string, geometryType: string) => {
		if (geometryType === "Curve") {
			// Toggle curve drawing target
			if (activeCurveForDrawing === key) {
				setActiveCurveForDrawing(null);
			} else {
				setActiveCurveForDrawing(key);
			}
		} else if (geometryType === "Camera") {
			if (attachedCameraKey !== key) {
				attachToCamera(key);
			}
		}
	};

	const editabilityTooltip = (key: string): string | null => {
		const geom = geometriesData[key];
		const owner = geom?.data?.owner as string | null | undefined;
		if (owner && owner !== userId) return "Owned by another user";
		if (userLock && userLock !== currentUserEmail)
			return "Room is being edited";
		if (superuserLock && owner !== userId) return "Room is locked";
		return null;
	};

	const columns: GridColDef[] = [
		{
			field: "key",
			headerName: "Key",
			flex: 1,
			minWidth: 150,
		},
		{
			field: "type",
			headerName: "Type",
			flex: 1,
			minWidth: 120,
		},
		{
			field: "selection",
			headerName: "Active",
			width: 70,
			renderCell: (params: GridRenderCellParams) => {
				const geometryType = params.row.type;
				const isCurve = geometryType === "Curve";
				const isCamera = geometryType === "Camera";

				// Only show for Curves and Cameras
				if (!isCurve && !isCamera) {
					return null;
				}

				const isSelected = isCurve
					? activeCurveForDrawing === params.row.key
					: attachedCameraKey === params.row.key;

				const tooltipText = isCurve
					? isSelected
						? "Deselect this curve for drawing"
						: "Select this curve for drawing"
					: isSelected
						? "Currently attached to this camera"
						: "Attach to this camera";

				return (
					<Tooltip title={tooltipText}>
						<IconButton
							size="small"
							onClick={(e) => {
								e.stopPropagation();
								handleSelectionToggle(params.row.key, geometryType);
							}}
							color={isSelected ? "primary" : "default"}
						>
							{isSelected ? (
								<RadioButtonCheckedIcon />
							) : (
								<RadioButtonUncheckedIcon />
							)}
						</IconButton>
					</Tooltip>
				);
			},
		},
		{
			field: "default",
			headerName: "",
			width: 50,
			sortable: false,
			renderCell: (params: GridRenderCellParams) => {
				const isCamera = params.row.type === "Camera";
				const isSessionCam = (params.row.key as string).startsWith("cam:");
				if (!isCamera || isSessionCam) return null;

				const isDefault = defaultCamera === params.row.key;
				return (
					<Tooltip
						title={isDefault ? "Unset default camera" : "Set as default camera"}
					>
						<span>
							<IconButton
								size="small"
								disabled={isSettingDefault}
								onClick={(e) => {
									e.stopPropagation();
									setDefaultCamera(isDefault ? null : params.row.key);
								}}
							>
								{isDefault ? <StarIcon color="warning" /> : <StarBorderIcon />}
							</IconButton>
						</span>
					</Tooltip>
				);
			},
		},
		{
			field: "active",
			headerName: "Visible",
			width: 80,
			renderCell: (params: GridRenderCellParams) => {
				const isActive = params.row.active !== false;
				const editable = canEdit(params.row.key);
				return (
					<Switch
						checked={isActive}
						disabled={!editable}
						onChange={(e) => {
							e.stopPropagation();
							handleActiveToggle(params.row.key, isActive);
						}}
						onClick={(e) => {
							e.stopPropagation();
						}}
						size="small"
					/>
				);
			},
		},
		{
			field: "actions",
			type: "actions",
			headerName: "Actions",
			width: 100,
			getActions: (params: GridRowParams) => {
				const editable = canEdit(params.row.key);
				const reason = editabilityTooltip(params.row.key);

				if (!editable) {
					return [
						<GridActionsCellItem
							key="locked"
							icon={
								<Tooltip title={reason ?? ""}>
									<LockIcon color="disabled" />
								</Tooltip>
							}
							label="Locked"
							disabled
							showInMenu={false}
						/>,
					];
				}

				return [
					<GridActionsCellItem
						key="edit"
						icon={<EditIcon />}
						label="Edit"
						onClick={() => handleEdit(params.row.key)}
						showInMenu={false}
					/>,
					<GridActionsCellItem
						key="delete"
						icon={<DeleteIcon />}
						label="Delete"
						onClick={() => handleDeleteClick(params.row.key)}
						showInMenu={false}
					/>,
				];
			},
		},
	];

	const filteredGeometries = useMemo(() => {
		if (!searchFilter) return geometries;
		const lowerFilter = searchFilter.toLowerCase();
		return geometries.filter(
			(g) =>
				g.key.toLowerCase().includes(lowerFilter) ||
				g.type.toLowerCase().includes(lowerFilter),
		);
	}, [geometries, searchFilter]);

	const rows = filteredGeometries.map((g, index) => ({
		id: index,
		key: g.key,
		type: g.type,
		active: geometriesData[g.key]?.data?.active !== false, // Default to true if undefined
	}));

	return (
		<Box sx={{ height: "100%", display: "flex", flexDirection: "column" }}>
			<Box sx={{ p: 2, pb: 1 }}>
				<TextField
					fullWidth
					size="small"
					placeholder="Search geometries..."
					value={searchFilter}
					onChange={(e) => setSearchFilter(e.target.value)}
					slotProps={{
						input: {
							startAdornment: (
								<InputAdornment position="start">
									<SearchIcon />
								</InputAdornment>
							),
						},
					}}
				/>
			</Box>
			<Box sx={{ flexGrow: 1, px: 2, pb: 2 }}>
				<DataGrid
					rows={rows}
					columns={columns}
					density="compact"
					disableRowSelectionOnClick
					slots={{
						noRowsOverlay: () => (
							<Box
								sx={{
									display: "flex",
									alignItems: "center",
									justifyContent: "center",
									height: "100%",
								}}
							>
								<Typography color="text.secondary">
									No geometries yet. Click "Add Geometry" to create one.
								</Typography>
							</Box>
						),
					}}
					sx={{
						border: "1px solid rgba(0, 0, 0, 0.12)",
						"& .MuiDataGrid-cell:focus": {
							outline: "none",
						},
						"& .MuiDataGrid-row:hover": {
							cursor: "pointer",
						},
						"& .MuiDataGrid-row.non-editable": {
							opacity: 0.6,
							cursor: "default",
						},
					}}
					getRowClassName={(params) =>
						!canEdit(params.row.key) ? "non-editable" : ""
					}
					onRowClick={(params) => {
						if (canEdit(params.row.key)) {
							handleEdit(params.row.key);
						}
					}}
				/>
			</Box>
			<DeleteConfirmDialog
				open={deleteDialogOpen}
				geometryKey={geometryToDelete}
				onConfirm={handleDeleteConfirm}
				onCancel={handleDeleteCancel}
			/>
		</Box>
	);
};

export default GeometryGrid;
