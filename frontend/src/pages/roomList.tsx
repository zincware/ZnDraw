import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import LockIcon from "@mui/icons-material/Lock";
import LockOpenIcon from "@mui/icons-material/LockOpen";
import StarIcon from "@mui/icons-material/Star";
import StarBorderIcon from "@mui/icons-material/StarBorder";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import Alert from "@mui/material/Alert";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import CircularProgress from "@mui/material/CircularProgress";
import Container from "@mui/material/Container";
import IconButton from "@mui/material/IconButton";
import TextField from "@mui/material/TextField";
import Tooltip from "@mui/material/Tooltip";
import Typography from "@mui/material/Typography";
import {
	DataGrid,
	type GridColDef,
	type GridRenderCellParams,
} from "@mui/x-data-grid";
import React, { useCallback, useState, useEffect } from "react";
import { useNavigate } from "react-router-dom";
import DropOverlay from "../components/DropOverlay";
import DuplicateRoomDialog from "../components/DuplicateRoomDialog";
import { TRAJECTORY_ACCEPT } from "../constants/fileTypes";
import { useDragAndDrop } from "../hooks/useDragAndDrop";
import { useSocketManager } from "../hooks/useSocketManager";
import {
	createRoom,
	setDefaultRoom,
	updateRoom,
	uploadTrajectory,
} from "../myapi/client";
import { useRoomsStore } from "../roomsStore";
import { useAppStore } from "../store";

export default function RoomListPage() {
	// Use individual selectors to prevent unnecessary re-renders
	const showSnackbar = useAppStore((state) => state.showSnackbar);

	// Connect to socket and join overview:public room
	useSocketManager({ isOverview: true });

	// Subscribe to rooms store (triggers re-render on changes)
	const allRooms = useRoomsStore((state) => state.roomsArray);
	const loading = useRoomsStore((state) => state.loading);
	const roomsError = useRoomsStore((state) => state.error);

	const [searchQuery, setSearchQuery] = useState<string>("");
	const [duplicateSource, setDuplicateSource] = useState<{
		roomId: string;
		description: string;
	} | null>(null);
	const navigate = useNavigate();

	// Drag and drop: create a new room, upload file(s), then navigate
	const handleFiles = useCallback(
		async (files: File[]) => {
			const newRoomId = crypto.randomUUID();
			try {
				await createRoom({ room_id: newRoomId });
				for (const file of files) {
					await uploadTrajectory(newRoomId, file);
				}
				showSnackbar(`Room created with ${files.length} file(s)`, "success");
				navigate(`/rooms/${newRoomId}`);
			} catch (error: any) {
				const detail =
					error?.response?.data?.detail || "Upload failed";
				showSnackbar(detail, "error");
			}
		},
		[navigate, showSnackbar],
	);

	const {
		isDragging,
		handleDragOver,
		handleDragEnter,
		handleDragLeave,
		handleDrop,
	} = useDragAndDrop(handleFiles);

	// File upload ref for button click
	const fileInputRef = React.useRef<HTMLInputElement>(null);

	const handleFileUploadClick = () => {
		fileInputRef.current?.click();
	};

	const handleFileInputChange = async (
		event: React.ChangeEvent<HTMLInputElement>,
	) => {
		const files = event.target.files;
		if (!files || files.length === 0) return;

		await handleFiles(Array.from(files));

		// Reset input
		if (fileInputRef.current) {
			fileInputRef.current.value = "";
		}
	};

	// Fetch rooms on mount (socket updates will keep them synced)
	useEffect(() => {
		useRoomsStore.getState().fetchRooms();
	}, []);

	// Filter rooms based on search query
	const rooms = searchQuery
		? allRooms.filter(
				(room) =>
					room.id.toLowerCase().includes(searchQuery.toLowerCase()) ||
					(room.description &&
						room.description.toLowerCase().includes(searchQuery.toLowerCase())),
			)
		: allRooms;

	const error = roomsError;
	const isError = !!roomsError;

	const handleUpdateDescription = async (
		roomId: string,
		description: string,
	) => {
		try {
			await updateRoom(roomId, { description: description || null });
			// Socket event will update Zustand store automatically
			showSnackbar("Description updated", "success");
		} catch (err) {
			showSnackbar("Failed to update description", "error");
		}
	};

	const handleToggleLock = async (roomId: string, currentLocked: boolean) => {
		try {
			await updateRoom(roomId, { locked: !currentLocked });
			// Socket event will update Zustand store automatically
			showSnackbar(currentLocked ? "Room unlocked" : "Room locked", "success");
		} catch (err) {
			showSnackbar("Failed to update lock status", "error");
		}
	};

	const handleToggleDefault = async (
		roomId: string,
		isCurrentlyDefault: boolean,
	) => {
		try {
			await setDefaultRoom(isCurrentlyDefault ? null : roomId);
			// Socket event will update Zustand store automatically
			showSnackbar(
				isCurrentlyDefault ? "Template cleared" : "Template set",
				"success",
			);
		} catch (err) {
			showSnackbar("Failed to update template", "error");
		}
	};

	const handleOpenDuplicateDialog = (roomId: string, description: string) => {
		setDuplicateSource({ roomId, description });
	};

	const handleOpenRoom = (roomId: string) => {
		navigate(`/rooms/${roomId}`);
	};

	const columns: GridColDef[] = [
		{
			field: "is_default",
			headerName: "",
			width: 50,
			sortable: true,
			renderCell: (params: GridRenderCellParams) => (
				<Tooltip title={params.value ? "Template room" : "Set as template"}>
					<IconButton
						size="small"
						onClick={() => handleToggleDefault(params.row.id, params.value)}
					>
						{params.value ? <StarIcon color="primary" /> : <StarBorderIcon />}
					</IconButton>
				</Tooltip>
			),
		},
		{
			field: "id",
			headerName: "Room ID",
			width: 250,
			sortable: true,
			renderCell: (params: GridRenderCellParams) => (
				<Box
					sx={{
						cursor: "pointer",
						color: "primary.main",
						textDecoration: "underline",
						"&:hover": {
							color: "primary.dark",
						},
					}}
					onClick={() => handleOpenRoom(params.value)}
				>
					{params.value}
				</Box>
			),
		},
		{
			field: "description",
			headerName: "Description",
			width: 300,
			editable: true,
			sortable: true,
		},
		{
			field: "frame_count",
			headerName: "Frames",
			width: 100,
			type: "number",
			sortable: true,
		},
		{
			field: "locked",
			headerName: "Lock",
			width: 80,
			renderCell: (params: GridRenderCellParams) => {
				const isLocked = params.value;

				if (isLocked) {
					return (
						<Tooltip title="Locked (immutable)">
							<IconButton
								size="small"
								onClick={() => handleToggleLock(params.row.id, isLocked)}
							>
								<LockIcon color="error" />
							</IconButton>
						</Tooltip>
					);
				} else {
					return (
						<Tooltip title="Unlocked (editable)">
							<IconButton
								size="small"
								onClick={() => handleToggleLock(params.row.id, isLocked)}
							>
								<LockOpenIcon color="success" />
							</IconButton>
						</Tooltip>
					);
				}
			},
		},
		{
			field: "actions",
			headerName: "Actions",
			width: 80,
			sortable: false,
			renderCell: (params: GridRenderCellParams) => (
				<Tooltip title="Duplicate room">
					<IconButton
						size="small"
						onClick={() =>
							handleOpenDuplicateDialog(
								params.row.id,
								params.row.description || "",
							)
						}
					>
						<ContentCopyIcon />
					</IconButton>
				</Tooltip>
			),
		},
	];

	if (loading) {
		return (
			<Container maxWidth="lg">
				<Box
					sx={{
						display: "flex",
						justifyContent: "center",
						alignItems: "center",
						minHeight: "100vh",
					}}
				>
					<CircularProgress />
				</Box>
			</Container>
		);
	}

	if (isError) {
		return (
			<Container maxWidth="lg">
				<Box sx={{ mt: 4 }}>
					<Alert severity="error">{error || "Unknown error"}</Alert>
				</Box>
			</Container>
		);
	}

	return (
		<Container
			maxWidth="lg"
			onDragOver={handleDragOver}
			onDragEnter={handleDragEnter}
			onDragLeave={handleDragLeave}
			onDrop={handleDrop}
		>
			<DropOverlay isDragging={isDragging} />

			{/* Hidden file input for button upload */}
			<input
				type="file"
				ref={fileInputRef}
				style={{ display: "none" }}
				onChange={handleFileInputChange}
				accept={TRAJECTORY_ACCEPT}
			/>

			<Box sx={{ mt: 4, mb: 4 }}>
				<Typography variant="h3" component="h1" gutterBottom>
					Room Management
				</Typography>
				<Typography variant="subtitle1" color="text.secondary" gutterBottom>
					Manage your rooms: lock, hide, set as template, or duplicate
				</Typography>

				{/* Search Bar */}
				<Box sx={{ mt: 3, mb: 2 }}>
					<TextField
						fullWidth
						size="small"
						label="Search rooms"
						placeholder="Search by metadata (file path, etc.) - supports regex"
						value={searchQuery}
						onChange={(e) => setSearchQuery(e.target.value)}
						variant="outlined"
					/>
				</Box>

				<Box sx={{ height: 600, width: "100%", mt: 3 }}>
					<DataGrid
						rows={rooms}
						columns={columns}
						initialState={{
							pagination: {
								paginationModel: { page: 0, pageSize: 10 },
							},
							sorting: {
								sortModel: [{ field: "is_default", sort: "desc" }],
							},
						}}
						pageSizeOptions={[5, 10, 25]}
						processRowUpdate={async (newRow, oldRow) => {
							if (newRow.description !== oldRow.description) {
								await handleUpdateDescription(newRow.id, newRow.description);
							}
							return newRow;
						}}
						onProcessRowUpdateError={(error) => {
							showSnackbar("Failed to update row", "error");
						}}
					/>
				</Box>

				<Box sx={{ mt: 2, display: "flex", gap: 2 }}>
					<Button
						variant="outlined"
						onClick={() => {
							const roomId = crypto.randomUUID();
							navigate(`/rooms/${roomId}?copy_from=@empty`);
						}}
					>
						Create New Empty Room
					</Button>
					<Button
						variant="contained"
						startIcon={<UploadFileIcon />}
						onClick={handleFileUploadClick}
					>
						Upload File
					</Button>
					</Box>
			</Box>

			{/* Duplicate Dialog */}
			<DuplicateRoomDialog
				open={duplicateSource !== null}
				sourceRoomId={duplicateSource?.roomId || ""}
				sourceDescription={duplicateSource?.description || ""}
				existingRoomIds={rooms.map((r) => r.id)}
				onClose={() => setDuplicateSource(null)}
			/>
		</Container>
	);
}
