import AddIcon from "@mui/icons-material/Add";
import NoteAddIcon from "@mui/icons-material/NoteAdd";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import { Box, IconButton, Tooltip } from "@mui/material";
import { useRef } from "react";
import { useNavigate } from "react-router-dom";
import { createRoom, uploadTrajectory } from "../myapi/client";
import { useAppStore } from "../store";
import { extractDetail } from "../utils/errors";

export function RoomsHeaderActions() {
	const navigate = useNavigate();
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const fileInputRef = useRef<HTMLInputElement | null>(null);

	// New room: no copy_from → server default template if set, else @empty.
	const onNewRoom = async () => {
		const id = crypto.randomUUID();
		try {
			await createRoom({ room_id: id });
			navigate(`/rooms/${id}`);
		} catch (err) {
			showSnackbar(extractDetail(err, "Failed to create room"), "error");
		}
	};

	// New empty room: copy_from=@none → zero frames, no default geometries.
	const onNewEmpty = async () => {
		const id = crypto.randomUUID();
		try {
			await createRoom({ room_id: id, copy_from: "@none" });
			navigate(`/rooms/${id}`);
		} catch (err) {
			showSnackbar(extractDetail(err, "Failed to create room"), "error");
		}
	};

	const onUploadClick = () => fileInputRef.current?.click();

	const onFileChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
		const files = e.target.files;
		if (!files || files.length === 0) return;
		const id = crypto.randomUUID();
		try {
			await createRoom({ room_id: id });
			for (const f of Array.from(files)) {
				await uploadTrajectory(id, f);
			}
			showSnackbar(`Room created with ${files.length} file(s)`, "success");
			navigate(`/rooms/${id}`);
		} catch (err) {
			showSnackbar(extractDetail(err, "Upload failed"), "error");
		} finally {
			if (fileInputRef.current) fileInputRef.current.value = "";
		}
	};

	return (
		<Box sx={{ display: "flex", gap: 0.25 }}>
			<Tooltip title="New room (from default template)">
				<IconButton size="small" data-testid="rooms-new" onClick={onNewRoom}>
					<AddIcon fontSize="small" />
				</IconButton>
			</Tooltip>
			<Tooltip title="New empty room (no template)">
				<IconButton
					size="small"
					data-testid="rooms-new-empty"
					onClick={onNewEmpty}
				>
					<NoteAddIcon fontSize="small" />
				</IconButton>
			</Tooltip>
			<Tooltip title="Upload file → new room">
				<IconButton
					size="small"
					data-testid="rooms-upload"
					onClick={onUploadClick}
				>
					<UploadFileIcon fontSize="small" />
				</IconButton>
			</Tooltip>
			<input
				type="file"
				ref={fileInputRef}
				style={{ display: "none" }}
				onChange={onFileChange}
				data-testid="rooms-upload-input"
			/>
		</Box>
	);
}
