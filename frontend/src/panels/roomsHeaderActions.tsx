import AddIcon from "@mui/icons-material/Add";
import UploadFileIcon from "@mui/icons-material/UploadFile";
import { Box, IconButton, Tooltip } from "@mui/material";
import { isAxiosError } from "axios";
import { useRef } from "react";
import { useNavigate } from "react-router-dom";
import { createRoom, uploadTrajectory } from "../myapi/client";
import { useAppStore } from "../store";

export function RoomsHeaderActions() {
	const navigate = useNavigate();
	const showSnackbar = useAppStore((s) => s.showSnackbar);
	const fileInputRef = useRef<HTMLInputElement | null>(null);

	const onNewEmpty = async () => {
		const id = crypto.randomUUID();
		try {
			await createRoom({ room_id: id });
			navigate(`/rooms/${id}`);
		} catch (err) {
			const detail = isAxiosError(err)
				? (err.response?.data?.detail ?? err.message)
				: "Failed to create room";
			showSnackbar(detail, "error");
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
			const detail = isAxiosError(err)
				? (err.response?.data?.detail ?? err.message)
				: "Upload failed";
			showSnackbar(detail, "error");
		} finally {
			if (fileInputRef.current) fileInputRef.current.value = "";
		}
	};

	return (
		<Box sx={{ display: "flex", gap: 0.25 }}>
			<Tooltip title="New empty room">
				<IconButton
					size="small"
					data-testid="rooms-new-empty"
					onClick={onNewEmpty}
				>
					<AddIcon fontSize="small" />
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
