import Button from "@mui/material/Button";
import Dialog from "@mui/material/Dialog";
import DialogActions from "@mui/material/DialogActions";
import DialogContent from "@mui/material/DialogContent";
import DialogTitle from "@mui/material/DialogTitle";
import TextField from "@mui/material/TextField";
import Typography from "@mui/material/Typography";
import { useState } from "react";
import { useNavigate } from "react-router-dom";
import { type Visibility, createRoom } from "../myapi/client";
import { validateRoomId } from "../utils/roomValidation";
import VisibilitySelector from "./VisibilitySelector";

interface DuplicateRoomDialogProps {
	open: boolean;
	sourceRoomId: string;
	sourceDescription: string;
	existingRoomIds: string[];
	onClose: () => void;
}

export default function DuplicateRoomDialog({
	open,
	sourceRoomId,
	sourceDescription,
	existingRoomIds,
	onClose,
}: DuplicateRoomDialogProps) {
	const navigate = useNavigate();
	const [newRoomId, setNewRoomId] = useState("");
	const [description, setDescription] = useState("");
	const [error, setError] = useState<string | null>(null);
	const [visibility, setVisibility] = useState<Visibility>("public");

	// Reset form when dialog opens
	const handleEnter = () => {
		setNewRoomId("");
		setDescription(`Copy of ${sourceDescription || sourceRoomId}`);
		setError(null);
		setVisibility("public");
	};

	const handleDuplicate = async () => {
		const validationError = validateRoomId(newRoomId, existingRoomIds);
		if (validationError) {
			setError(validationError);
			return;
		}

		try {
			const roomId = newRoomId || crypto.randomUUID();
			const result = await createRoom({
				room_id: roomId,
				copy_from: sourceRoomId,
				description,
				visibility,
			});
			onClose();
			navigate(`/rooms/${result.room_id}`);
		} catch (err) {
			setError(err instanceof Error ? err.message : "Failed to duplicate room");
		}
	};

	return (
		<Dialog
			open={open}
			onClose={onClose}
			maxWidth="sm"
			fullWidth
			TransitionProps={{ onEnter: handleEnter }}
		>
			<DialogTitle>Duplicate Room</DialogTitle>
			<DialogContent>
				<Typography
					variant="body2"
					color="text.secondary"
					sx={{ mb: 3, mt: 1 }}
				>
					Duplicating: {sourceDescription || sourceRoomId}
				</Typography>

				<TextField
					margin="dense"
					label="New Room ID (optional)"
					type="text"
					fullWidth
					variant="outlined"
					value={newRoomId}
					onChange={(e) => {
						const id = e.target.value;
						setNewRoomId(id);
						setError(validateRoomId(id, existingRoomIds));
					}}
					helperText={error || "Leave empty to auto-generate a unique ID"}
					error={!!error}
					sx={{ mb: 2 }}
				/>

				<TextField
					autoFocus
					margin="dense"
					label="Description for new room"
					type="text"
					fullWidth
					variant="outlined"
					value={description}
					onChange={(e) => setDescription(e.target.value)}
					sx={{ mb: 2 }}
				/>

				<VisibilitySelector
					value={visibility}
					onChange={setVisibility}
				/>
			</DialogContent>
			<DialogActions>
				<Button onClick={onClose}>Cancel</Button>
				<Button
					onClick={handleDuplicate}
					variant="contained"
					disabled={!!error}
				>
					Duplicate
				</Button>
			</DialogActions>
		</Dialog>
	);
}
