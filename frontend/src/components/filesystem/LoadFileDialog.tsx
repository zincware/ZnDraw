import {
	Autocomplete,
	Box,
	Button,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
	FormControl,
	FormControlLabel,
	Radio,
	RadioGroup,
	TextField,
	Typography,
} from "@mui/material";
import React, { useState } from "react";
import type { FilesystemFileItem, Room } from "../../myapi/client";
import { useAppStore } from "../../store";

export type RoomTarget =
	| { type: "current" }
	| { type: "new"; room_id: string; description: string }
	| { type: "existing"; room_id: string };

export interface LoadFileParams {
	path: string;
	start?: number;
	stop?: number;
	step?: number;
	room_target: RoomTarget;
}

interface LoadFileDialogProps {
	open: boolean;
	file: FilesystemFileItem | null;
	isLoading: boolean;
	rooms: Room[];
	onClose: () => void;
	onLoad: (params: LoadFileParams) => void;
}

/**
 * Dialog for loading a file from filesystem with slice parameters.
 */
export function LoadFileDialog({
	open,
	file,
	isLoading,
	rooms,
	onClose,
	onLoad,
}: LoadFileDialogProps) {
	const [sliceParams, setSliceParams] = useState<{
		start: string;
		stop: string;
		step: string;
	}>({ start: "", stop: "", step: "" });
	const [targetType, setTargetType] = useState<"current" | "new" | "existing">(
		"current",
	);
	const [newRoomId, setNewRoomId] = useState("");
	const [newRoomDescription, setNewRoomDescription] = useState("");
	const [selectedRoom, setSelectedRoom] = useState<Room | null>(null);

	// Reset state when dialog opens
	React.useEffect(() => {
		if (open) {
			setSliceParams({ start: "", stop: "", step: "" });
			setTargetType("current");
			setNewRoomId(crypto.randomUUID());
			setNewRoomDescription("");
			setSelectedRoom(null);
		}
	}, [open]);

	const validateSliceParams = (): string | null => {
		const start = sliceParams.start ? Number.parseInt(sliceParams.start) : null;
		const stop = sliceParams.stop ? Number.parseInt(sliceParams.stop) : null;
		const step = sliceParams.step ? Number.parseInt(sliceParams.step) : null;

		if (sliceParams.start && isNaN(start!)) {
			return "Start must be a valid integer";
		}
		if (sliceParams.stop && isNaN(stop!)) {
			return "Stop must be a valid integer";
		}
		if (sliceParams.step && isNaN(step!)) {
			return "Step must be a valid integer";
		}
		if (step !== null && step <= 0) {
			return "Step must be a positive integer (e.g., 1, 2, 5)";
		}
		if (start !== null && start < 0) {
			return "Start cannot be negative";
		}
		if (stop !== null && stop < 0) {
			return "Stop cannot be negative";
		}
		if (start !== null && stop !== null && start >= stop) {
			return "Start must be less than stop";
		}
		return null;
	};

	const showSnackbar = useAppStore((state) => state.showSnackbar);

	const handleLoad = () => {
		if (!file) return;

		const validationError = validateSliceParams();
		if (validationError) {
			showSnackbar(validationError, "warning");
			return;
		}

		if (targetType === "new" && !newRoomId.trim()) {
			showSnackbar("Please enter a room ID for the new room", "warning");
			return;
		}
		if (targetType === "existing" && !selectedRoom) {
			showSnackbar("Please select a room", "warning");
			return;
		}

		const room_target: RoomTarget =
			targetType === "new"
				? {
						type: "new",
						room_id: newRoomId.trim(),
						description: newRoomDescription.trim(),
					}
				: targetType === "existing"
					? { type: "existing", room_id: selectedRoom!.id }
					: { type: "current" };

		const params: LoadFileParams = {
			path: file.path,
			room_target,
			...(sliceParams.start && { start: Number.parseInt(sliceParams.start) }),
			...(sliceParams.stop && { stop: Number.parseInt(sliceParams.stop) }),
			...(sliceParams.step && { step: Number.parseInt(sliceParams.step) }),
		};

		onLoad(params);
	};

	if (!file) return null;

	return (
		<Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
			<DialogTitle>Load File from Remote Filesystem</DialogTitle>
			<DialogContent>
				<Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
					File: {file.name}
				</Typography>
				<Typography variant="body2" color="text.secondary" sx={{ mb: 2 }}>
					Path: {file.path}
				</Typography>

				<Typography variant="body2" sx={{ mb: 1, fontWeight: 500 }}>
					Frame Selection (optional)
				</Typography>
				<Typography
					variant="caption"
					color="text.secondary"
					sx={{ mb: 2, display: "block" }}
				>
					Load specific frames from trajectory files. Leave empty to load all
					frames.
				</Typography>

				<Box sx={{ display: "flex", gap: 2 }}>
					<TextField
						size="small"
						label="Start"
						type="number"
						value={sliceParams.start}
						onChange={(e) =>
							setSliceParams({ ...sliceParams, start: e.target.value })
						}
						helperText="First frame"
						sx={{ flex: 1 }}
						slotProps={{ htmlInput: { min: 0 } }}
					/>
					<TextField
						size="small"
						label="Stop"
						type="number"
						value={sliceParams.stop}
						onChange={(e) =>
							setSliceParams({ ...sliceParams, stop: e.target.value })
						}
						helperText="Last frame"
						sx={{ flex: 1 }}
						slotProps={{ htmlInput: { min: 0 } }}
					/>
					<TextField
						size="small"
						label="Step"
						type="number"
						value={sliceParams.step}
						onChange={(e) =>
							setSliceParams({ ...sliceParams, step: e.target.value })
						}
						helperText="Interval"
						sx={{ flex: 1 }}
						slotProps={{ htmlInput: { min: 1 } }}
					/>
				</Box>

				{(sliceParams.start || sliceParams.stop || sliceParams.step) && (
					<Typography
						variant="caption"
						color="primary.main"
						sx={{ mt: 1, display: "block" }}
					>
						Will load: [{sliceParams.start || "0"}:{sliceParams.stop || "end"}:
						{sliceParams.step || "1"}]
						{sliceParams.step &&
							Number.parseInt(sliceParams.step) > 1 &&
							` (every ${sliceParams.step} frame${
								Number.parseInt(sliceParams.step) > 1 ? "s" : ""
							})`}
					</Typography>
				)}

				<Typography variant="body2" sx={{ mt: 3, mb: 1, fontWeight: 500 }}>
					Target Room
				</Typography>
				<FormControl>
					<RadioGroup
						value={targetType}
						onChange={(e) =>
							setTargetType(
								e.target.value as "current" | "new" | "existing",
							)
						}
					>
						<FormControlLabel
							value="current"
							control={<Radio size="small" />}
							label="Current room"
						/>
						<FormControlLabel
							value="new"
							control={<Radio size="small" />}
							label="New room"
						/>
						<FormControlLabel
							value="existing"
							control={<Radio size="small" />}
							label="Existing room"
						/>
					</RadioGroup>
				</FormControl>
				{targetType === "new" && (
					<Box sx={{ display: "flex", flexDirection: "column", gap: 1.5, mt: 1 }}>
						<TextField
							size="small"
							fullWidth
							label="Room ID"
							value={newRoomId}
							onChange={(e) => setNewRoomId(e.target.value)}
							helperText="Auto-generated UUID — replace with a custom ID if desired"
						/>
						<TextField
							size="small"
							fullWidth
							label="Description (optional)"
							value={newRoomDescription}
							onChange={(e) => setNewRoomDescription(e.target.value)}
						/>
					</Box>
				)}
				{targetType === "existing" && (
					<Autocomplete
						size="small"
						options={rooms}
						value={selectedRoom}
						onChange={(_, value) => setSelectedRoom(value)}
						getOptionLabel={(room) =>
							room.description
								? `${room.description} (${room.id})`
								: room.id
						}
						renderInput={(params) => (
							<TextField {...params} label="Select room" />
						)}
						isOptionEqualToValue={(option, value) =>
							option.id === value.id
						}
						sx={{ mt: 1 }}
					/>
				)}
			</DialogContent>
			<DialogActions>
				<Button onClick={onClose}>Cancel</Button>
				<Button onClick={handleLoad} variant="contained" disabled={isLoading}>
					{isLoading ? "Loading..." : "Load"}
				</Button>
			</DialogActions>
		</Dialog>
	);
}
