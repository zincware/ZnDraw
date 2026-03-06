import {
	Alert,
	Box,
	Button,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
	TextField,
	Typography,
} from "@mui/material";
import React, { useState, useEffect } from "react";
import { useAppStore } from "../store";
import { changePassword } from "../utils/auth";

interface UserProfileDialogProps {
	open: boolean;
	onClose: () => void;
}

export default function UserProfileDialog({
	open,
	onClose,
}: UserProfileDialogProps) {
	const [error, setError] = useState<string | null>(null);
	const [success, setSuccess] = useState<string | null>(null);
	const [loading, setLoading] = useState(false);

	// Change password form
	const [newPassword, setNewPassword] = useState("");
	const [newPasswordConfirm, setNewPasswordConfirm] = useState("");

	// Use individual selectors to prevent unnecessary re-renders
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const userName = useAppStore((state) => state.user?.email ?? null);

	useEffect(() => {
		if (open) {
			clearForms();
		}
	}, [open]);

	const clearForms = () => {
		setNewPassword("");
		setNewPasswordConfirm("");
		setError(null);
		setSuccess(null);
	};

	const handleChangePassword = async () => {
		setError(null);
		setSuccess(null);

		if (!newPassword) {
			setError("New password is required");
			return;
		}

		if (newPassword !== newPasswordConfirm) {
			setError("Passwords do not match");
			return;
		}

		setLoading(true);
		try {
			await changePassword(newPassword);
			setSuccess("Password changed successfully!");
			showSnackbar("Password changed", "success");
			clearForms();
		} catch (err) {
			setError(err instanceof Error ? err.message : "Password change failed");
		} finally {
			setLoading(false);
		}
	};

	const handleClose = () => {
		if (!loading) {
			clearForms();
			onClose();
		}
	};

	return (
		<Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
			<DialogTitle>Manage Account</DialogTitle>
			<DialogContent>
				<Box sx={{ pt: 1, display: "flex", flexDirection: "column", gap: 2 }}>
					<Typography variant="body2" color="text.secondary">
						Account: <strong>{userName}</strong>
					</Typography>

					{error && (
						<Alert severity="error" onClose={() => setError(null)}>
							{error}
						</Alert>
					)}

					{success && (
						<Alert severity="success" onClose={() => setSuccess(null)}>
							{success}
						</Alert>
					)}

					<TextField
						autoFocus
						label="New Password"
						type="password"
						fullWidth
						value={newPassword}
						onChange={(e) => setNewPassword(e.target.value)}
						disabled={loading}
						autoComplete="new-password"
					/>
					<TextField
						label="Confirm New Password"
						type="password"
						fullWidth
						value={newPasswordConfirm}
						onChange={(e) => setNewPasswordConfirm(e.target.value)}
						onKeyPress={(e) => {
							if (e.key === "Enter" && !loading) {
								handleChangePassword();
							}
						}}
						disabled={loading}
						autoComplete="new-password"
					/>
				</Box>
			</DialogContent>
			<DialogActions sx={{ px: 3, pb: 2 }}>
				<Button onClick={handleClose} disabled={loading}>
					Cancel
				</Button>
				<Button
					onClick={handleChangePassword}
					variant="contained"
					disabled={loading || !newPassword || !newPasswordConfirm}
				>
					{loading ? "Changing..." : "Change Password"}
				</Button>
			</DialogActions>
		</Dialog>
	);
}
