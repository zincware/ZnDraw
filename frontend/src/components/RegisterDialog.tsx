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
import type React from "react";
import { useState } from "react";
import { connectWithAuth } from "../socket";
import { useAppStore } from "../store";
import { registerUser } from "../utils/auth";

interface RegisterDialogProps {
	open: boolean;
	onClose: () => void;
}

export default function RegisterDialog({ open, onClose }: RegisterDialogProps) {
	const [email, setEmail] = useState("");
	const [password, setPassword] = useState("");
	const [passwordConfirm, setPasswordConfirm] = useState("");
	const [error, setError] = useState<string | null>(null);
	const [loading, setLoading] = useState(false);

	// Use individual selectors to prevent unnecessary re-renders
	const setUser = useAppStore((state) => state.setUser);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const userEmail = useAppStore((state) => state.user?.email ?? null);

	const handleRegister = async () => {
		setError(null);

		if (!email || !email.trim()) {
			setError("Email is required");
			return;
		}

		if (!password) {
			setError("Password is required");
			return;
		}

		if (password !== passwordConfirm) {
			setError("Passwords do not match");
			return;
		}

		setLoading(true);
		try {
			await registerUser(email, password);
			const { user } = await connectWithAuth();

			setUser(user);

			showSnackbar(`Registered as ${user.email}`, "success");
			onClose();

			// Clear form
			setEmail("");
			setPassword("");
			setPasswordConfirm("");
		} catch (err) {
			setError(err instanceof Error ? err.message : "Registration failed");
		} finally {
			setLoading(false);
		}
	};

	const handleKeyDown = (event: React.KeyboardEvent) => {
		if (event.key === "Enter" && !loading) {
			if (email && password && passwordConfirm) {
				handleRegister();
			}
		}
	};

	const handleClose = () => {
		if (!loading) {
			setError(null);
			setEmail("");
			setPassword("");
			setPasswordConfirm("");
			onClose();
		}
	};

	return (
		<Dialog open={open} onClose={handleClose} maxWidth="xs" fullWidth>
			<DialogTitle>Register Account</DialogTitle>
			<DialogContent>
				<Box sx={{ pt: 1, display: "flex", flexDirection: "column", gap: 2 }}>
					<Typography variant="body2" color="text.secondary">
						Current temporary name: <strong>{userEmail}</strong>
					</Typography>
					<Typography variant="body2" color="text.secondary">
						Choose an email and password to register your account.
					</Typography>

					{error && (
						<Alert severity="error" onClose={() => setError(null)}>
							{error}
						</Alert>
					)}

					<TextField
						label="Email"
						type="email"
						value={email}
						onChange={(e) => setEmail(e.target.value)}
						onKeyDown={handleKeyDown}
						disabled={loading}
						fullWidth
						autoComplete="email"
					/>

					<TextField
						label="Password"
						type="password"
						value={password}
						onChange={(e) => setPassword(e.target.value)}
						onKeyDown={handleKeyDown}
						disabled={loading}
						fullWidth
						autoComplete="new-password"
					/>

					<TextField
						label="Confirm Password"
						type="password"
						value={passwordConfirm}
						onChange={(e) => setPasswordConfirm(e.target.value)}
						onKeyDown={handleKeyDown}
						disabled={loading}
						fullWidth
						autoComplete="new-password"
					/>
				</Box>
			</DialogContent>
			<DialogActions sx={{ px: 3, pb: 2 }}>
				<Button onClick={handleClose} disabled={loading}>
					Cancel
				</Button>
				<Button
					onClick={handleRegister}
					disabled={loading || !email || !password || !passwordConfirm}
					variant="contained"
				>
					{loading ? "Registering..." : "Register"}
				</Button>
			</DialogActions>
		</Dialog>
	);
}
