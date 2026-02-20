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
import { login } from "../utils/auth";

interface LoginDialogProps {
	open: boolean;
	onClose: () => void;
}

export default function LoginDialog({ open, onClose }: LoginDialogProps) {
	const [email, setEmail] = useState("");
	const [password, setPassword] = useState("");
	const [error, setError] = useState<string | null>(null);
	const [loading, setLoading] = useState(false);

	// Use individual selectors to prevent unnecessary re-renders
	const setUser = useAppStore((state) => state.setUser);
	const showSnackbar = useAppStore((state) => state.showSnackbar);

	const handleLoginAsGuest = async () => {
		setError(null);
		setLoading(true);
		try {
			await login();
			const { user } = await connectWithAuth();
			setUser(user);

			showSnackbar(`Logged in as ${user.email}`, "success");
			onClose();
		} catch (err) {
			setError(err instanceof Error ? err.message : "Login failed");
		} finally {
			setLoading(false);
		}
	};

	const handleLogin = async () => {
		if (!email.trim()) {
			setError("Email is required");
			return;
		}
		if (!password) {
			setError("Password is required");
			return;
		}

		setError(null);
		setLoading(true);
		try {
			await login(email, password);
			const { user } = await connectWithAuth();
			setUser(user);

			showSnackbar(
				user.is_superuser ? "Logged in as admin" : "Logged in successfully",
				"success",
			);
			onClose();

			// Clear form
			setEmail("");
			setPassword("");
		} catch (err) {
			setError(err instanceof Error ? err.message : "Login failed");
		} finally {
			setLoading(false);
		}
	};

	const handleKeyDown = (event: React.KeyboardEvent) => {
		if (event.key === "Enter" && !loading) {
			if (email && password) {
				handleLogin();
			}
		}
	};

	const handleClose = () => {
		if (!loading) {
			setError(null);
			setEmail("");
			setPassword("");
			onClose();
		}
	};

	return (
		<Dialog open={open} onClose={handleClose} maxWidth="xs" fullWidth>
			<DialogTitle>Login to ZnDraw</DialogTitle>
			<DialogContent>
				<Box sx={{ pt: 1, display: "flex", flexDirection: "column", gap: 2 }}>
					<Typography variant="body2" color="text.secondary">
						Login with your email and password, or continue as guest.
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
						autoComplete="current-password"
					/>
				</Box>
			</DialogContent>
			<DialogActions sx={{ px: 3, pb: 2, flexDirection: "column", gap: 1 }}>
				<Button
					onClick={handleLogin}
					disabled={loading || !email || !password}
					variant="contained"
					fullWidth
				>
					{loading ? "Logging in..." : "Login"}
				</Button>
				<Button
					onClick={handleLoginAsGuest}
					disabled={loading}
					variant="outlined"
					fullWidth
				>
					Continue as Guest
				</Button>
			</DialogActions>
		</Dialog>
	);
}
