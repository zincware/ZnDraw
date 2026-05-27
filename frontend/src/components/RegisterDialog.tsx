import RefreshIcon from "@mui/icons-material/Refresh";
import {
	Alert,
	Box,
	Button,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
	IconButton,
	InputAdornment,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import type React from "react";
import { useEffect, useRef, useState } from "react";
import { connectWithAuth } from "../socket";
import { useAppStore } from "../store";
import { registerUser } from "../utils/auth";

const DISPLAY_NAME_RE = /^[a-z][a-z0-9-]{2,63}$/;

async function fetchSuggestion(): Promise<string> {
	const resp = await fetch("/v1/users/available-display-name");
	if (!resp.ok) throw new Error("Failed to suggest a display name");
	const data = await resp.json();
	return data.display_name as string;
}

interface RegisterDialogProps {
	open: boolean;
	onClose: () => void;
}

export default function RegisterDialog({ open, onClose }: RegisterDialogProps) {
	const [email, setEmail] = useState("");
	const [displayName, setDisplayName] = useState("");
	const [password, setPassword] = useState("");
	const [passwordConfirm, setPasswordConfirm] = useState("");
	const [error, setError] = useState<string | null>(null);
	const [loading, setLoading] = useState(false);
	const userTouchedRef = useRef(false);

	// Use individual selectors to prevent unnecessary re-renders
	const setUser = useAppStore((state) => state.setUser);
	const showSnackbar = useAppStore((state) => state.showSnackbar);
	const userDisplayName = useAppStore(
		(state) => state.user?.display_name ?? null,
	);

	useEffect(() => {
		if (!open) {
			userTouchedRef.current = false;
			return;
		}
		let cancelled = false;
		fetchSuggestion()
			.then((suggested) => {
				if (cancelled || userTouchedRef.current) return;
				setDisplayName(suggested);
			})
			.catch(() => {});
		return () => {
			cancelled = true;
		};
	}, [open]);

	const regenerate = async () => {
		try {
			const suggested = await fetchSuggestion();
			userTouchedRef.current = false;
			setDisplayName(suggested);
		} catch {
			// Silent: keep the user's typed value.
		}
	};

	const handleRegister = async () => {
		setError(null);

		if (!email.trim()) {
			setError("Email is required");
			return;
		}

		if (!DISPLAY_NAME_RE.test(displayName)) {
			setError(
				"Display name must be 3–64 chars, lowercase, digits and hyphens, starting with a letter",
			);
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
			await registerUser(email, password, displayName);
			const { user } = await connectWithAuth();

			setUser(user);

			showSnackbar(`Registered as ${user.display_name}`, "success");
			onClose();

			// Clear form
			setEmail("");
			setDisplayName("");
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
			if (email && displayName && password && passwordConfirm) {
				handleRegister();
			}
		}
	};

	const handleClose = () => {
		if (!loading) {
			setError(null);
			setEmail("");
			setDisplayName("");
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
						Current temporary name: <strong>{userDisplayName}</strong>
					</Typography>
					<Typography variant="body2" color="text.secondary">
						Pick a display name (or accept the suggestion), then enter an email
						and password.
					</Typography>

					{error && (
						<Alert severity="error" onClose={() => setError(null)}>
							{error}
						</Alert>
					)}

					<TextField
						label="Display name"
						value={displayName}
						onChange={(e) => {
							userTouchedRef.current = true;
							setDisplayName(e.target.value);
						}}
						onKeyDown={handleKeyDown}
						disabled={loading}
						fullWidth
						autoComplete="off"
						InputProps={{
							endAdornment: (
								<InputAdornment position="end">
									<Tooltip title="Regenerate suggestion">
										<IconButton
											size="small"
											onClick={regenerate}
											disabled={loading}
										>
											<RefreshIcon fontSize="small" />
										</IconButton>
									</Tooltip>
								</InputAdornment>
							),
						}}
					/>

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
					disabled={
						loading || !email || !displayName || !password || !passwordConfirm
					}
					variant="contained"
				>
					{loading ? "Registering..." : "Register"}
				</Button>
			</DialogActions>
		</Dialog>
	);
}
