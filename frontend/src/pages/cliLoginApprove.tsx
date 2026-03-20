import CheckCircleOutlineIcon from "@mui/icons-material/CheckCircleOutline";
import DoNotDisturbIcon from "@mui/icons-material/DoNotDisturb";
import TerminalIcon from "@mui/icons-material/Terminal";
import {
	Alert,
	Box,
	Button,
	Card,
	CardContent,
	Chip,
	CircularProgress,
	Container,
	Stack,
	Typography,
} from "@mui/material";
import axios from "axios";
import { useEffect, useState } from "react";
import { useSearchParams } from "react-router-dom";
import { type UserInfo, acquireToken, getToken } from "../utils/auth";

type Status =
	| "loading"
	| "ready"
	| "approving"
	| "approved"
	| "denied"
	| "error";

export default function CliLoginApprovePage() {
	const [searchParams] = useSearchParams();
	const code = searchParams.get("code");
	const [status, setStatus] = useState<Status>("loading");
	const [error, setError] = useState<string | null>(null);
	const [user, setUser] = useState<UserInfo | null>(null);

	useEffect(() => {
		const init = async () => {
			try {
				const result = await acquireToken();
				setUser(result.user);
				setStatus("ready");
			} catch {
				setError("Failed to authenticate. Please log in first.");
				setStatus("error");
			}
		};
		init();
	}, []);

	const handleApprove = async () => {
		if (!code) return;
		setStatus("approving");
		setError(null);
		try {
			await axios.patch(
				`/v1/auth/cli-login/${encodeURIComponent(code)}`,
				null,
				{
					headers: { Authorization: `Bearer ${getToken()}` },
				},
			);
			setStatus("approved");
		} catch (err) {
			if (axios.isAxiosError(err) && err.response?.status === 410) {
				setError(
					"This login challenge has expired. Please try again from the CLI.",
				);
			} else {
				setError("Failed to approve. Please try again.");
			}
			setStatus("error");
		}
	};

	const handleDeny = async () => {
		if (!code) return;
		setStatus("approving");
		setError(null);
		try {
			await axios.delete(`/v1/auth/cli-login/${encodeURIComponent(code)}`, {
				headers: { Authorization: `Bearer ${getToken()}` },
			});
			setStatus("denied");
		} catch {
			setError("Failed to deny. The challenge may have already expired.");
			setStatus("error");
		}
	};

	if (!code) {
		return (
			<Container maxWidth="sm" sx={{ mt: 8 }}>
				<Alert severity="error">No login code provided in URL.</Alert>
			</Container>
		);
	}

	return (
		<Container maxWidth="sm" sx={{ mt: 8 }}>
			<Card variant="outlined">
				<CardContent>
					<Stack spacing={3} alignItems="center" sx={{ py: 2 }}>
						<TerminalIcon sx={{ fontSize: 48, color: "text.secondary" }} />
						<Typography variant="h5">Authorize CLI Access</Typography>
						<Typography
							variant="body2"
							color="text.secondary"
							textAlign="center"
						>
							A CLI client is requesting access to your account. Verify the code
							matches what you see in your terminal.
						</Typography>

						{user && (
							<Stack direction="row" spacing={1} alignItems="center">
								<Typography variant="body2" color="text.secondary">
									Logged in as
								</Typography>
								<Typography variant="body2" fontWeight={600}>
									{user.email}
								</Typography>
								<Chip
									label={user.is_superuser ? "admin" : "user"}
									size="small"
									color={user.is_superuser ? "warning" : "default"}
									variant="outlined"
								/>
							</Stack>
						)}

						<Box
							sx={{
								px: 4,
								py: 2,
								borderRadius: 2,
								bgcolor: "action.hover",
								fontFamily: "monospace",
								fontSize: "1.5rem",
								fontWeight: 700,
								letterSpacing: "0.15em",
							}}
						>
							{code}
						</Box>

						{error && (
							<Alert severity="error" sx={{ width: "100%" }}>
								{error}
							</Alert>
						)}

						{status === "loading" && <CircularProgress size={24} />}

						{status === "approved" && (
							<Alert
								severity="success"
								icon={<CheckCircleOutlineIcon />}
								sx={{ width: "100%" }}
							>
								CLI login approved. You can close this tab.
							</Alert>
						)}

						{status === "denied" && (
							<Alert
								severity="info"
								icon={<DoNotDisturbIcon />}
								sx={{ width: "100%" }}
							>
								CLI login denied. You can close this tab.
							</Alert>
						)}

						{(status === "ready" || status === "error") && (
							<Stack direction="row" spacing={2}>
								<Button
									variant="contained"
									color="primary"
									size="large"
									onClick={handleApprove}
									startIcon={<CheckCircleOutlineIcon />}
								>
									Approve
								</Button>
								<Button
									variant="outlined"
									color="error"
									size="large"
									onClick={handleDeny}
									startIcon={<DoNotDisturbIcon />}
								>
									Deny
								</Button>
							</Stack>
						)}

						{status === "approving" && <CircularProgress size={24} />}
					</Stack>
				</CardContent>
			</Card>
		</Container>
	);
}
