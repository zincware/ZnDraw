import DeleteIcon from "@mui/icons-material/Delete";
import PersonIcon from "@mui/icons-material/Person";
import SecurityIcon from "@mui/icons-material/Security";
import {
	Alert,
	Box,
	Button,
	Chip,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
	IconButton,
	Paper,
	Table,
	TableBody,
	TableCell,
	TableContainer,
	TableHead,
	TableRow,
	Tooltip,
	Typography,
} from "@mui/material";
import { useState, useEffect } from "react";
import {
	type AdminUser,
	deleteAdminUser,
	listAdminUsers,
	updateAdminUser,
} from "../myapi/client";

interface AdminPanelProps {
	open: boolean;
	onClose: () => void;
}

export default function AdminPanel({ open, onClose }: AdminPanelProps) {
	const [users, setUsers] = useState<AdminUser[]>([]);
	const [loading, setLoading] = useState(false);
	const [error, setError] = useState<string | null>(null);

	const fetchUsers = async () => {
		setLoading(true);
		setError(null);
		try {
			const data = await listAdminUsers();
			setUsers(data.items);
		} catch (err) {
			setError(err instanceof Error ? err.message : "Failed to load users");
		} finally {
			setLoading(false);
		}
	};

	useEffect(() => {
		if (open) {
			fetchUsers();
		}
	}, [open]);

	const handleToggleSuperuser = async (user: AdminUser) => {
		try {
			await updateAdminUser(user.id, { is_superuser: !user.is_superuser });
			await fetchUsers();
		} catch (err) {
			setError(
				err instanceof Error ? err.message : "Failed to update user role",
			);
		}
	};

	const handleDelete = async (user: AdminUser) => {
		if (
			!confirm(
				`Are you sure you want to delete user "${user.email}"? This cannot be undone.`,
			)
		) {
			return;
		}

		try {
			await deleteAdminUser(user.id);
			await fetchUsers();
		} catch (err) {
			setError(err instanceof Error ? err.message : "Failed to delete user");
		}
	};

	return (
		<Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
			<DialogTitle>User Management</DialogTitle>
			<DialogContent>
				<Box sx={{ mt: 1 }}>
					{error && (
						<Alert
							severity="error"
							onClose={() => setError(null)}
							sx={{ mb: 2 }}
						>
							{error}
						</Alert>
					)}

					{loading ? (
						<Typography>Loading users...</Typography>
					) : (
						<TableContainer component={Paper}>
							<Table>
								<TableHead>
									<TableRow>
										<TableCell>Email</TableCell>
										<TableCell>Role</TableCell>
										<TableCell align="right">Actions</TableCell>
									</TableRow>
								</TableHead>
								<TableBody>
									{users.map((user) => (
										<TableRow key={user.id}>
											<TableCell>{user.email}</TableCell>
											<TableCell>
												{user.is_superuser ? (
													<Chip label="Admin" color="error" size="small" />
												) : (
													<Chip label="User" color="primary" size="small" />
												)}
											</TableCell>
											<TableCell align="right">
												<Box
													sx={{
														display: "flex",
														gap: 0.5,
														justifyContent: "flex-end",
													}}
												>
													{user.is_superuser ? (
														<Tooltip title="Demote to user">
															<IconButton
																size="small"
																onClick={() => handleToggleSuperuser(user)}
																color="warning"
															>
																<PersonIcon />
															</IconButton>
														</Tooltip>
													) : (
														<Tooltip title="Promote to admin">
															<IconButton
																size="small"
																onClick={() => handleToggleSuperuser(user)}
																color="primary"
															>
																<SecurityIcon />
															</IconButton>
														</Tooltip>
													)}

													<Tooltip title="Delete user">
														<IconButton
															size="small"
															onClick={() => handleDelete(user)}
															color="error"
														>
															<DeleteIcon />
														</IconButton>
													</Tooltip>
												</Box>
											</TableCell>
										</TableRow>
									))}
								</TableBody>
							</Table>
						</TableContainer>
					)}
				</Box>
			</DialogContent>
			<DialogActions>
				<Button onClick={fetchUsers} disabled={loading}>
					Refresh
				</Button>
				<Button onClick={onClose}>Close</Button>
			</DialogActions>
		</Dialog>
	);
}
