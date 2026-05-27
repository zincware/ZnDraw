import AppBar from "@mui/material/AppBar";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Chip from "@mui/material/Chip";
import Container from "@mui/material/Container";
import IconButton from "@mui/material/IconButton";
import List from "@mui/material/List";
import ListItem from "@mui/material/ListItem";
import ListItemText from "@mui/material/ListItemText";
import Stack from "@mui/material/Stack";
import TextField from "@mui/material/TextField";
import Toolbar from "@mui/material/Toolbar";
import Typography from "@mui/material/Typography";
import ArrowBackIcon from "@mui/icons-material/ArrowBack";
import { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";
import { createGroup, type Group, listGroups } from "../myapi/client";

export default function GroupsPage() {
	const navigate = useNavigate();
	const [groups, setGroups] = useState<Group[]>([]);
	const [name, setName] = useState("");
	const [description, setDescription] = useState("");
	const [loading, setLoading] = useState(false);
	const [error, setError] = useState<string | null>(null);

	useEffect(() => {
		setLoading(true);
		listGroups()
			.then(setGroups)
			.catch((err) => setError(String(err)))
			.finally(() => setLoading(false));
	}, []);

	async function onCreate() {
		if (!name.trim()) return;
		setError(null);
		try {
			const g = await createGroup(name.trim(), description.trim() || undefined);
			setGroups((cur) => [...cur, g]);
			setName("");
			setDescription("");
		} catch (err) {
			setError(String(err));
		}
	}

	return (
		<Box>
			<AppBar position="static" color="default" elevation={1}>
				<Toolbar>
					<IconButton
						edge="start"
						aria-label="back"
						onClick={() => navigate("/")}
					>
						<ArrowBackIcon />
					</IconButton>
					<Typography variant="h6" sx={{ ml: 1 }}>
						Groups
					</Typography>
				</Toolbar>
			</AppBar>
			<Container maxWidth="md" sx={{ mt: 3 }}>
				<Typography variant="h5" gutterBottom>
					Create a group
				</Typography>
				<Stack direction="row" spacing={2} sx={{ mb: 3 }}>
					<TextField
						label="Name"
						value={name}
						onChange={(e) => setName(e.target.value)}
						size="small"
						inputProps={{ pattern: "[-a-zA-Z0-9_]+", maxLength: 64 }}
					/>
					<TextField
						label="Description (optional)"
						value={description}
						onChange={(e) => setDescription(e.target.value)}
						size="small"
						sx={{ flex: 1 }}
					/>
					<Button variant="contained" onClick={onCreate}>
						Create
					</Button>
				</Stack>

				{error && (
					<Typography color="error" sx={{ mb: 2 }}>
						{error}
					</Typography>
				)}
				{loading && <Typography>Loading...</Typography>}

				<Typography variant="h5" gutterBottom>
					My groups
				</Typography>
				{!loading && groups.length === 0 && (
					<Typography color="text.secondary">
						You're not in any groups yet.
					</Typography>
				)}
				<List>
					{groups.map((g) => (
						<ListItem key={g.id} divider>
							<ListItemText
								primary={g.name}
								secondary={g.description ?? undefined}
							/>
							{g.my_role && (
								<Chip
									label={g.my_role}
									size="small"
									color={g.my_role === "admin" ? "primary" : "default"}
								/>
							)}
						</ListItem>
					))}
				</List>
			</Container>
		</Box>
	);
}
