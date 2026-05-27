import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import DeleteIcon from "@mui/icons-material/Delete";
import Button from "@mui/material/Button";
import Dialog from "@mui/material/Dialog";
import DialogActions from "@mui/material/DialogActions";
import DialogContent from "@mui/material/DialogContent";
import DialogTitle from "@mui/material/DialogTitle";
import FormControl from "@mui/material/FormControl";
import IconButton from "@mui/material/IconButton";
import InputLabel from "@mui/material/InputLabel";
import List from "@mui/material/List";
import ListItem from "@mui/material/ListItem";
import ListItemText from "@mui/material/ListItemText";
import MenuItem from "@mui/material/MenuItem";
import Select, { type SelectChangeEvent } from "@mui/material/Select";
import Stack from "@mui/material/Stack";
import Typography from "@mui/material/Typography";
import { useEffect, useState } from "react";
import {
	createShareLink,
	listShareLinks,
	revokeShareLink,
	type ShareAccess,
	type ShareLink,
} from "../myapi/client";

interface Props {
	roomId: string;
	open: boolean;
	onClose: () => void;
}

export default function ShareDialog({ roomId, open, onClose }: Props) {
	const [links, setLinks] = useState<ShareLink[]>([]);
	const [access, setAccess] = useState<ShareAccess>("view");
	const [loading, setLoading] = useState(false);
	const [error, setError] = useState<string | null>(null);

	useEffect(() => {
		if (!open) return;
		let cancelled = false;
		setLoading(true);
		setError(null);
		listShareLinks(roomId)
			.then((items) => {
				if (!cancelled) setLinks(items);
			})
			.catch((err) => {
				if (!cancelled) setError(String(err));
			})
			.finally(() => {
				if (!cancelled) setLoading(false);
			});
		return () => {
			cancelled = true;
		};
	}, [open, roomId]);

	async function onCreate() {
		try {
			const link = await createShareLink(roomId, access);
			setLinks((cur) => [...cur, link]);
		} catch (err) {
			setError(String(err));
		}
	}

	async function onRevoke(id: string) {
		try {
			await revokeShareLink(roomId, id);
			setLinks((cur) => cur.filter((l) => l.id !== id));
		} catch (err) {
			setError(String(err));
		}
	}

	function shareUrl(token: string): string {
		const url = new URL(window.location.href);
		url.searchParams.set("share", token);
		return url.toString();
	}

	return (
		<Dialog open={open} onClose={onClose} fullWidth maxWidth="sm">
			<DialogTitle>Share Room</DialogTitle>
			<DialogContent>
				<Stack direction="row" spacing={2} sx={{ mb: 2, alignItems: "center" }}>
					<FormControl size="small" sx={{ minWidth: 140 }}>
						<InputLabel>Access</InputLabel>
						<Select
							value={access}
							label="Access"
							onChange={(e: SelectChangeEvent) =>
								setAccess(e.target.value as ShareAccess)
							}
						>
							<MenuItem value="view">View only</MenuItem>
							<MenuItem value="edit">View + edit</MenuItem>
						</Select>
					</FormControl>
					<Button variant="contained" onClick={onCreate}>
						Create link
					</Button>
				</Stack>

				{error && (
					<Typography color="error" variant="body2" sx={{ mb: 1 }}>
						{error}
					</Typography>
				)}
				{loading && <Typography variant="body2">Loading…</Typography>}

				{!loading && links.length === 0 && (
					<Typography variant="body2" color="text.secondary">
						No active share links.
					</Typography>
				)}

				<List>
					{links.map((l) => (
						<ListItem
							key={l.id}
							secondaryAction={
								<>
									<IconButton
										aria-label="copy"
										onClick={() =>
											navigator.clipboard.writeText(shareUrl(l.token))
										}
									>
										<ContentCopyIcon />
									</IconButton>
									<IconButton
										aria-label="revoke"
										onClick={() => onRevoke(l.id)}
									>
										<DeleteIcon />
									</IconButton>
								</>
							}
						>
							<ListItemText
								primary={l.access === "edit" ? "View + edit" : "View only"}
								secondary={shareUrl(l.token)}
							/>
						</ListItem>
					))}
				</List>
			</DialogContent>
			<DialogActions>
				<Button onClick={onClose}>Close</Button>
			</DialogActions>
		</Dialog>
	);
}
