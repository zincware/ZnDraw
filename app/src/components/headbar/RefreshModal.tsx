import {
	Button,
	Container,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
} from "@mui/material";

function getServerUrl() {
	const { protocol, host } = window.location;
	const basePath = import.meta.env.BASE_URL || "/";
	return `${protocol}//${host}${basePath}`;
}

export function RefreshModal({ show, onHide, room, token }) {
	const serverUrl = getServerUrl();
	const urlWithRoom = token ? `${serverUrl}token/${room}/${token}` : `${serverUrl}token/${room}`;
	const resetURL = `${serverUrl}reset`;

	return (
		<Dialog
			open={show}
			onClose={onHide}
			aria-labelledby="contained-modal-title-vcenter"
			maxWidth="lg"
			fullWidth
		>
			<DialogTitle id="contained-modal-title-vcenter">Reload Scene</DialogTitle>
			<DialogContent>
				<Container>
					<p>
						You are about the reload the scene. The current scene will be
						available as long as some client is connected, otherwise the data
						will be deleted. You can access the scene by visiting <br />
						<a href={urlWithRoom}>{urlWithRoom}</a>
					</p>
					<Button
						variant="outlined"
						onClick={() => navigator.clipboard.writeText(urlWithRoom)}
					>
						copy URL to clipboard
					</Button>
				</Container>
			</DialogContent>
			<DialogActions>
				<Button onClick={onHide}>Cancel</Button>
				<Button href={resetURL}>Create new Scene</Button>
			</DialogActions>
		</Dialog>
	);
}
