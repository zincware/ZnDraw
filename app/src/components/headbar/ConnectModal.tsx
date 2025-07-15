import {
	Button,
	Container,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
} from "@mui/material";
import { FaRegClipboard } from "react-icons/fa";
import Markdown from "react-markdown";

function getServerUrl() {
	const { protocol, host } = window.location;
	const basePath = import.meta.env.BASE_URL || "/";
	return `${protocol}//${host}${basePath}`;
}

export function ConnectModal({ show, onHide, token }) {
	// const url = window.location.href.replace(/\/$/, "");
	// for testing the socketio url is different and hard coded
	const serverUrl = getServerUrl();
	// const serverUrl = "http://localhost:1235";

	const pythonCode = `from zndraw import ZnDraw

vis = ZnDraw(url="${serverUrl}", token="${token}")`;

	const helpMD = `
You can connect a Python kernel to this ZnDraw session using

\`\`\`python
${pythonCode}
\`\`\`
  `;

	const copyToClipboard = () => {
		navigator.clipboard.writeText(pythonCode);
		onHide();
	};

	return (
		<Dialog
			open={show}
			onClose={onHide}
			aria-labelledby="contained-modal-title-vcenter"
			maxWidth="lg"
			fullWidth
		>
			<DialogTitle id="contained-modal-title-vcenter">
				Python Client
			</DialogTitle>
			<DialogContent>
				<Container>
					<Markdown>{helpMD}</Markdown>
				</Container>
			</DialogContent>
			<DialogActions>
				<Button onClick={copyToClipboard}>
					<FaRegClipboard /> Copy to Clipboard
				</Button>
				<Button onClick={onHide}>Close</Button>
			</DialogActions>
		</Dialog>
	);
}
