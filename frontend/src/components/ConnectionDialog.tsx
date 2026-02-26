import ContentCopyIcon from "@mui/icons-material/ContentCopy";
import Box from "@mui/material/Box";
import Button from "@mui/material/Button";
import Dialog from "@mui/material/Dialog";
import DialogContent from "@mui/material/DialogContent";
import DialogTitle from "@mui/material/DialogTitle";
import Typography from "@mui/material/Typography";
import { useColorScheme } from "@mui/material/styles";
import { useParams } from "react-router-dom";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import { oneLight } from "react-syntax-highlighter/dist/esm/styles/prism";
import { oneDark } from "react-syntax-highlighter/dist/esm/styles/prism";
import { useAppStore } from "../store";

interface ConnectionDialogProps {
	open: boolean;
	onClose: () => void;
}

const ConnectionDialog = ({ open, onClose }: ConnectionDialogProps) => {
	const { roomId } = useParams<{ roomId: string }>();
	const userName = useAppStore((state) => state.user?.email ?? null);
	const sessionId = useAppStore((state) => state.sessionId);
	const { mode } = useColorScheme();

	const connectionCode = `from zndraw import ZnDraw\n\nvis = ZnDraw(\n  url="${window.location.origin}/",\n  room="${roomId}",\n  user="${userName}"\n)`;

	const sessionCode = sessionId
		? `ses = vis.sessions["${sessionId}"]`
		: "# No active session";

	const handleCopyCode = () => {
		navigator.clipboard.writeText(connectionCode);
	};

	const handleCopySessionCode = () => {
		navigator.clipboard.writeText(sessionCode);
	};

	return (
		<Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
			<DialogTitle>Connect to this ZnDraw session</DialogTitle>
			<DialogContent>
				<Typography variant="body2" sx={{ mb: 2 }}>
					You can connect a Python kernel to this ZnDraw session using:
				</Typography>
				<Box sx={{ position: "relative" }}>
					<SyntaxHighlighter
						language="python"
						style={mode === "dark" ? oneDark : oneLight}
						customStyle={{
							margin: 0,
							borderRadius: "4px",
						}}
					>
						{connectionCode}
					</SyntaxHighlighter>
					<Button
						startIcon={<ContentCopyIcon />}
						onClick={handleCopyCode}
						sx={{
							mt: 1,
						}}
						variant="outlined"
						size="small"
					>
						Copy to clipboard
					</Button>
				</Box>
				<Typography variant="body2" sx={{ mt: 3, mb: 2 }}>
					Access this browser session:
				</Typography>
				<Box sx={{ position: "relative" }}>
					<SyntaxHighlighter
						language="python"
						style={mode === "dark" ? oneDark : oneLight}
						customStyle={{
							margin: 0,
							borderRadius: "4px",
						}}
					>
						{sessionCode}
					</SyntaxHighlighter>
					<Button
						startIcon={<ContentCopyIcon />}
						onClick={handleCopySessionCode}
						sx={{
							mt: 1,
						}}
						variant="outlined"
						size="small"
					>
						Copy to clipboard
					</Button>
				</Box>
			</DialogContent>
		</Dialog>
	);
};

export default ConnectionDialog;
