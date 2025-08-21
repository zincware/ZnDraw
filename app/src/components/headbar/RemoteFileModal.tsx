import {
	Button,
	Container,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
	FormControl,
	Grid,
	InputLabel,
	MenuItem,
	Select,
	TextField,
} from "@mui/material";
import { useEffect, useRef, useState } from "react";
import { MdOutlineAutoGraph } from "react-icons/md";
import { SiMoleculer } from "react-icons/si";
import Markdown from "react-markdown";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import {
	oneDark,
	oneLight,
} from "react-syntax-highlighter/dist/esm/styles/prism";
import { socket } from "../../socket";

export function RemoteFileModal({ show, onHide, colorMode }: any) {
	const [availableNodes, setAvailableNodes] = useState<string[]>([]);
	const [selectedNodeAttributes, setSelectedNodeAttributes] = useState<
		string[][]
	>([]);
	const [loading, setLoading] = useState(false);
	const [selectedNode, setSelectedNode] = useState("");
	const remoteRepoRef = useRef("");
	const remoteRevRef = useRef("");

	function fetchAvailableNodes() {
		setLoading(true);
		socket.emit(
			"zntrack:list-stages",
			{ remote: remoteRepoRef.current, rev: remoteRevRef.current },
			(data: string[]) => {
				setAvailableNodes(data);
				setLoading(false);
			},
		);
	}

	useEffect(() => {
		if (selectedNode) {
			setLoading(true);
			socket.emit(
				"zntrack:inspect-stage",
				{
					remote: remoteRepoRef.current,
					rev: remoteRevRef.current,
					name: selectedNode,
				},
				(data: string[][]) => {
					setLoading(false);
					setSelectedNodeAttributes(data);
				},
			);
		}
	}, [selectedNode]);

	function loadFrames(attribute: string) {
		socket.emit("zntrack:load-frames", {
			remote: remoteRepoRef.current,
			rev: remoteRevRef.current,
			name: selectedNode,
			attribute: attribute,
		});
		onHide();
	}
	function loadFigure(attribute: string) {
		socket.emit("zntrack:load-figure", {
			remote: remoteRepoRef.current,
			rev: remoteRevRef.current,
			name: selectedNode,
			attribute: attribute,
		});
		onHide();
	}

	return (
		<Dialog
			open={show}
			onClose={onHide}
			aria-labelledby="contained-modal-title-vcenter"
			maxWidth="lg"
			fullWidth
		>
			<DialogTitle id="contained-modal-title-vcenter">
				Load Data via DVC and ZnTrack
			</DialogTitle>
			<DialogContent>
				<Container>
					<TextField
						label="Repo"
						placeholder="Repository path"
						variant="outlined"
						fullWidth
						sX={{ mb: 3 }}
						onChange={(e) => {
							remoteRepoRef.current = e.target.value;
						}}
					/>
					<TextField
						label="Rev"
						placeholder="Revision"
						variant="outlined"
						fullWidth
						sX={{ mb: 3 }}
						onChange={(e) => {
							remoteRevRef.current = e.target.value;
						}}
					/>
					<Button onClick={fetchAvailableNodes}>Show available data</Button>
					<hr />
					{loading && <p>Loading...</p>}

					{availableNodes.length > 0 && (
						<FormControl fullWidth sx={{ mt: 2 }}>
							<InputLabel>Choose...</InputLabel>
							<Select
								value={selectedNode}
								onChange={(e) => setSelectedNode(e.target.value)}
								label="Choose..."
							>
								{availableNodes.map((node) => (
									<MenuItem key={node} value={node}>
										{node}
									</MenuItem>
								))}
							</Select>
						</FormControl>
					)}

					{selectedNodeAttributes.length > 0 &&
						selectedNodeAttributes.map((attr, index) => (
							<Grid container key={index} justifyContent="center">
								<Grid item xs={12}>
									{" "}
									<Markdown
										children={`~~~python\n${attr.join(": ")}`}
										components={{
											code(props) {
												const { children, className, node, ...rest } = props;
												const match = /language-(\w+)/.exec(className || "");
												return match ? (
													<SyntaxHighlighter
														{...rest}
														PreTag="div"
														language={match[1]}
														style={colorMode === "light" ? oneLight : oneDark}
													>
														{String(children).replace(/\n$/, "")}
													</SyntaxHighlighter>
												) : (
													<code {...rest} className={className}>
														{children}
													</code>
												);
											},
										}}
									/>{" "}
								</Grid>
								<Grid item md="auto">
									<Button
										variant="contained"
										onClick={() => loadFrames(attr[0])}
									>
										<SiMoleculer /> Load Frames
									</Button>
								</Grid>
								<Grid item md="auto">
									<Button
										variant="contained"
										onClick={() => loadFigure(attr[0])}
									>
										<MdOutlineAutoGraph /> Load Figure
									</Button>
								</Grid>
							</Grid>
						))}
				</Container>
			</DialogContent>
			<DialogActions>
				<Button onClick={onHide}>Cancel</Button>
			</DialogActions>
		</Dialog>
	);
}
