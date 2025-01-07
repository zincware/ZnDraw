import type React from "react";
import {
	SetStateAction,
	forwardRef,
	useEffect,
	useMemo,
	useRef,
	useState,
} from "react";
import {
	Button,
	Card,
	Col,
	Container,
	Form,
	InputGroup,
	Modal,
	Nav,
	Navbar,
	Row,
	ToggleButton,
} from "react-bootstrap";
import { version } from "../../package.json";
import { MdOutlineAutoGraph } from "react-icons/md";
import { SiMoleculer } from "react-icons/si";
import Select from "react-select";
import rehypeKatex from "rehype-katex";
import rehypeRaw from "rehype-raw";
import remarkBreaks from "remark-breaks";
import remarkGfm from "remark-gfm";
import remarkMath from "remark-math";
import "katex/dist/katex.min.css"; // `rehype-katex` does not import the CSS for you
import {
	FaCloudDownloadAlt,
	FaCode,
	FaDownload,
	FaFilm,
	FaHandSparkles,
	FaLock,
	FaLockOpen,
	FaMoon,
	FaPlus,
	FaRegClipboard,
	FaRocket,
	FaSave,
	FaSun,
	FaTerminal,
	FaUpload,
} from "react-icons/fa";
import { GrHelpBook } from "react-icons/gr";
import Markdown from "react-markdown";
import { Rnd } from "react-rnd";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import {
	oneDark,
	oneLight,
} from "react-syntax-highlighter/dist/esm/styles/prism";
import * as znsocket from "znsocket";
import { client, socket } from "../socket";
import { BtnTooltip } from "./tooltips";

import {
	FaArrowRotateRight,
	FaFileCirclePlus,
	FaPencil,
} from "react-icons/fa6";
import { MdAddChart, MdExitToApp } from "react-icons/md";
import { TbPlugConnected } from "react-icons/tb";

function getServerUrl() {
	const { protocol, host } = window.location;
	const basePath = import.meta.env.BASE_URL || "/";
	return `${protocol}//${host}${basePath}`;
}

function ConsoleWindow({
	messages,
	setConsoleShow,
	token,
	setMessages,
	colorMode,
	step,
	selection,
}: {
	messages: string[];
	setConsoleShow: any;
	token: string;
	setMessages: any;
	colorMode: string;
	step: number;
	selection: Set<number>;
}) {
	const [showTime, setShowTime] = useState(false);
	const [chatInput, setChatInput] = useState<object>({});
	const [showDropdown, setShowDropdown] = useState(false);
	const chatInputRef = useRef(null);
	const scrollRef = useRef(null);

	useEffect(() => {
		if (scrollRef.current) {
			scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
		}
	}, [messages]);

	const handleChatInputChange = (e: any) => {
		// log the selectionStart
		setChatInput({
			msg: e.target.value,
			time: new Date().toLocaleTimeString(),
		});
		if (e.target.value.endsWith("!!")) {
			setShowDropdown(true);
		} else {
			setShowDropdown(false);
		}
	};

	const handleSendMessage = () => {
		setMessages([...messages, chatInput]);
		if (chatInputRef.current) {
			chatInputRef.current.value = "";
		}
	};

	const handleKeyPress = (e) => {
		if (e.key === "Enter" && !e.shiftKey) {
			e.preventDefault();
			handleSendMessage();
		}
	};

	const [isEditing, setIsEditing] = useState(null);
	const [tempMsg, setTempMsg] = useState("");

	const handleEdit = (key, msg) => {
		setIsEditing(key);
		setTempMsg(msg);
	};

	const handleSave = (key: number) => {
		// TODO: do we want to have an edited state?
		const newMessages = messages.map((line, idx) =>
			idx === key ? { ...line, msg: tempMsg } : line,
		);
		setMessages(newMessages);
		setIsEditing(null); // Exit editing mode
	};

	const handleEditKeyPress = (e, idx) => {
		if (e.key === "Enter" && !e.shiftKey) {
			e.preventDefault(); // Prevents adding a new line
			handleSave(idx); // Calls the save function
		}
	};

	return (
		<>
			<Rnd
				default={{
					x: window.innerWidth / 2 - 400,
					y: window.innerHeight / 2 - 300,
					width: 380,
					height: 280,
				}}
				minHeight={200}
				minWidth={200}
				style={{
					zIndex: 1000,
					padding: 0,
					margin: 0,
				}}
			>
				<Card
					style={{
						margin: 0,
						padding: 0,
						width: "100%",
						height: "100%",
					}}
					// ref={cardRef}
				>
					<Card.Header className="d-flex justify-content-between align-items-center">
						<Card.Title>Chat</Card.Title>
						<div className="d-flex align-items-center">
							<Form.Check
								type="switch"
								id="show-time-switch"
								label="Show Time"
								checked={showTime}
								onChange={() => {
									setShowTime(!showTime);
								}}
								className="me-2"
							/>
							<Button variant="close" onClick={() => setConsoleShow(false)} />
						</div>
					</Card.Header>

					{/* Message Body with Optional Timestamp */}
					<Card.Body className="text-start overflow-y-auto" ref={scrollRef}>
						{messages.map((line, idx) => (
							<div key={idx} className="mb-2">
								{/* Row for timestamp and edit icons */}
								<div className="d-flex justify-content-between align-items-center">
									{showTime && (
										<div className="d-flex align-items-center">
											<span className="text-muted">{line.time}</span>
											{isEditing === idx ? (
												<FaSave
													onClick={() => handleSave(idx)}
													className="text-muted ms-2"
													style={{ cursor: "pointer" }}
												/>
											) : (
												<FaPencil
													onClick={() => handleEdit(idx, line.msg)}
													className="text-muted ms-2"
													style={{ cursor: "pointer" }}
												/>
											)}
										</div>
									)}
								</div>

								{/* Row for message content or editable input */}
								<div>
									{isEditing !== idx ? (
										<Markdown
											remarkPlugins={[remarkMath, remarkGfm, remarkBreaks]}
											rehypePlugins={[rehypeKatex, rehypeRaw]} // security risk with rehypeRaw?
											children={line.msg}
											components={{
												code(props) {
													const { children, className, ...rest } = props;
													const match = /language-(\w+)/.exec(className || "");
													return match ? (
														<SyntaxHighlighter
															{...rest}
															PreTag="div"
															children={String(children).replace(/\n$/, "")}
															language={match[1]}
															style={colorMode === "light" ? oneLight : oneDark}
														/>
													) : (
														<code {...rest} className={className}>
															{children}
														</code>
													);
												},
											}}
										/>
									) : (
										<InputGroup>
											<Form.Control
												as="textarea"
												value={tempMsg}
												rows={4}
												onChange={(e) => setTempMsg(e.target.value)}
												onKeyDown={(e) => handleEditKeyPress(e, idx)}
											/>
										</InputGroup>
									)}
								</div>
							</div>
						))}
					</Card.Body>

					<Card.Footer>
						<InputGroup>
							<Form.Control
								as="textarea"
								rows={1}
								placeholder="Type a message..."
								onInput={handleChatInputChange}
								onKeyDown={handleKeyPress}
								ref={chatInputRef}
							/>
							<Button variant="primary" onClick={handleSendMessage}>
								Send
							</Button>
						</InputGroup>
					</Card.Footer>
				</Card>
			</Rnd>
			{showDropdown && (
				<ChatInsertModal
					show={showDropdown}
					onHide={() => setShowDropdown(false)}
					chatInputRef={chatInputRef}
					step={step}
					selection={selection}
				/>
			)}
		</>
	);
}

function ChatInsertModal({ show, onHide, chatInputRef, step, selection }: any) {
	const options = [
		{ value: "step", label: "step" },
		{ value: "selection", label: "selection" },
	];

	// TODO: support in new / unconnected room here

	const handleSelectChange = (selectedOption: any) => {
		chatInputRef.current.value = chatInputRef.current.value.slice(0, -2);
		const basePath = `${window.location.origin}${window.location.pathname}`
			.replace(/\/+$/, "")
			.replace(/\/token\/.*/, ""); // replace everything after /token/

		if (selectedOption.value === "step") {
			chatInputRef.current.value = `${chatInputRef.current.value}[step ${step}](${basePath}/?step=${step})`;
		} else if (selectedOption.value === "selection") {
			if (selection.size === 0) {
				chatInputRef.current.value = `${chatInputRef.current.value}[${selectedOption.value}](${basePath}/?selection=null)`;
			} else {
				chatInputRef.current.value = `${chatInputRef.current.value}[${selectedOption.value}](${basePath}/?selection=${Array.from(selection)})`;
			}
		}
		// trigger the change event
		chatInputRef.current.dispatchEvent(
			new InputEvent("input", { bubbles: true }),
		);
		onHide();
	};

	return (
		<Modal
			show={show}
			aria-labelledby="contained-modal-title-vcenter"
			size="lg"
		>
			<Modal.Header closeButton>
				<Modal.Title id="contained-modal-title-vcenter">
					ZnDraw Chat Insert
				</Modal.Title>
			</Modal.Header>
			<Modal.Body>
				<Container>
					Share your current view:
					<Select
						options={options}
						onChange={handleSelectChange}
						// menuIsOpen={true}
						placeholder="Choose..."
					/>
				</Container>
			</Modal.Body>
			<Modal.Footer>
				<Button onClick={onHide}>Close</Button>
			</Modal.Footer>
		</Modal>
	);
}

function HelpModel(props: any) {
	const helpMD = `
- **play / pause**: \`keypress space\`
- **frame forwards / backwards**: \`keypress ▶\\◀\`
- **jump forwards / backwards**: \`keypress ▲\\▼\`
- **center camera around selected particle**: \`keypress C\`
- **select multiple particles**: \`keydown shift\`
- **show particle index**: \`keydown I\`
- **add bookmark at current step**: \`keypress B\`
- **jump between bookmarks**: \`shift + keypress ▶\\◀\`
- **remove bookmark**: \`shift + mouse click\`
- **toggle drawing mode**: \`keypress X\`
- **select all particles**: \`ctrl + A\`
- **remove selected particles**: \`backspace\`
- **remove line point**: \`shift + backspace\`
- **reset camera to original position**: \`keypress O\`
- **rotate camera perpendicular to the screen**: \`keypress R\` or \`keypress ctrl + R\`
- (experimental) **Enter molecule editor (none, translate, rotate)**: \`keypress E\`
`;

	return (
		<Modal {...props} aria-labelledby="contained-modal-title-vcenter" size="lg">
			<Modal.Header closeButton>
				<Modal.Title id="contained-modal-title-vcenter">
					ZnDraw Help
				</Modal.Title>
			</Modal.Header>
			<Modal.Body>
				<Container>
					<Markdown>{helpMD}</Markdown>
				</Container>
			</Modal.Body>
			<Modal.Footer>
				<Button onClick={props.onHide}>Close</Button>
			</Modal.Footer>
		</Modal>
	);
}

function ConnectModal({ show, onHide, room }) {
	// const url = window.location.href.replace(/\/$/, "");
	// for testing the socketio url is different and hard coded
	const serverUrl = getServerUrl();
	// const serverUrl = "http://localhost:1235";

	const pythonCode = `from zndraw import ZnDraw

vis = ZnDraw(url="${serverUrl}", token="${room}")`;

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
		<Modal
			show={show}
			onHide={onHide}
			aria-labelledby="contained-modal-title-vcenter"
			size="lg"
		>
			<Modal.Header closeButton>
				<Modal.Title id="contained-modal-title-vcenter">
					Python Client
				</Modal.Title>
			</Modal.Header>
			<Modal.Body>
				<Container>
					<Markdown>{helpMD}</Markdown>
				</Container>
			</Modal.Body>
			<Modal.Footer>
				<Button onClick={copyToClipboard}>
					<FaRegClipboard /> Copy to Clipboard
				</Button>
				<Button onClick={onHide}>Close</Button>
			</Modal.Footer>
		</Modal>
	);
}

function RefreshModal({ show, onHide, room }) {
	const serverUrl = getServerUrl();
	const urlWithRoom = `${serverUrl}token/${room}`;
	const resetURL = `${serverUrl}reset`;

	return (
		<Modal
			show={show}
			onHide={onHide}
			aria-labelledby="contained-modal-title-vcenter"
			size="lg"
		>
			<Modal.Header closeButton>
				<Modal.Title id="contained-modal-title-vcenter">
					Reload Scene
				</Modal.Title>
			</Modal.Header>
			<Modal.Body>
				<Container>
					<p>
						You are about the reload the scene. The current scene will be
						available as long as some client is connected, otherwise the data
						will be deleted. You can access the scene by visiting <br />
						<a href={urlWithRoom}>{urlWithRoom}</a>
					</p>
					<Button
						variant="outline-primary"
						onClick={() => navigator.clipboard.writeText(urlWithRoom)}
					>
						copy URL to clipboard
					</Button>
				</Container>
			</Modal.Body>
			<Modal.Footer>
				<Button onClick={onHide}>Cancel</Button>
				<Button href={resetURL}>Create new Scene</Button>
			</Modal.Footer>
		</Modal>
	);
}

function RemoteFileModal({ show, onHide, colorMode }: any) {
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
		<Modal
			show={show}
			onHide={onHide}
			aria-labelledby="contained-modal-title-vcenter"
			size="lg"
		>
			<Modal.Header closeButton>
				<Modal.Title id="contained-modal-title-vcenter">
					Load Data via DVC and ZnTrack
				</Modal.Title>
			</Modal.Header>
			<Modal.Body>
				<Container>
					<InputGroup className="mb-3">
						<InputGroup.Text id="basic-addon1">Repo</InputGroup.Text>
						<Form.Control
							placeholder="Repository path"
							aria-label="repositoryPath"
							aria-describedby="basic-addon1"
							onChange={(e) => {
								remoteRepoRef.current = e.target.value;
							}}
						/>
					</InputGroup>
					<InputGroup className="mb-3">
						<InputGroup.Text id="basic-addon1">Rev</InputGroup.Text>
						<Form.Control
							// defaultValue={"HEAD"}
							placeholder="Revision"
							aria-label="repository"
							aria-describedby="basic-addon1"
							onChange={(e) => {
								remoteRevRef.current = e.target.value;
							}}
						/>
					</InputGroup>
					<Button onClick={fetchAvailableNodes}>Show available data</Button>
					<hr />
					{loading && <p>Loading...</p>}

					{availableNodes.length > 0 && (
						<ul>
							<Select
								options={availableNodes.map((node) => ({
									value: node,
									label: node,
								}))}
								onChange={(selectedOption: any) => {
									setSelectedNode(selectedOption.value);
								}}
								// menuIsOpen={true}
								placeholder="Choose..."
							/>
						</ul>
					)}

					{selectedNodeAttributes.length > 0 &&
						selectedNodeAttributes.map((attr, index) => (
							<Row key={index} className="justify-content-center">
								<Col>
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
								</Col>
								<Col md="auto">
									<Button variant="primary" onClick={() => loadFrames(attr[0])}>
										<SiMoleculer /> Load Frames
									</Button>
								</Col>
								<Col md="auto">
									<Button variant="primary" onClick={() => loadFigure(attr[0])}>
										<MdOutlineAutoGraph /> Load Figure
									</Button>
								</Col>
							</Row>
						))}
				</Container>
			</Modal.Body>
			<Modal.Footer>
				<Button onClick={onHide}>Cancel</Button>
			</Modal.Footer>
		</Modal>
	);
}

interface TutorialModalProps {
	show: boolean;
	onHide: () => void;
	url: string;
}

const TutorialModal: React.FC<TutorialModalProps> = ({ show, onHide, url }) => {
	return (
		<Modal show={show} onHide={onHide} size="xl" dialogClassName="custom-modal">
			<Modal.Header closeButton>
				<Modal.Title>ZnDraw Tutorial</Modal.Title>
			</Modal.Header>
			<Modal.Body className="modal-body-custom">
				<iframe
					src={url}
					id="tutorialIframe"
					allowFullScreen
					className="iframe-custom"
				/>
			</Modal.Body>
		</Modal>
	);
};

function SiMGenButtons({
	visible,
	token,
}: {
	visible: boolean;
	token: string;
}) {
	const [disabledBtn, setDisabledBtn] = useState<boolean>(false);
	const queueRef = useRef<any>(null);

	useEffect(() => {
		const queue = new znsocket.Dict({
			client: client,
			key: `queue:${token}:modifier`,
		});
		queueRef.current = queue;

		queue.length().then((length: any) => {
			setDisabledBtn(length > 0);
		});
		queue.onRefresh(async (x: any) => {
			const length = await queue.length();
			setDisabledBtn(length > 0);
		});

		return () => {
			queue.offRefresh();
		};
	}, [token]);

	const runConnect = () => {
		if (queueRef.current) {
			queueRef.current.Connect = {};
			socket.emit("room:worker:run");
			setDisabledBtn(true);
		}
	};

	const runGenerate = () => {
		if (queueRef.current) {
			queueRef.current.SiMGenDemo = {};
			socket.emit("room:worker:run");
			setDisabledBtn(true);
		}
	};

	const createNewCanvas = () => {
		if (queueRef.current) {
			queueRef.current.NewCanvas = {};
			socket.emit("room:worker:run");
			setDisabledBtn(true);
		}
	};

	return (
		<>
			{visible && (
				<>
					<BtnTooltip text="Connect selected atoms (shift click)">
						<Button
							variant="success"
							className="mx-1"
							onClick={runConnect}
							disabled={disabledBtn}
						>
							<TbPlugConnected /> Connect
						</Button>
					</BtnTooltip>
					<BtnTooltip text="Run SiMGen molecular generation">
						<Button
							variant="success"
							className="mx-1"
							onClick={runGenerate}
							disabled={disabledBtn}
						>
							<FaRocket /> Generate
						</Button>
					</BtnTooltip>
					<BtnTooltip text="Replace scene with empty canvas and enter drawing mode">
						<Button
							variant="success"
							className="mx-1"
							onClick={createNewCanvas}
							disabled={disabledBtn}
						>
							<FaFileCirclePlus /> New Canvas
						</Button>
					</BtnTooltip>
				</>
			)}
		</>
	);
}

const FileUpload = forwardRef((props, ref) => {
	const fileInputRef = useRef(null);

	const handleClick = () => {
		fileInputRef.current.click();
	};

	const handleFileChange = (event) => {
		const file = event.target.files[0];
		const formData = new FormData();
		formData.append("file", file);

		const basePath = import.meta.env.BASE_URL || "/";
		fetch(`${basePath}upload`, {
			method: "POST",
			body: formData,
		});

		if (ref) {
			// Check if ref is provided
			ref.current = fileInputRef.current; // Forward the ref to the underlying DOM element
		}
	};

	return (
		<div>
			<input
				type="file"
				ref={fileInputRef}
				onChange={handleFileChange}
				style={{ display: "none" }} // Hide the file input visually
			/>
			<BtnTooltip text="Upload">
				<Button
					variant="outline-primary"
					className="mx-1"
					onClick={handleClick}
				>
					<FaUpload />
				</Button>
			</BtnTooltip>
		</div>
	);
});

interface HeadBarProps {
	room: string;
	colorMode: string;
	handleColorMode: any;
	setIsDrawing: any;
	setGeometries: any;
	setPoints: any;
	isDrawing: boolean;
	tutorialURL: string;
	showSiMGen: boolean;
	modifierQueue: number;
	isAuthenticated: boolean;
	roomLock: boolean;
	setAddPlotsWindow: any;
	messages: any[];
	setMessages: any;
	token: string;
	step: number;
	selection: Set<number>;
}

const HeadBar = ({
	room,
	colorMode,
	handleColorMode,
	setIsDrawing,
	setGeometries,
	setPoints,
	isDrawing,
	tutorialURL,
	showSiMGen,
	modifierQueue,
	isAuthenticated,
	roomLock,
	setAddPlotsWindow,
	messages,
	setMessages,
	token,
	step,
	selection,
}: HeadBarProps) => {
	const [helpModalShow, setHelpModalShow] = useState(false);
	const [connectModalShow, setConnectModalShow] = useState(false);
	const [refreshModalShow, setRefreshModalShow] = useState(false);
	const [tutorialModalShow, setTutorialModalShow] = useState(false);
	const [consoleShow, setConsoleShow] = useState(false);
	const [remoteFileModalShow, setRemoteFileModalShow] = useState(false);
	const [zntrackAvailable, setZntrackAvailable] = useState(false);

	useEffect(() => {
		socket.emit("zntrack:available", (available: boolean) => {
			setZntrackAvailable(available);
		});
	}, []);

	useEffect(() => {
		setConsoleShow(showSiMGen);
	}, [showSiMGen]);

	const handleRemovePointsGeometries = () => {
		console.log("remove points and geometries");
		socket.emit("room:geometry:set", []);
		socket.emit("room:point:set", []);
		setGeometries([]);
		setPoints([]);
	};
	const basePath = useMemo(() => import.meta.env.BASE_URL || "/", []);

	return (
		<>
			<Navbar
				expand="md"
				className="bg-body-tertiary fixed-top"
				style={{ height: 50 }}
			>
				<Container fluid>
					<Navbar.Brand>
						<Button
							className="px-0 py-0 btn-lg"
							variant="tertiary"
							href="https://github.com/zincware/zndraw"
							target="_blank"
						>
							ZnDraw {version}
						</Button>
						{showSiMGen && (
							<>
								<Button className="px-0 py-0 btn-lg" variant="tertiary">
									+
								</Button>
								<Button
									className="px-0 py-0 btn-lg"
									variant="tertiary"
									href="https://github.com/RokasEl/simgen"
									target="_blank"
								>
									SiMGen
								</Button>
							</>
						)}
					</Navbar.Brand>
					<Navbar.Toggle aria-controls="basic-navbar-nav" />
					<Navbar.Collapse id="basic-navbar-nav">
						<Nav className="me-auto">
							<BtnTooltip text="Reset Scene">
								<Button
									variant="outline-danger"
									className="mx-1"
									onClick={() => {
										setRefreshModalShow((prev: boolean) => !prev);
									}}
								>
									<FaArrowRotateRight />
								</Button>
							</BtnTooltip>
							<BtnTooltip text="Activate Drawing Tool">
								<ToggleButton
									variant="outline-primary"
									className="mx-1"
									value="1"
									id="toggle-drawing"
									active={isDrawing}
									onClick={(e) => {
										setIsDrawing((prev: boolean) => !prev);
									}}
								>
									<FaPencil />
								</ToggleButton>
							</BtnTooltip>
							<BtnTooltip text="Remove all guiding points and geometries">
								<Button
									variant="outline-primary"
									className="mx-1"
									onClick={handleRemovePointsGeometries}
								>
									<FaHandSparkles />
								</Button>
							</BtnTooltip>
							<SiMGenButtons visible={showSiMGen} token={token} />
						</Nav>
						<Nav className="ms-auto">
							{tutorialURL && (
								<BtnTooltip text="Tutorial">
									<Button
										variant="warning"
										className="mx-1"
										onClick={() => setTutorialModalShow(true)}
									>
										Tutorial <FaFilm />
									</Button>
								</BtnTooltip>
							)}
							<BtnTooltip text="Access console">
								<ToggleButton
									variant="outline-primary"
									className="mx-1"
									value="1"
									id="toggle-console"
									active={consoleShow}
									onClick={() => {
										setConsoleShow((prev: boolean) => !prev);
									}}
								>
									<FaTerminal />
								</ToggleButton>
							</BtnTooltip>
							<BtnTooltip text="Python access">
								<Button
									variant="outline-primary"
									className="mx-1"
									onClick={() => setConnectModalShow(true)}
								>
									<FaCode />
								</Button>
							</BtnTooltip>
							<BtnTooltip text="Add plots window">
								<Button
									variant="outline-primary"
									className="mx-1"
									onClick={() => setAddPlotsWindow((prev: number) => prev + 1)}
								>
									<MdAddChart />
								</Button>
							</BtnTooltip>
							<FileUpload />
							<BtnTooltip text="Download">
								<Button
									variant="outline-primary"
									className="mx-1"
									href={`${basePath}download`}
									target="_blank"
								>
									<FaDownload />
								</Button>
							</BtnTooltip>
							{zntrackAvailable && (
								<BtnTooltip text="Open File via DVC">
									<Button
										variant="outline-primary"
										className="mx-1"
										onClick={() => {
											setRemoteFileModalShow(true);
										}}
									>
										<FaCloudDownloadAlt />
									</Button>
								</BtnTooltip>
							)}
							<BtnTooltip text="Help">
								<Button
									variant="outline-primary"
									className="mx-1"
									onClick={() => setHelpModalShow(true)}
								>
									<GrHelpBook />
								</Button>
							</BtnTooltip>
							<BtnTooltip text="Switch Colormode">
								<Button
									variant="outline-primary"
									className="mx-1"
									onClick={handleColorMode}
								>
									{colorMode === "light" ? <FaSun /> : <FaMoon />}
								</Button>
							</BtnTooltip>
							{isAuthenticated && (
								<>
									<BtnTooltip
										text={roomLock ? "Unlock this room" : "Lock this room"}
									>
										<Button
											variant="outline-danger"
											className="mx-1"
											onClick={() => {
												socket.emit("room:lock:set", !roomLock);
											}}
										>
											{roomLock ? <FaLock /> : <FaLockOpen />}
										</Button>
									</BtnTooltip>
									{/* <Button variant="outline-primary" className="mx-1">
                <FaUsers />
              </Button> */}
									<BtnTooltip text="Close ZnDraw">
										<Button
											variant="outline-danger"
											className="mx-1"
											onClick={() => {
												socket.emit("shutdown");
												close();
											}}
										>
											<MdExitToApp />
										</Button>
									</BtnTooltip>
								</>
							)}
						</Nav>
					</Navbar.Collapse>
				</Container>
			</Navbar>
			<HelpModel show={helpModalShow} onHide={() => setHelpModalShow(false)} />
			<ConnectModal
				show={connectModalShow}
				onHide={() => setConnectModalShow(false)}
				room={room}
			/>
			<RefreshModal
				show={refreshModalShow}
				onHide={() => setRefreshModalShow(false)}
				room={room}
			/>
			<TutorialModal
				show={tutorialModalShow}
				onHide={() => setTutorialModalShow(false)}
				url={tutorialURL}
			/>
			<RemoteFileModal
				show={remoteFileModalShow}
				onHide={() => setRemoteFileModalShow(false)}
				colorMode={colorMode}
			/>
			{consoleShow && (
				<ConsoleWindow
					messages={messages}
					setConsoleShow={setConsoleShow}
					token={token}
					setMessages={setMessages}
					colorMode={colorMode}
					step={step}
					selection={selection}
				/>
			)}
		</>
	);
};

export default HeadBar;
