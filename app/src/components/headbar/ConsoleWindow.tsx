import { Close as CloseIcon } from "@mui/icons-material";
import {
	Box,
	Button,
	Card,
	CardActions,
	CardContent,
	CardHeader,
	Container,
	Dialog,
	DialogActions,
	DialogContent,
	DialogTitle,
	FormControl,
	FormControlLabel,
	IconButton,
	InputLabel,
	MenuItem,
	Select,
	Switch,
	TextField,
} from "@mui/material";
import type React from "react";
import { useEffect, useRef, useState } from "react";
import rehypeKatex from "rehype-katex";
import rehypeRaw from "rehype-raw";
import remarkBreaks from "remark-breaks";
import remarkGfm from "remark-gfm";
import remarkMath from "remark-math";
import "katex/dist/katex.min.css"; // `rehype-katex` does not import the CSS for you
import { FaSave } from "react-icons/fa";
import Markdown from "react-markdown";
import { Rnd } from "react-rnd";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import {
	oneDark,
	oneLight,
} from "react-syntax-highlighter/dist/esm/styles/prism";

import { FaPencil } from "react-icons/fa6";

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
		<Dialog
			open={show}
			onClose={onHide}
			aria-labelledby="contained-modal-title-vcenter"
			maxWidth="lg"
			fullWidth
		>
			<DialogTitle id="contained-modal-title-vcenter">
				ZnDraw Chat Insert
			</DialogTitle>
			<DialogContent>
				<Container>
					Share your current view:
					<FormControl fullWidth sx={{ mt: 2 }}>
						<InputLabel>Choose...</InputLabel>
						<Select
							value=""
							onChange={(e) => handleSelectChange({ value: e.target.value })}
							label="Choose..."
						>
							{options.map((option) => (
								<MenuItem key={option.value} value={option.value}>
									{option.label}
								</MenuItem>
							))}
						</Select>
					</FormControl>
				</Container>
			</DialogContent>
			<DialogActions>
				<Button onClick={onHide}>Close</Button>
			</DialogActions>
		</Dialog>
	);
}

export function ConsoleWindow({
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
					<CardHeader
						title="Chat"
						action={
							<Box sx={{ display: "flex", alignItems: "center" }}>
								<FormControlLabel
									control={
										<Switch
											checked={showTime}
											onChange={() => setShowTime(!showTime)}
										/>
									}
									label="Show Time"
									sx={{ mr: 2 }}
								/>
								<IconButton onClick={() => setConsoleShow(false)}>
									<CloseIcon />
								</IconButton>
							</Box>
						}
					/>

					{/* Message Body with Optional Timestamp */}
					<CardContent
						sx={{ textAlign: "start", overflowY: "auto" }}
						ref={scrollRef}
					>
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
										<TextField
											multiline
											rows={4}
											value={tempMsg}
											onChange={(e) => setTempMsg(e.target.value)}
											onKeyDown={(e) => handleEditKeyPress(e, idx)}
											variant="outlined"
											fullWidth
										/>
									)}
								</div>
							</div>
						))}
					</CardContent>

					<CardActions sx={{ display: "flex", gap: 1 }}>
						<TextField
							multiline
							rows={1}
							placeholder="Type a message..."
							onInput={handleChatInputChange}
							onKeyDown={handleKeyPress}
							ref={chatInputRef}
							variant="outlined"
							fullWidth
							size="small"
						/>
						<Button variant="contained" onClick={handleSendMessage}>
							Send
						</Button>
					</CardActions>
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
