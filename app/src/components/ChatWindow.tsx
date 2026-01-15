import { useAppStore } from "../store";
import { remarkFrameLink } from "../utils/remark-frame-link.js";

import { useState, useRef, useEffect, memo, useCallback } from "react";
import {
	Box,
	TextField,
	IconButton,
	Typography,
	Paper,
	CircularProgress,
	LinearProgress,
	Divider,
	useMediaQuery,
} from "@mui/material";
import SendIcon from "@mui/icons-material/Send";
import CloseIcon from "@mui/icons-material/Close";
import EditIcon from "@mui/icons-material/Edit";
import CheckIcon from "@mui/icons-material/Check";
import CancelIcon from "@mui/icons-material/Cancel";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
import FullscreenIcon from "@mui/icons-material/Fullscreen";
import FullscreenExitIcon from "@mui/icons-material/FullscreenExit";
import ReactMarkdown from "react-markdown";
import { useParams } from "react-router-dom";
import {
	useChatMessages,
	useSendMessage,
	useEditMessage,
} from "../hooks/useChat";
import type { ChatMessage } from "../types/chat";
import { Rnd } from "react-rnd";
import TextareaAutosize from "react-textarea-autosize";
import { format } from "date-fns";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import { oneLight } from "react-syntax-highlighter/dist/esm/styles/prism";
import { oneDark } from "react-syntax-highlighter/dist/esm/styles/prism";
import { FrameReference } from "./FrameReference";
import { useColorScheme, useTheme } from "@mui/material/styles";
import { LAYOUT_CONSTANTS } from "../constants/layout";
import { MoleculePreview } from "./shared/MoleculePreview";

import "katex/dist/katex.min.css";

const CHAT_WIDTH = 400;
const CHAT_HEIGHT = 500;
const CHAT_PADDING = 20;

interface ChatWindowProps {
	open: boolean;
	onClose: () => void;
}

const markdownPlugins = [remarkMath, remarkFrameLink];
const htmlPlugins = [rehypeKatex];

type ProgressColor =
	| "primary"
	| "secondary"
	| "success"
	| "error"
	| "warning"
	| "info";

interface ProgressConfig {
	value?: number;
	min: number;
	max: number;
	description?: string;
	color: ProgressColor;
}

const parseProgressConfig = (content: string): ProgressConfig => {
	const config: ProgressConfig = {
		min: 0,
		max: 100,
		color: "primary",
	};

	const lines = content.trim().split("\n");
	for (const line of lines) {
		const colonIndex = line.indexOf(":");
		if (colonIndex === -1) continue;

		const key = line.substring(0, colonIndex).trim().toLowerCase();
		const rawValue = line.substring(colonIndex + 1).trim();

		switch (key) {
			case "value": {
				const numValue = parseFloat(rawValue);
				if (!Number.isNaN(numValue)) config.value = numValue;
				break;
			}
			case "min": {
				const minValue = parseFloat(rawValue);
				if (!Number.isNaN(minValue)) config.min = minValue;
				break;
			}
			case "max": {
				const maxValue = parseFloat(rawValue);
				if (!Number.isNaN(maxValue)) config.max = maxValue;
				break;
			}
			case "description":
				config.description = rawValue;
				break;
			case "color": {
				const validColors: ProgressColor[] = [
					"primary",
					"secondary",
					"success",
					"error",
					"warning",
					"info",
				];
				if (validColors.includes(rawValue as ProgressColor)) {
					config.color = rawValue as ProgressColor;
				}
				break;
			}
		}
	}

	return config;
};

const ProgressRenderer = ({ content }: { content: string }) => {
	const config = parseProgressConfig(content);

	const percentage =
		config.value !== undefined
			? ((config.value - config.min) / (config.max - config.min)) * 100
			: null;

	return (
		<Box sx={{ my: 1, display: "flex", flexDirection: "column", gap: 0.5 }}>
			{config.description && (
				<Typography variant="body2" color="text.secondary">
					{config.description}
				</Typography>
			)}
			{percentage !== null ? (
				<>
					<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
						<LinearProgress
							variant="determinate"
							value={Math.min(100, Math.max(0, percentage))}
							color={config.color}
							sx={{ height: 6, borderRadius: 3, flexGrow: 1 }}
						/>
						<Typography
							variant="caption"
							color="text.secondary"
							sx={{ minWidth: 40, textAlign: "right" }}
						>
							{Math.round(percentage)}%
						</Typography>
					</Box>
					<Typography variant="caption" color="text.secondary">
						{config.value} / {config.max}
						{config.min !== 0 && ` (min: ${config.min})`}
					</Typography>
				</>
			) : (
				<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
					<CircularProgress size={20} color={config.color} />
					<Typography variant="caption" color="text.secondary">
						In progress...
					</Typography>
				</Box>
			)}
		</Box>
	);
};

const SmilesRenderer = ({ content }: { content: string }) => {
	const smiles = content.trim();

	if (!smiles) {
		return null;
	}

	return (
		<Box sx={{ my: 1 }}>
			<MoleculePreview smiles={smiles} size="medium" throttleMs={0} />
		</Box>
	);
};

const MemoizedMarkdown = memo(ReactMarkdown);

/** Calculate initial position for the chat window */
const getInitialPosition = () => {
	const availableHeight =
		window.innerHeight -
		LAYOUT_CONSTANTS.APPBAR_HEIGHT -
		LAYOUT_CONSTANTS.PROGRESSBAR_HEIGHT;
	return {
		x: Math.max(CHAT_PADDING, window.innerWidth - CHAT_WIDTH - CHAT_PADDING),
		y: Math.max(
			LAYOUT_CONSTANTS.APPBAR_HEIGHT + CHAT_PADDING,
			LAYOUT_CONSTANTS.APPBAR_HEIGHT +
				availableHeight -
				CHAT_HEIGHT -
				CHAT_PADDING,
		),
	};
};

/** Clamp position to keep window within viewport bounds */
const clampPosition = (x: number, y: number, width: number, height: number) => {
	const maxX = Math.max(CHAT_PADDING, window.innerWidth - width - CHAT_PADDING);
	const minY = LAYOUT_CONSTANTS.APPBAR_HEIGHT + CHAT_PADDING;
	const maxY = Math.max(
		minY,
		window.innerHeight -
			LAYOUT_CONSTANTS.PROGRESSBAR_HEIGHT -
			height -
			CHAT_PADDING,
	);

	return {
		x: Math.min(Math.max(CHAT_PADDING, x), maxX),
		y: Math.min(Math.max(minY, y), maxY),
	};
};

const ChatWindow = ({ open, onClose }: ChatWindowProps) => {
	// Use individual selectors to prevent unnecessary re-renders
	const setCurrentFrame = useAppStore((state) => state.setCurrentFrame);
	const resetChatUnread = useAppStore((state) => state.resetChatUnread);
	const userName = useAppStore((state) => state.userName);
	const { roomId } = useParams<{ roomId: string }>();
	const [messageInput, setMessageInput] = useState("");
	const [editingMessageId, setEditingMessageId] = useState<string | null>(null);
	const [editContent, setEditContent] = useState("");
	const messagesEndRef = useRef<HTMLDivElement>(null);
	const scrollContainerRef = useRef<HTMLDivElement>(null);
	const { mode } = useColorScheme();
	const theme = useTheme();

	// Fullscreen and mobile handling
	const isMobile = useMediaQuery(theme.breakpoints.down("sm"));
	const [isFullscreen, setIsFullscreen] = useState(false);
	const [position, setPosition] = useState(getInitialPosition);
	const [size, setSize] = useState({ width: CHAT_WIDTH, height: CHAT_HEIGHT });
	const rndRef = useRef<Rnd>(null);

	// Auto-fullscreen on mobile
	useEffect(() => {
		setIsFullscreen(isMobile);
	}, [isMobile]);

	// Handle viewport resize - keep chat within bounds
	useEffect(() => {
		const handleResize = () => {
			if (isFullscreen) return; // Don't adjust position in fullscreen mode

			const clamped = clampPosition(
				position.x,
				position.y,
				size.width,
				size.height,
			);
			if (clamped.x !== position.x || clamped.y !== position.y) {
				setPosition(clamped);
				// Update Rnd component position
				rndRef.current?.updatePosition(clamped);
			}
		};

		window.addEventListener("resize", handleResize);
		return () => window.removeEventListener("resize", handleResize);
	}, [position, size, isFullscreen]);

	// Reset position when opening chat (ensures proper initial position)
	useEffect(() => {
		if (open && !isFullscreen) {
			const initialPos = getInitialPosition();
			setPosition(initialPos);
			rndRef.current?.updatePosition(initialPos);
		}
	}, [open, isFullscreen]);

	const { data, fetchNextPage, hasNextPage, isFetchingNextPage, isLoading } =
		useChatMessages(roomId || "", 50);
	const sendMessage = useSendMessage(roomId || "");
	const editMessage = useEditMessage(roomId || "");

	useEffect(() => {
		if (open) {
			document.body.style.overflow = "hidden";
		} else {
			document.body.style.overflow = "unset";
		}
		return () => {
			document.body.style.overflow = "unset";
		};
	}, [open]);

	const scrollToBottom = () => {
		if (scrollContainerRef.current) {
			scrollContainerRef.current.scrollTop =
				scrollContainerRef.current.scrollHeight;
		}
	};

	useEffect(() => {
		if (open) {
			// Use timeout to ensure DOM is updated before scrolling
			setTimeout(scrollToBottom, 100);
			// Reset unread count when chat is opened
			resetChatUnread();
		}
	}, [open, data, resetChatUnread]); // Rerun when opened or when new data arrives

	const handleScroll = (e: React.UIEvent<HTMLDivElement>) => {
		if (scrollContainerRef.current) {
			const { scrollTop } = scrollContainerRef.current;
			if (scrollTop < 50 && hasNextPage && !isFetchingNextPage) {
				fetchNextPage();
			}
		}
	};

	const handleWheel = (e: React.WheelEvent<HTMLDivElement>) => {
		if (scrollContainerRef.current) {
			const { scrollTop, scrollHeight, clientHeight } =
				scrollContainerRef.current;
			// If scrolling up at the top or scrolling down at the bottom, prevent bubbling
			if (
				(e.deltaY < 0 && scrollTop === 0) ||
				(e.deltaY > 0 && scrollTop + clientHeight >= scrollHeight - 1)
			) {
				// Allow scroll if there's more content to load upwards
				if (e.deltaY < 0 && hasNextPage) {
					return;
				}
				e.preventDefault();
				e.stopPropagation();
			}
		}
	};

	const handleSendMessage = async () => {
		if (!messageInput.trim()) return;
		try {
			await sendMessage.mutateAsync(messageInput);
			setMessageInput("");
		} catch (error) {
			console.error("Failed to send message:", error);
		}
	};

	const handleStartEdit = (message: ChatMessage) => {
		setEditingMessageId(message.id);
		setEditContent(message.content);
	};

	const handleSaveEdit = async () => {
		if (!editingMessageId || !editContent.trim()) return;
		try {
			await editMessage.mutateAsync({
				messageId: editingMessageId,
				content: editContent,
			});
			setEditingMessageId(null);
			setEditContent("");
		} catch (error) {
			console.error("Failed to edit message:", error);
		}
	};

	const handleCancelEdit = () => {
		setEditingMessageId(null);
		setEditContent("");
	};

	const allMessages =
		data?.pages
			.flatMap((page) => page.messages)
			.sort((a, b) => a.createdAt - b.createdAt) || [];

	if (!open) {
		return null;
	}

	// Fullscreen mode - render as fixed overlay
	if (isFullscreen) {
		return (
			<Paper
				elevation={4}
				sx={{
					position: "fixed",
					top: LAYOUT_CONSTANTS.APPBAR_HEIGHT,
					left: 0,
					right: 0,
					bottom: LAYOUT_CONSTANTS.PROGRESSBAR_HEIGHT,
					zIndex: 1301,
					display: "flex",
					flexDirection: "column",
					overflow: "hidden",
					borderRadius: 0,
				}}
			>
				{/* Header */}
				<Box
					sx={{
						p: 1.4,
						display: "flex",
						alignItems: "center",
						justifyContent: "space-between",
						borderBottom: 1,
						borderColor: "divider",
					}}
				>
					<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
						<Typography variant="h6">Chat</Typography>
						<IconButton size="small" sx={{ color: "text.secondary" }}>
							<HelpOutlineIcon fontSize="small" />
						</IconButton>
					</Box>
					<Box>
						{!isMobile && (
							<IconButton
								onClick={() => setIsFullscreen(false)}
								size="small"
								title="Exit fullscreen"
							>
								<FullscreenExitIcon />
							</IconButton>
						)}
						<IconButton onClick={onClose} size="small">
							<CloseIcon />
						</IconButton>
					</Box>
				</Box>

				{/* Messages Container */}
				<Box
					ref={scrollContainerRef}
					onScroll={handleScroll}
					onWheel={handleWheel}
					sx={{
						flexGrow: 1,
						overflowY: "auto",
						p: 2,
						display: "flex",
						flexDirection: "column",
						gap: 1,
					}}
				>
					{isLoading && (
						<Box sx={{ display: "flex", justifyContent: "center", p: 2 }}>
							<CircularProgress size={24} />
						</Box>
					)}
					{isFetchingNextPage && (
						<Box sx={{ display: "flex", justifyContent: "center", p: 1 }}>
							<CircularProgress size={20} />
						</Box>
					)}

					{allMessages.map((message) => {
						const isOwnMessage = message.author.id === userName;
						const isEditing = editingMessageId === message.id;

						return (
							<Paper
								key={message.id}
								sx={{
									p: 1.5,
									backgroundColor: isOwnMessage
										? mode === "dark"
											? "primary.dark"
											: "primary.light"
										: "background.paper",
									color: isOwnMessage
										? mode === "dark"
											? "primary.contrastText"
											: "text.primary"
										: "text.primary",
									alignSelf: isOwnMessage ? "flex-end" : "flex-start",
									maxWidth: "80%",
								}}
								elevation={1}
							>
								<Box
									sx={{
										display: "flex",
										justifyContent: "space-between",
										alignItems: "center",
										mb: 0.5,
									}}
								>
									<Typography variant="caption" color="text.secondary">
										{message.author.id} -{" "}
										{format(new Date(message.createdAt), "HH:mm")}
									</Typography>
									{isOwnMessage && !isEditing && (
										<IconButton
											size="small"
											onClick={() => handleStartEdit(message)}
										>
											<EditIcon fontSize="small" />
										</IconButton>
									)}
								</Box>

								{isEditing ? (
									<Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
										<TextField
											fullWidth
											size="small"
											value={editContent}
											onChange={(e) => setEditContent(e.target.value)}
											multiline
											maxRows={4}
										/>
										<IconButton
											size="small"
											onClick={handleSaveEdit}
											color="primary"
										>
											<CheckIcon fontSize="small" />
										</IconButton>
										<IconButton size="small" onClick={handleCancelEdit}>
											<CancelIcon fontSize="small" />
										</IconButton>
									</Box>
								) : (
									<>
										<Box sx={{ "& p": { margin: 0 }, wordBreak: "break-word" }}>
											<MemoizedMarkdown
												remarkPlugins={markdownPlugins}
												rehypePlugins={htmlPlugins}
												components={
													{
														code({ node, className, children, ...props }) {
															const match = /language-(\w+)/.exec(
																className || "",
															);
															const language = match?.[1];
															const content = String(children).replace(
																/\n$/,
																"",
															);

															if (language === "progress") {
																return <ProgressRenderer content={content} />;
															}

															if (language === "smiles") {
																return <SmilesRenderer content={content} />;
															}

															return match ? (
																<SyntaxHighlighter
																	style={
																		mode === "dark"
																			? (oneDark as {
																					[key: string]: React.CSSProperties;
																				})
																			: (oneLight as {
																					[key: string]: React.CSSProperties;
																				})
																	}
																	language={language}
																	PreTag="div"
																	{...props}
																>
																	{content}
																</SyntaxHighlighter>
															) : (
																<code className={className} {...props}>
																	{children}
																</code>
															);
														},

														frameLink: ({ frame }: { frame: number }) => {
															if (typeof frame === "number") {
																return <FrameReference frame={frame} />;
															}
															return null;
														},
													} as Record<string, React.ComponentType<any>>
												}
											>
												{message.content}
											</MemoizedMarkdown>
										</Box>
										{message.isEdited && (
											<Typography
												variant="caption"
												color="text.secondary"
												sx={{ fontStyle: "italic", mt: 0.5 }}
											>
												(edited)
											</Typography>
										)}
									</>
								)}
							</Paper>
						);
					})}
					<div ref={messagesEndRef} />
				</Box>

				<Divider />

				{/* Input Area */}
				<Box sx={{ p: 2, backgroundColor: "background.default" }}>
					<Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
						<TextareaAutosize
							minRows={1}
							maxRows={6}
							placeholder="Type a message..."
							value={messageInput}
							onChange={(e) => setMessageInput(e.target.value)}
							onKeyPress={(e) => {
								if (e.key === "Enter" && !e.shiftKey) {
									e.preventDefault();
									handleSendMessage();
								}
							}}
							style={{
								width: "100%",
								resize: "none",
								padding: "8.5px 14px",
								borderRadius: "4px",
								border: `1px solid ${theme.palette.divider}`,
								backgroundColor: theme.palette.background.paper,
								color: theme.palette.text.primary,
								fontFamily: "inherit",
								fontSize: "inherit",
							}}
						/>
						<IconButton
							color="primary"
							onClick={handleSendMessage}
							disabled={!messageInput.trim() || sendMessage.isPending}
						>
							<SendIcon />
						</IconButton>
					</Box>
				</Box>
			</Paper>
		);
	}

	// Normal draggable/resizable mode
	return (
		<Rnd
			ref={rndRef}
			position={position}
			size={size}
			onDragStop={(e, d) => {
				const clamped = clampPosition(d.x, d.y, size.width, size.height);
				setPosition(clamped);
			}}
			onResizeStop={(e, direction, ref, delta, pos) => {
				const newWidth = parseInt(ref.style.width, 10);
				const newHeight = parseInt(ref.style.height, 10);
				const clamped = clampPosition(pos.x, pos.y, newWidth, newHeight);
				setSize({ width: newWidth, height: newHeight });
				setPosition(clamped);
			}}
			minWidth={300}
			minHeight={400}
			bounds="window"
			dragHandleClassName="drag-handle"
			style={{ zIndex: 1301 }}
		>
			<Paper
				elevation={4}
				sx={{
					width: "100%",
					height: "100%",
					display: "flex",
					flexDirection: "column",
					overflow: "hidden",
				}}
			>
				{/* Header */}
				<Box
					className="drag-handle"
					sx={{
						p: 1.4,
						display: "flex",
						alignItems: "center",
						justifyContent: "space-between",
						borderBottom: 1,
						borderColor: "divider",
						cursor: "move",
					}}
				>
					<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
						<Typography variant="h6">Chat</Typography>
						<IconButton size="small" sx={{ color: "text.secondary" }}>
							<HelpOutlineIcon fontSize="small" />
						</IconButton>
					</Box>
					<Box>
						<IconButton
							onClick={() => setIsFullscreen(true)}
							size="small"
							title="Fullscreen"
						>
							<FullscreenIcon />
						</IconButton>
						<IconButton onClick={onClose} size="small">
							<CloseIcon />
						</IconButton>
					</Box>
				</Box>

				{/* Messages Container */}
				<Box
					ref={scrollContainerRef}
					onScroll={handleScroll}
					onWheel={handleWheel}
					sx={{
						flexGrow: 1,
						overflowY: "auto",
						p: 2,
						display: "flex",
						flexDirection: "column",
						gap: 1,
					}}
				>
					{isLoading && (
						<Box sx={{ display: "flex", justifyContent: "center", p: 2 }}>
							<CircularProgress size={24} />
						</Box>
					)}
					{isFetchingNextPage && (
						<Box sx={{ display: "flex", justifyContent: "center", p: 1 }}>
							<CircularProgress size={20} />
						</Box>
					)}

					{allMessages.map((message) => {
						const isOwnMessage = message.author.id === userName;
						const isEditing = editingMessageId === message.id;

						return (
							<Paper
								key={message.id}
								sx={{
									p: 1.5,
									backgroundColor: isOwnMessage
										? mode === "dark"
											? "primary.dark"
											: "primary.light"
										: "background.paper",
									color: isOwnMessage
										? mode === "dark"
											? "primary.contrastText"
											: "text.primary"
										: "text.primary",
									alignSelf: isOwnMessage ? "flex-end" : "flex-start",
									maxWidth: "80%",
								}}
								elevation={1}
							>
								<Box
									sx={{
										display: "flex",
										justifyContent: "space-between",
										alignItems: "center",
										mb: 0.5,
									}}
								>
									<Typography variant="caption" color="text.secondary">
										{message.author.id} -{" "}
										{format(new Date(message.createdAt), "HH:mm")}
									</Typography>
									{isOwnMessage && !isEditing && (
										<IconButton
											size="small"
											onClick={() => handleStartEdit(message)}
										>
											<EditIcon fontSize="small" />
										</IconButton>
									)}
								</Box>

								{isEditing ? (
									<Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
										<TextField
											fullWidth
											size="small"
											value={editContent}
											onChange={(e) => setEditContent(e.target.value)}
											multiline
											maxRows={4}
										/>
										<IconButton
											size="small"
											onClick={handleSaveEdit}
											color="primary"
										>
											<CheckIcon fontSize="small" />
										</IconButton>
										<IconButton size="small" onClick={handleCancelEdit}>
											<CancelIcon fontSize="small" />
										</IconButton>
									</Box>
								) : (
									<>
										<Box sx={{ "& p": { margin: 0 }, wordBreak: "break-word" }}>
											<MemoizedMarkdown
												remarkPlugins={markdownPlugins}
												rehypePlugins={htmlPlugins}
												components={
													{
														code({ node, className, children, ...props }) {
															const match = /language-(\w+)/.exec(
																className || "",
															);
															const language = match?.[1];
															const content = String(children).replace(
																/\n$/,
																"",
															);

															if (language === "progress") {
																return <ProgressRenderer content={content} />;
															}

															if (language === "smiles") {
																return <SmilesRenderer content={content} />;
															}

															return match ? (
																<SyntaxHighlighter
																	style={
																		mode === "dark"
																			? (oneDark as {
																					[key: string]: React.CSSProperties;
																				})
																			: (oneLight as {
																					[key: string]: React.CSSProperties;
																				})
																	}
																	language={language}
																	PreTag="div"
																	{...props}
																>
																	{content}
																</SyntaxHighlighter>
															) : (
																<code className={className} {...props}>
																	{children}
																</code>
															);
														},

														frameLink: ({ frame }: { frame: number }) => {
															if (typeof frame === "number") {
																return <FrameReference frame={frame} />;
															}
															return null;
														},
													} as Record<string, React.ComponentType<any>>
												}
											>
												{message.content}
											</MemoizedMarkdown>
										</Box>
										{message.isEdited && (
											<Typography
												variant="caption"
												color="text.secondary"
												sx={{ fontStyle: "italic", mt: 0.5 }}
											>
												(edited)
											</Typography>
										)}
									</>
								)}
							</Paper>
						);
					})}
					<div ref={messagesEndRef} />
				</Box>

				<Divider />

				{/* Input Area */}
				<Box sx={{ p: 2, backgroundColor: "background.default" }}>
					<Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
						<TextareaAutosize
							minRows={1}
							maxRows={6}
							placeholder="Type a message..."
							value={messageInput}
							onChange={(e) => setMessageInput(e.target.value)}
							onKeyPress={(e) => {
								if (e.key === "Enter" && !e.shiftKey) {
									e.preventDefault();
									handleSendMessage();
								}
							}}
							style={{
								width: "100%",
								resize: "none",
								padding: "8.5px 14px",
								borderRadius: "4px",
								border: `1px solid ${theme.palette.divider}`,
								backgroundColor: theme.palette.background.paper,
								color: theme.palette.text.primary,
								fontFamily: "inherit",
								fontSize: "inherit",
							}}
						/>
						<IconButton
							color="primary"
							onClick={handleSendMessage}
							disabled={!messageInput.trim() || sendMessage.isPending}
						>
							<SendIcon />
						</IconButton>
					</Box>
				</Box>
			</Paper>
		</Rnd>
	);
};

// Memoize to prevent unnecessary re-renders
export default memo(ChatWindow);
