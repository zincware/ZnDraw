import CancelIcon from "@mui/icons-material/Cancel";
import CheckIcon from "@mui/icons-material/Check";
import EditIcon from "@mui/icons-material/Edit";
import SendIcon from "@mui/icons-material/Send";
import {
	Box,
	CircularProgress,
	Divider,
	IconButton,
	LinearProgress,
	Paper,
	TextField,
	Typography,
} from "@mui/material";
import { useColorScheme, useTheme } from "@mui/material/styles";
import { format } from "date-fns";
import { memo, useCallback, useEffect, useRef, useState } from "react";
import ReactMarkdown from "react-markdown";
import { useParams } from "react-router-dom";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import {
	oneDark,
	oneLight,
} from "react-syntax-highlighter/dist/esm/styles/prism";
import TextareaAutosize from "react-textarea-autosize";
import rehypeKatex from "rehype-katex";
import remarkMath from "remark-math";
import { FrameReference } from "../components/FrameReference";
import { MoleculePreview } from "../components/shared/MoleculePreview";
import {
	useChatMessages,
	useEditMessage,
	useSendMessage,
} from "../hooks/useChat";
import { socket } from "../socket";
import { useAppStore } from "../store";
import type { ChatMessage } from "../types/chat";
import { remarkFrameLink } from "../utils/remark-frame-link.js";

import "katex/dist/katex.min.css";

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
				const numValue = Number.parseFloat(rawValue);
				if (!Number.isNaN(numValue)) config.value = numValue;
				break;
			}
			case "min": {
				const minValue = Number.parseFloat(rawValue);
				if (!Number.isNaN(minValue)) config.min = minValue;
				break;
			}
			case "max": {
				const maxValue = Number.parseFloat(rawValue);
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

const ChatPanel = () => {
	// Use individual selectors to prevent unnecessary re-renders
	const userName = useAppStore((state) => state.user?.email ?? null);
	const typingUsers = useAppStore((state) => state.typingUsers);
	const { roomId } = useParams<{ roomId: string }>();
	const [messageInput, setMessageInput] = useState("");
	const [editingMessageId, setEditingMessageId] = useState<number | null>(null);
	const [editContent, setEditContent] = useState("");
	const messagesEndRef = useRef<HTMLDivElement>(null);
	const scrollContainerRef = useRef<HTMLDivElement>(null);
	const typingTimeoutRef = useRef<ReturnType<typeof setTimeout> | null>(null);
	const { mode } = useColorScheme();
	const theme = useTheme();

	const { data, fetchNextPage, hasNextPage, isFetchingNextPage, isLoading } =
		useChatMessages(roomId || "", 50);
	const sendMessage = useSendMessage(roomId || "");
	const editMessage = useEditMessage(roomId || "");

	const scrollToBottom = useCallback(() => {
		if (scrollContainerRef.current) {
			scrollContainerRef.current.scrollTop =
				scrollContainerRef.current.scrollHeight;
		}
	}, []);

	// Scroll on new data, regardless of visibility.
	// biome-ignore lint/correctness/useExhaustiveDependencies: data is an intentional tick trigger — the effect scrolls when new messages arrive
	useEffect(() => {
		setTimeout(scrollToBottom, 100);
	}, [data, scrollToBottom]);

	// Reset unread count only while the panel is visible. ChatPanel renders
	// inside SidebarZone, so visibility is "the chat icon is the active tool
	// in some sidebar". Re-run on data changes and on activation so the badge
	// clears when the user switches back.
	const isChatActive = useAppStore(
		(s) =>
			s.activeLeft === "chat" ||
			s.activeRight === "chat" ||
			s.activeBottom === "chat",
	);
	// biome-ignore lint/correctness/useExhaustiveDependencies: data is an intentional tick trigger — the effect runs on new messages to clear the badge
	useEffect(() => {
		if (isChatActive) {
			useAppStore.getState().resetChatUnread();
		}
	}, [data, isChatActive]);

	const handleScroll = (_e: React.UIEvent<HTMLDivElement>) => {
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

	const emitTypingStart = useCallback(() => {
		if (!roomId) return;
		socket.emit("typing_start", { room_id: roomId });
		// Clear previous stop timeout
		if (typingTimeoutRef.current) clearTimeout(typingTimeoutRef.current);
		// Auto-stop after 3s of inactivity
		typingTimeoutRef.current = setTimeout(() => {
			socket.emit("typing_stop", { room_id: roomId });
			typingTimeoutRef.current = null;
		}, 3000);
	}, [roomId]);

	const emitTypingStop = useCallback(() => {
		if (!roomId) return;
		if (typingTimeoutRef.current) {
			clearTimeout(typingTimeoutRef.current);
			typingTimeoutRef.current = null;
		}
		socket.emit("typing_stop", { room_id: roomId });
	}, [roomId]);

	// Clean up typing timeout on unmount
	useEffect(() => {
		return () => {
			if (typingTimeoutRef.current) clearTimeout(typingTimeoutRef.current);
		};
	}, []);

	const handleSendMessage = async () => {
		if (!messageInput.trim()) return;
		emitTypingStop();
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
			.flatMap((page) => page.items)
			.sort((a, b) => a.created_at.localeCompare(b.created_at)) || [];

	const typingNames = Array.from(typingUsers);
	const typingText =
		typingNames.length === 1
			? `${typingNames[0]} is typing...`
			: typingNames.length > 1
				? `${typingNames.join(", ")} are typing...`
				: null;

	const markdownComponents = {
		code({ node: _node, className, children, ...props }: any) {
			const match = /language-(\w+)/.exec(className || "");
			const language = match?.[1];
			const content = String(children).replace(/\n$/, "");

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
							? (oneDark as { [key: string]: React.CSSProperties })
							: (oneLight as { [key: string]: React.CSSProperties })
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
	} as Record<string, React.ComponentType<any>>;

	return (
		<Box
			data-testid="chat-panel"
			sx={{
				display: "flex",
				flexDirection: "column",
				height: "100%",
				overflow: "hidden",
			}}
		>
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
					const isOwnMessage = message.email === userName;
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
									{message.email} -{" "}
									{format(new Date(message.created_at), "HH:mm")}
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
											components={markdownComponents}
										>
											{message.content}
										</MemoizedMarkdown>
									</Box>
									{message.updated_at !== null && (
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

			{/* Typing Indicator */}
			{typingText && (
				<Box sx={{ px: 2, py: 0.5 }}>
					<Typography
						variant="caption"
						color="text.secondary"
						sx={{ fontStyle: "italic" }}
					>
						{typingText}
					</Typography>
				</Box>
			)}

			{/* Input Area */}
			<Box sx={{ p: 2, backgroundColor: "background.default" }}>
				<Box sx={{ display: "flex", gap: 1, alignItems: "center" }}>
					<TextareaAutosize
						data-testid="chat-input"
						minRows={1}
						maxRows={6}
						placeholder="Type a message..."
						value={messageInput}
						onChange={(e) => {
							setMessageInput(e.target.value);
							if (e.target.value.trim()) emitTypingStart();
						}}
						onKeyDown={(e) => {
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
						data-testid="chat-send"
						color="primary"
						onClick={handleSendMessage}
						disabled={!messageInput.trim() || sendMessage.isPending}
					>
						<SendIcon />
					</IconButton>
				</Box>
			</Box>
		</Box>
	);
};

// Memoize to prevent unnecessary re-renders
export default memo(ChatPanel);
