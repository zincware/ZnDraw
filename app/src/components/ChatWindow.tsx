import { useAppStore } from "../store";
import { remarkFrameLink } from "../utils/remark-frame-link.js";

import { useState, useRef, useEffect, memo } from "react";
import {
  Box,
  TextField,
  IconButton,
  Typography,
  Paper,
  CircularProgress,
  Divider,
} from "@mui/material";
import SendIcon from "@mui/icons-material/Send";
import CloseIcon from "@mui/icons-material/Close";
import EditIcon from "@mui/icons-material/Edit";
import CheckIcon from "@mui/icons-material/Check";
import CancelIcon from "@mui/icons-material/Cancel";
import HelpOutlineIcon from "@mui/icons-material/HelpOutline";
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

import "katex/dist/katex.min.css";

interface ChatWindowProps {
  open: boolean;
  onClose: () => void;
}

const markdownPlugins = [remarkMath, remarkFrameLink];
const htmlPlugins = [rehypeKatex];

const MemoizedMarkdown = memo(ReactMarkdown);

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

  return (
    <Rnd
      default={{
        x: window.innerWidth - 420,
        y: 20,
        width: 400,
        height: 600,
      }}
      minWidth={300}
      minHeight={400}
      bounds=".drag-boundary-container"
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
          <IconButton onClick={onClose} size="small">
            <CloseIcon />
          </IconButton>
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
                        components={{
                          code({ node, className, children, ...props }) {
                            const match = /language-(\w+)/.exec(
                              className || "",
                            );
                            return match ? (
                              <SyntaxHighlighter
                                style={mode === "dark" ? oneDark : oneLight}
                                language={match[1]}
                                PreTag="div"
                                {...props}
                              >
                                {String(children).replace(/\n$/, "")}
                              </SyntaxHighlighter>
                            ) : (
                              <code className={className} {...props}>
                                {children}
                              </code>
                            );
                          },

                          // 👇 map frameLink directly
                          frameLink({ frame }) {
                            // Destructure 'frame' directly from props
                            if (typeof frame === "number") {
                              return <FrameReference frame={frame} />;
                            }
                            // Optional: return a fallback if frame is missing for some reason
                            return null;
                          },

                          // fallback for real links
                          // a({ node, children, ...props }) {
                          //   return (
                          //     <a {...props} style={{ color: "#1976d2" }}>
                          //       {children}
                          //     </a>
                          //   );
                          // },
                        }}
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
