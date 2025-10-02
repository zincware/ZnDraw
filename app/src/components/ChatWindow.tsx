import { useState, useRef, useEffect } from 'react';
import {
  Drawer,
  Box,
  TextField,
  IconButton,
  Typography,
  Paper,
  CircularProgress,
  Divider,
  Chip,
} from '@mui/material';
import SendIcon from '@mui/icons-material/Send';
import CloseIcon from '@mui/icons-material/Close';
import EditIcon from '@mui/icons-material/Edit';
import CheckIcon from '@mui/icons-material/Check';
import CancelIcon from '@mui/icons-material/Cancel';
import ReactMarkdown from 'react-markdown';
import { useParams } from 'react-router-dom';
import { useChatMessages, useSendMessage, useEditMessage } from '../hooks/useChat';
import type { ChatMessage } from '../types/chat';

interface ChatWindowProps {
  open: boolean;
  onClose: () => void;
}

const ChatWindow = ({ open, onClose }: ChatWindowProps) => {
  const { roomId, userId } = useParams<{ roomId: string; userId: string }>();
  const [messageInput, setMessageInput] = useState('');
  const [editingMessageId, setEditingMessageId] = useState<string | null>(null);
  const [editContent, setEditContent] = useState('');
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const scrollContainerRef = useRef<HTMLDivElement>(null);

  const { data, fetchNextPage, hasNextPage, isFetchingNextPage, isLoading } = useChatMessages(
    roomId || '',
    30
  );
  const sendMessage = useSendMessage(roomId || '');
  const editMessage = useEditMessage(roomId || '');

  // Auto-scroll to bottom on new messages
  useEffect(() => {
    if (messagesEndRef.current && !isFetchingNextPage) {
      messagesEndRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  }, [data, isFetchingNextPage]);

  // Handle infinite scroll
  const handleScroll = () => {
    if (!scrollContainerRef.current) return;
    const { scrollTop } = scrollContainerRef.current;

    // Load more when scrolled near the top
    if (scrollTop < 100 && hasNextPage && !isFetchingNextPage) {
      fetchNextPage();
    }
  };

  const handleSendMessage = async () => {
    if (!messageInput.trim()) return;

    try {
      await sendMessage.mutateAsync(messageInput);
      setMessageInput('');
    } catch (error) {
      console.error('Failed to send message:', error);
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
      setEditContent('');
    } catch (error) {
      console.error('Failed to edit message:', error);
    }
  };

  const handleCancelEdit = () => {
    setEditingMessageId(null);
    setEditContent('');
  };

  // Flatten all messages from all pages
  const allMessages = data?.pages.flatMap((page) => page.messages) || [];

  return (
    <Drawer
      anchor="right"
      open={open}
      onClose={onClose}
      sx={{
        '& .MuiDrawer-paper': {
          width: 400,
          maxWidth: '100%',
        },
      }}
    >
      <Box sx={{ display: 'flex', flexDirection: 'column', height: '100%' }}>
        {/* Header */}
        <Box
          sx={{
            p: 2,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            borderBottom: 1,
            borderColor: 'divider',
          }}
        >
          <Typography variant="h6">Chat</Typography>
          <IconButton onClick={onClose} size="small">
            <CloseIcon />
          </IconButton>
        </Box>

        {/* Messages Container */}
        <Box
          ref={scrollContainerRef}
          onScroll={handleScroll}
          sx={{
            flexGrow: 1,
            overflowY: 'auto',
            p: 2,
            display: 'flex',
            flexDirection: 'column',
            gap: 1,
          }}
        >
          {/* Loading indicator for initial load */}
          {isLoading && (
            <Box sx={{ display: 'flex', justifyContent: 'center', p: 2 }}>
              <CircularProgress size={24} />
            </Box>
          )}

          {/* Load more indicator */}
          {isFetchingNextPage && (
            <Box sx={{ display: 'flex', justifyContent: 'center', p: 1 }}>
              <CircularProgress size={20} />
            </Box>
          )}

          {/* Messages */}
          {allMessages.map((message) => {
            const isOwnMessage = message.author.id === userId;
            const isEditing = editingMessageId === message.id;

            return (
              <Paper
                key={message.id}
                sx={{
                  p: 1.5,
                  backgroundColor: isOwnMessage ? 'primary.light' : 'background.paper',
                  alignSelf: isOwnMessage ? 'flex-end' : 'flex-start',
                  maxWidth: '80%',
                }}
                elevation={1}
              >
                <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 0.5 }}>
                  <Typography variant="caption" color="text.secondary">
                    {message.author.id}
                  </Typography>
                  {isOwnMessage && !isEditing && (
                    <IconButton size="small" onClick={() => handleStartEdit(message)}>
                      <EditIcon fontSize="small" />
                    </IconButton>
                  )}
                </Box>

                {isEditing ? (
                  <Box sx={{ display: 'flex', gap: 1, alignItems: 'center' }}>
                    <TextField
                      fullWidth
                      size="small"
                      value={editContent}
                      onChange={(e) => setEditContent(e.target.value)}
                      multiline
                      maxRows={4}
                    />
                    <IconButton size="small" onClick={handleSaveEdit} color="primary">
                      <CheckIcon fontSize="small" />
                    </IconButton>
                    <IconButton size="small" onClick={handleCancelEdit}>
                      <CancelIcon fontSize="small" />
                    </IconButton>
                  </Box>
                ) : (
                  <>
                    <Box sx={{ '& p': { margin: 0 } }}>
                      <ReactMarkdown>{message.content}</ReactMarkdown>
                    </Box>
                    {message.isEdited && (
                      <Typography variant="caption" color="text.secondary" sx={{ fontStyle: 'italic', mt: 0.5 }}>
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
        <Box sx={{ p: 2 }}>
          <Box sx={{ display: 'flex', gap: 1 }}>
            <TextField
              fullWidth
              size="small"
              placeholder="Type a message... (Markdown supported)"
              value={messageInput}
              onChange={(e) => setMessageInput(e.target.value)}
              onKeyPress={(e) => {
                if (e.key === 'Enter' && !e.shiftKey) {
                  e.preventDefault();
                  handleSendMessage();
                }
              }}
              multiline
              maxRows={4}
            />
            <IconButton
              color="primary"
              onClick={handleSendMessage}
              disabled={!messageInput.trim() || sendMessage.isPending}
            >
              <SendIcon />
            </IconButton>
          </Box>
          <Typography variant="caption" color="text.secondary" sx={{ mt: 1, display: 'block' }}>
            Press Enter to send, Shift+Enter for new line
          </Typography>
        </Box>
      </Box>
    </Drawer>
  );
};

export default ChatWindow;
