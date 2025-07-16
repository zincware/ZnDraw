import { Close as CloseIcon, Edit as EditIcon, Save as SaveIcon } from "@mui/icons-material";
import {
    Box,
    Button,
    Card,
    CardActions,
    CardContent,
    CardHeader,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    FormControl,
    IconButton,
    InputLabel,
    MenuItem,
    Select,
    Stack,
    TextField,
    Typography,
} from "@mui/material";
import "katex/dist/katex.min.css";
import React, { memo, useCallback, useEffect, useRef, useState } from "react";
import Markdown from "react-markdown";
import { Rnd } from "react-rnd";
import { Prism as SyntaxHighlighter } from "react-syntax-highlighter";
import { oneDark, oneLight } from "react-syntax-highlighter/dist/esm/styles/prism";
import rehypeKatex from "rehype-katex";
import rehypeRaw from "rehype-raw";
import remarkBreaks from "remark-breaks";
import remarkGfm from "remark-gfm";
import remarkMath from "remark-math";

// --- Type Definitions ---

interface Message {
    msg: string;
    time: string;
}

interface ConsoleWindowProps {
    messages: Message[];
    setMessages: React.Dispatch<React.SetStateAction<Message[]>>;
    setConsoleShow: (show: boolean) => void;
    token: string;
    colorMode: "light" | "dark";
    step: number;
    selection: Set<number>;
}

// --- Helper & Child Components ---

const CodeBlock = memo(({ language, value, theme }: any) => (
    <SyntaxHighlighter PreTag="div" language={language} style={theme === "light" ? oneLight : oneDark}>
        {String(value).replace(/\n$/, "")}
    </SyntaxHighlighter>
));

const MessageItem = memo(({ message, index, onSave, colorMode }: { message: Message; index: number; onSave: (index: number, newMsg: string) => void; colorMode: "light" | "dark"; }) => {
    const [isEditing, setIsEditing] = useState(false);
    const [tempMsg, setTempMsg] = useState(message?.msg || "");

    useEffect(() => {
        setTempMsg(message?.msg || "");
    }, [message]);

    const handleSave = () => {
        onSave(index, tempMsg);
        setIsEditing(false);
    };

    const handleKeyDown = (e: React.KeyboardEvent) => {
        if (e.key === "Enter" && !e.shiftKey) {
            e.preventDefault();
            e.stopPropagation(); 
            handleSave();
        }
    };

    if (!message) return null;

    return (
        <Box mb={2} sx={{ textAlign: 'left' }}>
            <Stack direction="row" justifyContent="flex-start" alignItems="center">
                <Typography variant="caption" color="text.secondary">{message.time}</Typography>
                <IconButton size="small" sx={{ ml: 1 }} onClick={() => isEditing ? handleSave() : setIsEditing(true)}>
                    {isEditing ? <SaveIcon fontSize="inherit" /> : <EditIcon fontSize="inherit" />}
                </IconButton>
            </Stack>
            {isEditing ? (
                <TextField multiline fullWidth variant="outlined" value={tempMsg} onChange={(e) => setTempMsg(e.target.value)} onKeyDown={handleKeyDown} />
            ) : (
                <Markdown
                    remarkPlugins={[remarkMath, remarkGfm, remarkBreaks]}
                    rehypePlugins={[rehypeKatex, rehypeRaw]}
                    components={{
                        code: ({ node, className, children, ...props }) => {
                            const match = /language-(\w+)/.exec(className || "");
                            return match ? <CodeBlock language={match[1]} value={String(children)} theme={colorMode} /> : <code className={className} {...props}>{children}</code>;
                        },
                    }}
                >
                    {message.msg}
                </Markdown>
            )}
        </Box>
    );
});

// --- Main Component ---

export function ConsoleWindow({ messages = [], setMessages, setConsoleShow, token, colorMode, step, selection }: ConsoleWindowProps) {
    const [chatInput, setChatInput] = useState("");
    const [showInsertModal, setShowInsertModal] = useState(false);
    const scrollRef = useRef<HTMLDivElement>(null);

    useEffect(() => {
        if (scrollRef.current) {
            scrollRef.current.scrollTop = scrollRef.current.scrollHeight;
        }
    }, [messages]);
    
    const handleInsert = (textToInsert: string) => {
        setChatInput(prev => prev.slice(0, -2) + textToInsert);
    };

    const handleSendMessage = useCallback(() => {
        if (chatInput.trim() === "") return;
        const newMessage: Message = { msg: chatInput, time: new Date().toLocaleTimeString("de-DE") };
        setMessages(prev => [...prev, newMessage]);
        setChatInput("");
    }, [chatInput, setMessages]);

    const handleKeyPress = (e: React.KeyboardEvent) => {
        if (e.key === "Enter" && !e.shiftKey) {
            e.preventDefault();
            handleSendMessage();
        }
    };

    const handleSaveEdit = useCallback((index: number, newMsg: string) => {
        setMessages(prev =>
            prev.map((msg, idx) => (idx === index ? { ...msg, msg: newMsg } : msg))
        );
    }, [setMessages]);

    return (
        <>
            <Rnd
                default={{ x: window.innerWidth / 2 - 200, y: 100, width: 400, height: 500 }}
                minHeight={300} minWidth={300}
                bounds="window" style={{ zIndex: 1050 }}
                dragHandleClassName="drag-handle"
            >
                <Card sx={{ width: "100%", height: "100%", display: "flex", flexDirection: "column" }}>
                    <CardHeader
                        className="drag-handle"
                        title="Chat"
                        action={<IconButton onClick={() => setConsoleShow(false)}><CloseIcon /></IconButton>}
                        sx={{ cursor: "move", borderBottom: 1, borderColor: 'divider' }}
                    />
                    <CardContent sx={{ flexGrow: 1, overflowY: "auto", p: 2 }} ref={scrollRef}>
                        {messages.filter(Boolean).map((message, idx) => (
                            <MessageItem key={idx} index={idx} message={message} onSave={handleSaveEdit} colorMode={colorMode} />
                        ))}
                    </CardContent>
                    <CardActions sx={{ p: 2, borderTop: 1, borderColor: 'divider' }}>
                        <TextField
                            multiline
                            maxRows={4}
                            placeholder="Type a message..."
                            value={chatInput}
                            onChange={(e) => setChatInput(e.target.value)}
                            onKeyDown={handleKeyPress}
                            variant="outlined"
                            fullWidth
                            size="small"
                        />
                        <Button variant="contained" onClick={handleSendMessage} sx={{ ml: 1 }}>Send</Button>
                    </CardActions>
                </Card>
            </Rnd>
            <ChatInsertModal
                show={showInsertModal}
                onHide={() => setShowInsertModal(false)}
                onInsert={handleInsert}
                step={step}
                selection={selection}
            />
        </>
    );
}

// This component was unchanged but is included for completeness.
const ChatInsertModal = ({ show, onHide, onInsert, step, selection }: any) => {
    const handleSelectChange = (value: string) => {
        const basePath = `${window.location.origin}${window.location.pathname}`.replace(/\/+$/, "").replace(/\/token\/.*/, "");
        let textToInsert = "";

        if (value === "step") {
            textToInsert = `[step ${step}](${basePath}/?step=${step})`;
        } else if (value === "selection") {
            const selectionString = selection.size > 0 ? Array.from(selection).join(",") : "null";
            textToInsert = `[selection](${basePath}/?selection=${selectionString})`;
        }

        if (textToInsert) onInsert(textToInsert);
        onHide();
    };

    return (
        <Dialog open={show} onClose={onHide} fullWidth maxWidth="xs">
            <DialogTitle>Insert ZnDraw Link</DialogTitle>
            <DialogContent>
                <FormControl fullWidth sx={{ mt: 1 }}>
                    <InputLabel id="insert-select-label">Insert...</InputLabel>
                    <Select labelId="insert-select-label" value="" label="Insert..." onChange={(e) => handleSelectChange(e.target.value)}>
                        <MenuItem value="step">Current Step</MenuItem>
                        <MenuItem value="selection">Current Selection</MenuItem>
                    </Select>
                </FormControl>
            </DialogContent>
            <DialogActions>
                <Button onClick={onHide}>Close</Button>
            </DialogActions>
        </Dialog>
    );
};
