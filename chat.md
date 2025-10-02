# Data Modeling in Redis

Sorted Set (for order): `room:<room_id>:chat:index`
This set will store message IDs.
The score of each member will be its creation timestamp (e.g., a Unix timestamp with milliseconds). This keeps all message IDs sorted by time automatically.


Hash (for data): `room:<room_id>:chat:data`

This hash will store the actual message objects.
The field will be the unique message ID.

Consider using
- react-syntax-highlighter
- react-markdown
- remark-gfm
- remark-math
- remark-breaks
- rehype-katex
- rehype-raw

The value will be a JSON string representing the message object. This allows for fast lookups and edits of a specific message without scanning the whole list.

```json
{
  "id": "msg_uuid_12345",
  "roomId": "myroom_1",
  "author": "sid",
  "content": "Hello world! This is **markdown**.",
  "createdAt": 1727906231000, // Unix timestamp in ms
  "updatedAt": 1727906231000 // Same as createdAt initially
}
```

# Socket events
- `chat:message:create`, payload: `{ string: content }`
Create a new message object (generate UUID, add author info, set timestamps).
Save the message to Redis (HSET for data, ZADD for the index).
Broadcast the full message object to the room via `chat:message:new`.
- `chat:message:edit`, payload: `{ string: messageId, string: newContent }`
Update the content and updatedAt fields.
Broadcast the full, updated message object via `chat:message:updated`.

- `chat:message:new`
Payload: The full new message object.
Clients listening will add this message to the bottom of their chat list.

- `chat:message:updated`
Payload: The full updated message object.
Clients listening will find the message by ID and update its content and timestamp.

# REST Endpoint (for History & Infinite Scroll)

GET /api/rooms/<string:room_id>/chat/messages

## Query Parameters:
limit (optional, default: 30): How many messages to fetch.
before (optional): The createdAt timestamp of the last message the client has. If not provided, it fetches the most recent page.


# TypeScript 
- use `react-query` for fetching history
- create a new hook
- use exiting zustand appstore if appropriate
- use socketmanger

previous implementation you can use as reference, but adapt the code for good practices and the new design.
```ts
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
```

# Python

Implement the chat via the `vis = ZnDraw` object as follows:

```python
vis = ZnDraw(...)
vis.log("Hello world!") # chat message
```

# Tests
Test creating new chat messages and updating them.
Access the chat history via the REST endpoint.
Create a new file tests/test_chat.py.
Look at the oder tests/*.py for examples.
Use pytest and keep each test function small and focused on one aspect of the chat functionality.

