# Chat Feature Implementation Plan

## Overview
Implement a real-time chat system for ZnDraw rooms with markdown support, message editing, and infinite scroll.

## Architecture Decisions

### 1. Message ID Strategy
- **Format**: `msg_<room_id>_<counter>` (e.g., `msg_room1_42`)
- **Rationale**: Smaller than UUIDs, readable, sequential
- **Implementation**: Use Redis INCR `room:<room_id>:chat:counter`

### 3. Data Modeling in Redis

#### Message Counter
```
Key: room:<room_id>:chat:counter
Type: String (integer)
Purpose: Atomic counter for message IDs
```

#### Message Index (Sorted Set)
```
Key: room:<room_id>:chat:index
Type: Sorted Set
Members: message IDs (e.g., "msg_room1_42")
Scores: Unix timestamp in milliseconds
Purpose: Fast range queries for pagination
```

#### Message Data (Hash)
```
Key: room:<room_id>:chat:data
Type: Hash
Fields: message IDs
Values: JSON strings
Purpose: Fast lookups and updates
```

#### Message Object Schema
```json
{
  "id": "msg_room1_42",
  "roomId": "room1",
  "authorSid": "abc123",
  "authorUserId": "user1",
  "content": "Hello **world**!",
  "createdAt": 1727906231000,
  "updatedAt": 1727906231000,
  "isEdited": false
}
```

### 4. Pagination Strategy
- **Default limit**: 30 messages
- **Max limit**: 100 messages
- **Direction**: Descending by default (newest first)
- **Parameters**:
  - `limit`: Number of messages (1-100)
  - `before`: Timestamp - get messages before this time
  - `after`: Timestamp - get messages after this time
- **Response Metadata**:
  - `hasMore`: Boolean indicating more messages available
  - `totalCount`: Total messages in room
  - `oldestTimestamp`: Timestamp of oldest message returned
  - `newestTimestamp`: Timestamp of newest message returned

## Implementation Tasks

### Phase 1: Backend Foundation

#### Task 1.1: Socket Event Handlers (src/zndraw/app/events.py)
```python
@socketio.on("chat:message:create")
def handle_chat_message_create(data):
    """
    Payload: { "content": "message text" }
    Returns: { "success": bool, "message": Message | None, "error": str | None }
    """
    # 3. Generate message ID using INCR
    # 4. Create message object with metadata
    # 5. Store in Redis (HSET + ZADD)
    # 6. Emit 'chat:message:new' to room
    # 7. Return success response

@socketio.on("chat:message:edit")
def handle_chat_message_edit(data):
    """
    Payload: { "messageId": "msg_room_42", "content": "new text" }
    Returns: { "success": bool, "message": Message | None, "error": str | None }
    """
    # 2. Fetch message from Redis
    # 5. Update message (content, updatedAt, isEdited=true)
    # 6. Store in Redis (HSET)
    # 7. Emit 'chat:message:updated' to room
    # 8. Return success response
```

**Implementation Details**:
- Use `get_project_room_from_session(sid)` for room validation
- Error responses: Include descriptive error messages
- Emit to room using `socketio.emit(..., room=room_id)`

#### Task 1.2: REST API Endpoint (src/zndraw/app/routes.py)
```python
@main.route("/api/rooms/<string:room_id>/chat/messages", methods=["GET"])
def get_chat_messages(room_id: str):
    """
    Query Parameters:
        - limit (int): Number of messages (default: 30, max: 100)
        - before (int): Get messages before this timestamp
        - after (int): Get messages after this timestamp

    Returns:
    {
        "messages": [Message],
        "metadata": {
            "hasMore": bool,
            "totalCount": int,
            "oldestTimestamp": int | null,
            "newestTimestamp": int | null
        }
    }
    """
    # 1. Validate room exists
    # 2. Parse and validate query parameters
    # 3. Use ZREVRANGEBYSCORE for range query
    # 4. Fetch message data with HMGET
    # 5. Calculate metadata
    # 6. Return JSON response
```

**Implementation Details**:
- Use `ZCARD` for total count
- Use `ZREVRANGEBYSCORE` for descending order
- Use `ZRANGEBYSCORE` for ascending order (with "after")
- Efficient bulk fetch with pipeline: `r.pipeline()` + `HMGET`

#### Task 1.3: Helper Functions (src/zndraw/app/chat_utils.py - NEW FILE)
```python
def create_message(redis_client, room_id: str, author_sid: str,
                   author_user_id: str, content: str) -> dict:
    """Create and store a new message."""

def get_message(redis_client, room_id: str, message_id: str) -> dict | None:
    """Fetch a message by ID."""

def update_message(redis_client, room_id: str, message_id: str,
                   content: str) -> dict:
    """Update message content."""
```

### Phase 2: Python Client Integration

#### Task 2.1: Client Methods (src/zndraw/client.py)
```python
class ZnDrawClient:
    def log(self, message: str) -> dict:
        """
        Send a chat message to the room.

        Args:
            message: Chat message content (markdown supported)

        Returns:
            dict: Response with success status and message data
        """
        response = self.sio.call("chat:message:create", {"content": message})
        return response

    def edit_message(self, message_id: str, new_content: str) -> dict:
        """Edit a previously sent message."""
        response = self.sio.call("chat:message:edit", {
            "messageId": message_id,
            "content": new_content
        })
        return response

    def get_messages(self, limit: int = 30, before: int = None) -> list[dict]:
        """Fetch chat message history."""
        params = {"limit": limit}
        if before:
            params["before"] = before
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/chat/messages",
            params=params
        )
        response.raise_for_status()
        return response.json()
```

### Phase 3: Frontend Implementation

#### Task 3.1: TypeScript Types (app/src/types/chat.ts - NEW FILE)
```typescript
export interface ChatMessage {
  id: string;
  roomId: string;
  authorSid: string;
  authorUserId: string;
  content: string;
  createdAt: number;
  updatedAt: number;
  isEdited: boolean;
}

export interface ChatMessagesResponse {
  messages: ChatMessage[];
  metadata: {
    hasMore: boolean;
    totalCount: number;
    oldestTimestamp: number | null;
    newestTimestamp: number | null;
  };
}
```

#### Task 3.2: React Query Hook (app/src/hooks/useChat.ts - NEW FILE)
```typescript
import { useQuery, useMutation, useQueryClient, useInfiniteQuery } from '@tanstack/react-query';
import { socket } from '../socket';

// Fetch chat messages with infinite scroll support
export const useChatMessages = (room: string, limit: number = 30) => {
  return useInfiniteQuery({
    queryKey: ['chat', room],
    queryFn: async ({ pageParam }) => {
      const params = new URLSearchParams();
      params.append('limit', limit.toString());
      if (pageParam) {
        params.append('before', pageParam.toString());
      }
      const response = await fetch(`/api/rooms/${room}/chat/messages?${params}`);
      if (!response.ok) throw new Error('Failed to fetch messages');
      return response.json();
    },
    getNextPageParam: (lastPage) =>
      lastPage.metadata.hasMore ? lastPage.metadata.oldestTimestamp : undefined,
    enabled: !!room,
    staleTime: 1000 * 60, // 1 minute
  });
};

// Send new message
export const useSendMessage = (room: string) => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: async (content: string) => {
      return new Promise((resolve, reject) => {
        socket.emit('chat:message:create', { content }, (response) => {
          if (response.success) resolve(response.message);
          else reject(new Error(response.error));
        });
      });
    },
    onSuccess: () => {
      // Messages are added via socket event, so we just invalidate
      queryClient.invalidateQueries({ queryKey: ['chat', room] });
    },
  });
};

// Edit message
export const useEditMessage = (room: string) => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: async ({ messageId, content }: { messageId: string; content: string }) => {
      return new Promise((resolve, reject) => {
        socket.emit('chat:message:edit', { messageId, content }, (response) => {
          if (response.success) resolve(response.message);
          else reject(new Error(response.error));
        });
      });
    },
    onSuccess: () => {
      // Messages are updated via socket event
      queryClient.invalidateQueries({ queryKey: ['chat', room] });
    },
  });
};
```

#### Task 3.3: Socket Integration (app/src/hooks/useSocketManager.ts - UPDATE)
```typescript
// Add to existing useSocketManager hook
function onChatMessageNew(data: ChatMessage) {
  queryClient.setQueryData(['chat', room], (oldData: any) => {
    if (!oldData) return oldData;
    // Add new message to the end of the first page
    const newPages = [...oldData.pages];
    newPages[0] = {
      ...newPages[0],
      messages: [...newPages[0].messages, data],
      metadata: {
        ...newPages[0].metadata,
        totalCount: newPages[0].metadata.totalCount + 1,
      }
    };
    return { ...oldData, pages: newPages };
  });
}

function onChatMessageUpdated(data: ChatMessage) {
  queryClient.setQueryData(['chat', room], (oldData: any) => {
    if (!oldData) return oldData;
    // Find and update the message
    const newPages = oldData.pages.map(page => ({
      ...page,
      messages: page.messages.map(msg =>
        msg.id === data.id ? data : msg
      ),
    }));
    return { ...oldData, pages: newPages };
  });
}

socket.on('chat:message:new', onChatMessageNew);
socket.on('chat:message:updated', onChatMessageUpdated);
```

#### Task 3.4: Chat Component (app/src/components/ChatWindow.tsx - NEW FILE)
- Use the previous implementation as reference
- Key improvements:
  - Remove rehypeRaw, use sanitization
  - Add infinite scroll with react-query
  - Optimistic updates for sending messages
  - Show edit indicator on edited messages
  - Better error handling and loading states
  - Proper TypeScript types
  - Accessibility improvements

### Phase 4: Testing

#### Task 4.1: Backend Tests (tests/test_chat.py - NEW FILE)
```python
def test_create_chat_message(server):
    """Test creating a new chat message via socket"""

def test_edit_chat_message(server):
    """Test editing an existing message"""

def test_get_chat_messages_rest(server):
    """Test REST API for fetching messages"""

def test_chat_pagination(server):
    """Test pagination with before/after parameters"""

def test_chat_room_isolation(server):
    """Test messages are isolated per room"""
```

#### Task 4.2: Integration Tests
```python
def test_chat_python_client(server):
    """Test vis.log() method"""
```

## Implementation Order

1. **Backend Core** (Day 1)
   - [ ] Create chat_utils.py with helper functions
   - [ ] Implement socket event handlers
   - [ ] Implement REST API endpoint
   - [ ] Add basic tests

2. **Python Client** (Day 1-2)
   - [ ] Add log() method to ZnDrawClient
   - [ ] Add edit_message() method
   - [ ] Add get_messages() method
   - [ ] Test with pytest

3. **Frontend Foundation** (Day 2)
   - [ ] Create TypeScript types
   - [ ] Create useChat hook
   - [ ] Update useSocketManager
   - [ ] Test query/mutation logic

4. **UI Components** (Day 3)
   - [ ] Create ChatWindow component
   - [ ] Implement infinite scroll
   - [ ] Add markdown rendering
   - [ ] Add edit functionality
   - [ ] Polish UI/UX

5. **Testing & Polish** (Day 3-4)
   - [ ] Complete test coverage
   - [ ] Handle edge cases
   - [ ] Performance optimization

## Performance Considerations

1. **Redis Optimization**
   - Use pipelining for bulk operations
   - Consider message TTL or max count (e.g., 10k messages/room)
   - Index cleanup for old messages

2. **Frontend Optimization**
   - Virtual scrolling for large message lists (react-window)
   - Debounce edit operations
   - Optimistic updates for better UX
   - Lazy load markdown components

3. **Network Optimization**
   - Cache query results appropriately

## Future Enhancements (Not in Scope)

- Message reactions/emoji
- File attachments
- Message threading
- Search functionality
- Mentions/notifications
- Read receipts
- Message deletion
- User presence indicators
- Typing indicators

-----------------------------------------------

## REVIEW
This is an exceptionally well-structured and detailed implementation plan. It demonstrates a strong understanding of full-stack development, from data modeling to frontend state management. It's comprehensive, clear, and provides a solid roadmap for development.

My review will focus on validating the excellent choices made and offering a few minor refinements to further enhance the robustness and efficiency of the proposed system.

-----

### Overall Assessment

**Excellent.** This is a professional-grade plan that is clear, actionable, and technically sound. The architecture is robust, the task breakdown is logical, and the considerations for performance and testing are thorough. You are in a great position to start implementation.

-----

### Architectural Review

#### **Confirmations (Excellent Choices)**

  * **Redis Data Modeling:** The use of a `counter` for IDs, a `Sorted Set` for the index, and a `Hash` for the data is the ideal Redis pattern for this use case. It provides atomicity, fast ordered queries, and efficient individual message lookups/updates.
  * **API Separation:** The division of labor is perfect: REST for fetching paginated historical data and WebSockets for real-time creation and updates.
  * **Pagination Strategy:** Cursor-based pagination using timestamps (`before`/`after`) is the correct approach for a real-time feed. The metadata (`hasMore`, `totalCount`, etc.) is well-designed and provides the frontend with everything it needs.
  * **Full-Stack Cohesion:** The plan clearly connects the dots from the backend data schema to the Python client methods to the TypeScript types and React Query hooks. This is a sign of a well-thought-out system.

#### **Suggestions for Refinement**

1.  **Message Object Schema (`authorSid` vs. `authorUserId`):**

      * **Observation:** The schema includes both `authorSid` and `authorUserId`. A Socket.IO `sid` is ephemeral; it changes every time a user reconnects (e.g., refreshes the page). The `authorUserId` (from the `Flask-Login` session) is the persistent identifier for that user's session.
      * **Suggestion:** **Remove `authorSid` from the schema.** It provides little value once the message is stored and can be misleading. The `authorUserId` is the single source of truth for who sent the message. You could enhance the `author` object to be more structured, which is more extensible for the future.
      * **Revised Schema:**
        ```json
        {
          "id": "msg_room1_42",
          "author": {
            "id": "user1" // The stable session user ID
            // "name": "User 1" // You could add a display name in the future
          },
          "content": "Hello **world**!",
          "createdAt": 1727906231000,
          "updatedAt": 1727906231000,
          "isEdited": false
        }
        ```

2.  **REST API Response Metadata (`totalCount`):**

      * **Observation:** Including `totalCount` requires a `ZCARD` call on every request. For rooms with hundreds of thousands of messages, this is an O(1) operation and very fast, but it's still an extra command.
      * **Suggestion:** This is a minor point, but consider if `totalCount` is needed on every paginated response. Often, just `hasMore` is sufficient for the UI. If you do want to show a total count, you could fetch it once on initial load and then update it client-side as new messages arrive. Your current plan is perfectly fine, just something to be aware of for hyper-optimization.

-----

### Implementation Review

#### **Backend**

  * **Task 1.1 (`chat:message:edit`):** The implementation steps are missing a critical check.
      * **Suggestion:** Add a step: **"4. Verify `current_user.user_id` matches the message's `authorUserId`."** This is a crucial authorization check to prevent one user from editing another user's messages.

#### **Python Client**

  * **Task 2.1:** The use of `sio.call` is excellent. It simplifies the client-side logic by turning an async event with a callback into a synchronous-style request/response, which is perfect for a library.

#### **Frontend**

  * **Task 3.2 (React Query Hooks):** This is very well-designed, but there's a common "gotcha" that can be refined for better performance.

      * **Observation:** In `useSendMessage` and `useEditMessage`, the `onSuccess` handler calls `queryClient.invalidateQueries`. This will trigger a full refetch of the latest page of messages via the REST API.
      * **Suggestion:** **Remove `invalidateQueries` from the mutations.** Your `useSocketManager` is already set up to listen for the `chat:message:new` and `chat:message:updated` broadcast events. These listeners will update the React Query cache in real-time using `setQueryData`. This is more efficient than invalidating and refetching. The mutation's job is to *send* the command; the socket listener's job is to handle the resulting state update. This decouples the actions cleanly.

  * **Task 3.3 (Socket Integration):** The cache update logic is almost perfect, but the page order needs to be considered carefully.

      * **Observation:** The `onChatMessageNew` function adds the new message to `newPages[0].messages`. In `useInfiniteQuery`, the pages are typically stored in the order they were fetched. When fetching backwards in time (`before`), `page[0]` is the *oldest* page of messages currently in memory. New messages should appear at the very end of the chat.
      * **Suggestion:** New messages should be added to the **last page** in the cache, as that page contains the most recent data.
      * **Revised Logic:**
        ```typescript
        function onChatMessageNew(data: ChatMessage) {
          queryClient.setQueryData(['chat', room], (oldData: any) => {
            if (!oldData) return oldData;
            const newPages = [...oldData.pages];
            const lastPageIndex = newPages.length - 1; // Target the last page
            newPages[lastPageIndex] = {
              ...newPages[lastPageIndex],
              messages: [...newPages[lastPageIndex].messages, data],
              // You might not need to update metadata here if you're not using it reactively
            };
            return { ...oldData, pages: newPages };
          });
        }
        ```

### Final Verdict

This is a stellar plan. The suggested refinements are minor and intended to harden an already excellent architecture. You have a clear path forward. Proceed with confidence.

## Phase 5: Frontend Polish and Fixes

This phase addresses several UI/UX issues in the chat frontend to make it fully functional and user-friendly.

### Task 5.1: Make Chat Window Draggable and Resizable
- **Problem**: The chat window is fixed and hidden behind other UI elements.
- **Solution**: Integrate `react-rnd` to allow users to move and resize the chat window.
- **Steps**:
  1. Install `react-rnd` and its types: `bun add react-rnd @types/react-rnd`.
  2. In `ChatWindow.tsx`, wrap the main component container with the `<Rnd>` component.
  3. Configure a default size and position for the window to ensure it's visible on load.
  4. Style the window to indicate it's draggable (e.g., with a header bar).

### Task 5.2: Fix Message Submission and Input Resizing
- **Problem**: Messages cannot be sent from the frontend, and the input area is not suitable for long messages.
- **Solution**: Debug the submission logic and replace the standard `textarea` with an auto-sizing one.
- **Steps**:
  1. Install `react-textarea-autosize`: `bun add react-textarea-autosize`.
  2. In `ChatWindow.tsx`, replace the `<textarea>` with the `TextareaAutosize` component.
  3. Investigate the `onSubmit` handler and the `useSendMessage` mutation hook to identify and fix the submission issue. Ensure `socket.emit` is called correctly.
  4. Ensure pressing "Enter" submits the message and "Shift+Enter" creates a new line.

### Task 5.3: Correct Message Order and Add Timestamps
- **Problem**: Messages are displayed newest-first, but a chat interface should have the newest messages at the bottom. Timestamps are missing.
- **Solution**: Reverse the message order for display and add formatted timestamps. Implement auto-scrolling.
- **Steps**:
  1. In `ChatWindow.tsx`, process the data from `useInfiniteQuery`. Flatten the `pages` array and then reverse the resulting `messages` array to display oldest messages first.
  2. For each message, display the `createdAt` timestamp. Format it readably (e.g., "HH:mm"). Consider installing `date-fns` for robust date formatting.
  3. Implement an auto-scroll mechanism. Use a `ref` on the message container and a `useEffect` hook that scrolls to the bottom whenever the messages array changes.

### Task 5.4: Enable Markdown, Math, and Syntax Highlighting
- **Problem**: Markdown, LaTeX-style math, and code syntax highlighting are not being rendered.
- **Solution**: Install and correctly configure the necessary libraries for markdown rendering.
- **Steps**:
  1. Install dependencies: `bun add react-markdown remark-math rehype-katex react-syntax-highlighter katex`. Also add `@types/react-syntax-highlighter`.
  2. In `ChatWindow.tsx`, use the `ReactMarkdown` component to render message content.
  3. Pass the required plugins to `ReactMarkdown`: `remarkPlugins={[remarkMath]}`, `rehypePlugins={[rehypeKatex]}`.
  4. Provide a custom component for code blocks (`<code>`) that uses `react-syntax-highlighter` to render highlighted code.
  5. Ensure the Katex CSS is imported into the project for math formulas to render correctly. This can be done in `index.html` or `index.css`.