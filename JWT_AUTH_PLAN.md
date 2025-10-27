# JWT Authentication Implementation Plan

## Executive Summary

Refactor authentication from **client-generated UUIDs** to **server-issued JWT tokens** for proper security.

**Key Changes:**
- Add `/api/login` endpoint for JWT generation
- Validate JWT in REST API endpoints
- Validate JWT in WebSocket `connect` handler
- Reject unauthenticated connections
- Keep client IDs for tracking (but server-issued)

---

## 1. Current State Analysis

### Security Issues
```python
# app/src/utils/clientId.ts:22
clientId = crypto.randomUUID()  # ❌ Client-generated = insecure
```

**Problems:**
- Any client can generate arbitrary UUIDs
- No server-side session control
- No authentication, only identification
- Client IDs can be stolen/reused

### Current Flow
1. Browser generates client ID (UUID)
2. Client sends X-Client-ID header to `/api/rooms/{room}/join`
3. Server accepts any UUID without validation
4. WebSocket connects with client ID in auth
5. Server trusts the client ID

---

## 2. Target Architecture

### JWT-Based Authentication

**JWT Claims:**
```json
{
  "sub": "client-uuid",        // Subject: unique client ID (server-generated)
  "userName": "John",          // User's display name
  "iat": 1234567890,           // Issued at timestamp
  "exp": 1234654290,           // Expiration timestamp (7 days)
  "jti": "token-uuid"          // JWT ID for revocation support (optional)
}
```

**Token Expiration Strategy:**
- Tokens do not expire (stateless sessions)
- Future: Add optional expiration via config if needed

### New Flow

```
┌─────────┐                    ┌────────┐
│ Browser │                    │ Server │
└────┬────┘                    └───┬────┘
     │                             │
     │ POST /api/login             │
     │ {"userName": "John"}        │
     │────────────────────────────>│
     │                             │ Generate clientId (UUID)
     │                             │ Create JWT token
     │                             │ Store client in Redis
     │                             │
     │ {"token": "jwt...",         │
     │  "clientId": "uuid"}        │
     │<────────────────────────────│
     │                             │
     │ Store token in localStorage │
     │                             │
     │ POST /api/rooms/X/join      │
     │ Authorization: Bearer jwt   │
     │────────────────────────────>│
     │                             │ Validate JWT
     │                             │ Extract clientId from claims
     │                             │
     │ {"status": "ok", ...}       │
     │<────────────────────────────│
     │                             │
     │ WebSocket connect           │
     │ auth: {token: "jwt..."}     │
     │=============================>│
     │                             │ Validate JWT
     │                             │ Extract clientId from claims
     │                             │ Associate SID with clientId
     │                             │
     │ {"status": "ok"}            │
     │<=============================│
```

---

## 3. Implementation Details

### 3.1 Secret Key Management

**File:** `src/zndraw/server.py:113`

**Current:**
```python
app.config["SECRET_KEY"] = "your_secret_key"  # ❌ Hardcoded
```

**New:**
```python
app.config["SECRET_KEY"] = os.getenv(
    "FLASK_SECRET_KEY",
    "dev-secret-key-change-in-production"  # Default for development
)
```

---


### 3.3 JWT Utility Functions

**File:** `src/zndraw/auth.py` (NEW)

```python
"""JWT authentication utilities."""
import logging
import time
import typing as t
import uuid
from functools import wraps

import jwt
from flask import current_app, request

log = logging.getLogger(__name__)


class AuthError(Exception):
    """Authentication error."""
    def __init__(self, message: str, status_code: int = 401):
        self.message = message
        self.status_code = status_code
        super().__init__(self.message)


def create_jwt_token(client_id: str, user_name: str | None) -> str:
    """Create JWT token for authenticated client.

    Parameters
    ----------
    client_id : str
        Unique client identifier (UUID)
    user_name : str | None
        Optional display name of the user

    Returns
    -------
    str
        Encoded JWT token
    """
    secret_key = current_app.config["SECRET_KEY"]
    algorithm = current_app.config.get("JWT_ALGORITHM", "HS256")


    payload = {
        "sub": client_id,          # Subject: client ID
        "userName": user_name,      # Display name
        "jti": str(uuid.uuid4()),   # JWT ID for revocation (future use)
    }

    token = jwt.encode(payload, secret_key, algorithm=algorithm)
    log.info(f"Created JWT for client {client_id}")

    return token


def decode_jwt_token(token: str) -> dict:
    """Decode and validate JWT token.

    Parameters
    ----------
    token : str
        JWT token string

    Returns
    -------
    dict
        Decoded payload with claims

    Raises
    ------
    AuthError
        If token is invalid, expired, or malformed
    """
    secret_key = current_app.config["SECRET_KEY"]
    algorithm = current_app.config.get("JWT_ALGORITHM", "HS256")

    try:
        payload = jwt.decode(token, secret_key, algorithms=[algorithm])
        return payload
    except jwt.InvalidTokenError as e:
        raise AuthError(f"Invalid token: {str(e)}", 401)


def extract_token_from_request() -> str | None:
    """Extract JWT token from request Authorization header.

    Expects: Authorization: Bearer <token>

    Returns
    -------
    str | None
        Token string if found, None otherwise
    """
    auth_header = request.headers.get("Authorization")
    if auth_header and auth_header.startswith("Bearer "):
        return auth_header[7:]  # Remove "Bearer " prefix
    return None


def get_current_client() -> dict:
    """Get current authenticated client from request.

    Returns
    -------
    dict
        JWT payload with clientId and userName

    Raises
    ------
    AuthError
        If no token found or token is invalid
    """
    token = extract_token_from_request()
    if not token:
        raise AuthError("No authentication token provided", 401)

    payload = decode_jwt_token(token)
    return {
        "clientId": payload["sub"],
        "userName": payload["userName"],
    }


def require_auth(f):
    """Decorator to require JWT authentication for route.

    Usage:
        @app.route("/api/protected")
        @require_auth
        def protected_route():
            client = get_current_client()
            return {"clientId": client["clientId"]}
    """
    @wraps(f)
    def decorated_function(*args, **kwargs):
        try:
            get_current_client()  # Validate token
            return f(*args, **kwargs)
        except AuthError as e:
            return {"error": e.message}, e.status_code
    return decorated_function
```

---

### 3.4 Login Endpoint

**File:** `src/zndraw/app/routes.py` (ADD NEW ENDPOINT)

```python
@main.route("/api/login", methods=["POST"])
def login():
    """Authenticate user and issue JWT token.

    Request:
        {
            "userName": "John Doe"  // Required
        }

    Response:
        {
            "status": "ok",
            "token": "eyJhbGc...",     // JWT token
            "clientId": "uuid-string",  // Server-generated client ID
            "expiresIn": 604800         // Seconds until expiration
        }
    """
    from zndraw.auth import create_jwt_token

    data = request.get_json() or {}
    user_name = data.get("userName")

    if not user_name or not user_name.strip():
        return {"error": "userName is required"}, 400

    # Generate server-side client ID
    import uuid
    client_id = str(uuid.uuid4())

    # Create JWT token
    token = create_jwt_token(client_id, user_name)

    # Store client metadata in Redis
    r = current_app.extensions["redis"]
    client_key = f"client:{client_id}"
    r.hset(client_key, mapping={
        "userName": user_name,
        "createdAt": time.time(),
        "lastLogin": time.time(),
    })

    log.info(f"User '{user_name}' logged in with client ID: {client_id}")

    return {
        "status": "ok",
        "token": token,
        "clientId": client_id,
    }
```

---

### 3.5 Update Join Endpoint

**File:** `src/zndraw/app/routes.py` (MODIFY EXISTING)

```python
@main.route("/api/rooms/<string:room_id>/join", methods=["POST"])
def join_room(room_id):
    """Join a room (requires JWT authentication).

    Headers:
        Authorization: Bearer <jwt-token> (required)

    Request:
        {
            "description": "optional room description",
            "copyFrom": "optional-source-room",
            "allowCreate": true
        }

    Response:
        {
            "status": "ok",
            "roomId": "room-name",
            "frameCount": 0,
            "step": 0,
            "created": true
        }
    """
    from zndraw.auth import get_current_client, AuthError

    # Authenticate request (JWT required)
    try:
        client = get_current_client()
        client_id = client["clientId"]
        user_name = client["userName"]
    except AuthError as e:
        return {"error": e.message}, e.status_code

    data = request.get_json() or {}

    # Validate room ID
    if ":" in room_id:
        return {"error": "Room ID cannot contain ':' character"}, 400

    description = data.get("description")
    copy_from = data.get("copyFrom")
    allow_create = data.get("allowCreate", True)

    r = current_app.extensions["redis"]

    # Check if room exists
    room_exists = r.exists(f"room:{room_id}:current_frame")

    if not room_exists and not allow_create:
        return {
            "status": "not_found",
            "message": f"Room '{room_id}' does not exist"
        }, 404

    # Create room if needed
    created = False
    if not room_exists:
        # Room creation logic (same as before)
        # ... (initialize room, geometries, settings, etc.)
        created = True

    # Update client metadata
    client_key = f"client:{client_id}"
    r.hset(client_key, mapping={
        "currentRoom": room_id,
        "lastActivity": time.time(),
    })

    # Add client to room's client set
    r.sadd(f"room:{room_id}:clients", client_id)

    # Get room metadata
    frame_count = r.zcard(f"room:{room_id}:trajectory:indices")
    current_frame = int(r.get(f"room:{room_id}:current_frame") or 0)

    log.info(f"Client {client_id} ({user_name}) joined room: {room_id}")

    return {
        "status": "ok",
        "roomId": room_id,
        "frameCount": frame_count,
        "step": current_frame,
        "created": created,
    }
```

---

### 3.6 Update WebSocket Connect Handler

**File:** `src/zndraw/app/events.py` (MODIFY EXISTING)

```python
@socketio.on("connect")
def handle_connect(auth):
    """Handle socket connection with JWT authentication.

    Auth payload:
        {
            "token": "jwt-token-string"
        }
    """
    from flask_socketio import join_room
    from flask_socketio import ConnectionRefusedError
    from zndraw.auth import decode_jwt_token, AuthError

    sid = request.sid
    r = current_app.extensions["redis"]

    # Get JWT token from auth
    token = auth.get("token") if auth else None

    if not token:
        log.warning(f"Client {sid} connected without JWT token")
        raise ConnectionRefusedError("Authentication token required")

    # Validate JWT token
    try:
        payload = decode_jwt_token(token)
        client_id = payload["sub"]
        user_name = payload["userName"]
    except AuthError as e:
        log.error(f"Client {sid} authentication failed: {e.message}")
        raise ConnectionRefusedError(e.message)

    # Verify client exists in Redis
    client_key = f"client:{client_id}"
    if not r.exists(client_key):
        log.error(f"Client {client_id} not found in Redis")
        raise ConnectionRefusedError("Client not registered. Call /api/login first.")

    # Update client's current SID
    r.hset(client_key, mapping={
        "currentSid": sid,
        "lastActivity": time.time(),
    })

    # Register SID -> clientId mapping
    r.set(f"sid:{sid}", client_id)

    # Get client's current room (if any)
    current_room = r.hget(client_key, "currentRoom")

    if current_room:
        # Rejoin room after reconnection
        join_room(f"room:{current_room}")
        join_room(f"user:{user_name}")
        log.info(f"Client {client_id} ({user_name}) reconnected to room {current_room} (sid: {sid})")
    else:
        log.info(f"Client {client_id} ({user_name}) connected but not in any room yet (sid: {sid})")

    return {"status": "ok", "clientId": client_id}
```

---

### 3.7 Update Python Client

**File:** `src/zndraw/zndraw.py` (MODIFY)

```python
class ZnDraw(MutableSequence):
    def __init__(
        self,
        url: str = "http://localhost:5000",
        room: str = "default",
        user_name: str = "PythonClient",
        # ... other params
    ):
        self.url = url.rstrip("/")
        self.room = room
        self._user_name = user_name
        self._client_id = None
        self._jwt_token = None  # NEW

        # Step 1: Login to get JWT token
        self._login()

        # Step 2: Create API manager with token
        self.api = APIManager(
            url=self.url,
            room=self.room,
            client_id=self._client_id,
            jwt_token=self._jwt_token  # NEW
        )

        # Step 3: Join room (authenticated)
        response_data = self.api.join_room(
            user_id=self._user_name,
            # ... other params
        )

        # Step 4: Create socket manager with token
        self.socket = SocketManager(
            zndraw_instance=self,
            jwt_token=self._jwt_token  # NEW
        )
        self.connect()

    def _login(self):
        """Authenticate and get JWT token via APIManager."""
        # Login logic is now in APIManager for better separation of concerns
        data = self.api.login(user_name=self._user_name)
        self._jwt_token = data["token"]
        self._client_id = data["clientId"]
        log.info(f"Logged in as {self._user_name} (client: {self._client_id})")
```

**File:** `src/zndraw/api_manager.py` (MODIFY)

```python
class APIManager:
    def __init__(
        self,
        url: str,
        room: str,
        client_id: str,
        jwt_token: str | None = None,  # NEW - optional during init
    ):
        self.url = url
        self.room = room
        self.client_id = client_id
        self.jwt_token = jwt_token  # NEW

    def login(self, user_name: str) -> dict:
        """Authenticate and get JWT token.

        Parameters
        ----------
        user_name : str
            Username for authentication

        Returns
        -------
        dict
            {"token": str, "clientId": str}

        Raises
        ------
        RuntimeError
            If login fails
        """
        response = requests.post(
            f"{self.url}/api/login",
            json={"userName": user_name}
        )

        if response.status_code != 200:
            raise RuntimeError(f"Login failed: {response.text}")

        data = response.json()
        # Update internal state
        self.jwt_token = data["token"]
        self.client_id = data["clientId"]
        return data

    def join_room(self, user_id: str, description=None, copy_from=None):
        """Join room with JWT authentication."""
        payload = {}
        if description:
            payload["description"] = description
        if copy_from:
            payload["copyFrom"] = copy_from

        headers = {
            "Authorization": f"Bearer {self.jwt_token}"  # NEW
        }

        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/join",
            json=payload,
            headers=headers
        )

        if response.status_code != 200:
            raise RuntimeError(f"Join failed: {response.text}")

        return response.json()
```

**File:** `src/zndraw/socket_manager.py` (MODIFY)

```python
class SocketManager:
    def __init__(self, zndraw_instance: "ZnDraw", jwt_token: str):  # MODIFIED
        self.zndraw = zndraw_instance
        self.jwt_token = jwt_token  # NEW
        self.sio = socketio.Client()
        self._register_handlers()

    def connect(self):
        """Connect to server with JWT authentication."""
        if self.sio.connected:
            print("Already connected.")
            return

        # Connect with JWT token
        self.sio.connect(
            self.zndraw.url,
            auth={"token": self.jwt_token},  # MODIFIED
            wait=True
        )
```

---

### 3.8 Update Browser Client

**File:** `app/src/utils/auth.ts` (NEW)

```typescript
/**
 * JWT Authentication utilities for browser client.
 *
 * AUTO-LOGIN STRATEGY:
 * - On first visit, generate a random UUID username
 * - Call /api/login automatically to get JWT token
 * - Store both username and token in localStorage
 * - This avoids username collisions while maintaining seamless UX
 */

const TOKEN_KEY = 'zndraw_jwt_token';
const CLIENT_ID_KEY = 'zndraw_client_id';
const USERNAME_KEY = 'zndraw_username';

export interface LoginResponse {
  status: string;
  token: string;
  clientId: string;
}

/**
 * Get or generate a unique username for this browser.
 * Generates a UUID-based username on first use.
 */
function getOrCreateUsername(): string {
  let username = localStorage.getItem(USERNAME_KEY);

  if (!username) {
    // Generate UUID-based username: "user-abc123..."
    const uuid = crypto.randomUUID();
    username = `user-${uuid.slice(0, 8)}`;
    localStorage.setItem(USERNAME_KEY, username);
    console.log('Generated new username:', username);
  }

  return username;
}

/**
 * Login and get JWT token from server.
 * If userName is not provided, uses the stored/generated username.
 */
export async function login(userName?: string): Promise<LoginResponse> {
  const finalUserName = userName || getOrCreateUsername();

  const response = await fetch('/api/login', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ userName: finalUserName }),
  });

  if (!response.ok) {
    throw new Error(`Login failed: ${response.statusText}`);
  }

  const data = await response.json();

  // Store token and client ID in localStorage
  localStorage.setItem(TOKEN_KEY, data.token);
  localStorage.setItem(CLIENT_ID_KEY, data.clientId);
  if (userName) {
    // If explicit username provided, update stored username
    localStorage.setItem(USERNAME_KEY, userName);
  }

  console.log('Logged in successfully:', data.clientId);

  return data;
}

/**
 * Ensure user is authenticated. Auto-login if needed.
 * Call this before making authenticated requests.
 */
export async function ensureAuthenticated(): Promise<void> {
  const token = getToken();
  if (!token) {
    await login(); // Auto-login with generated username
  }
}

/**
 * Get stored JWT token.
 */
export function getToken(): string | null {
  return localStorage.getItem(TOKEN_KEY);
}

/**
 * Get stored client ID.
 */
export function getClientId(): string | null {
  return localStorage.getItem(CLIENT_ID_KEY);
}

/**
 * Clear stored authentication data.
 */
export function logout(): void {
  localStorage.removeItem(TOKEN_KEY);
  localStorage.removeItem(CLIENT_ID_KEY);
  console.log('Logged out');
}

/**
 * Check if user is authenticated (has token).
 */
export function isAuthenticated(): boolean {
  return getToken() !== null;
}
```

**File:** `app/src/myapi/client.ts` (MODIFY)

```typescript
import axios from "axios";
import { decode } from "@msgpack/msgpack";
import { getToken } from "../utils/auth";  // NEW

// ... existing code ...

const apiClient = axios.create({});

// Add interceptor to include JWT token in all requests
apiClient.interceptors.request.use((config) => {
  const token = getToken();
  if (token) {
    config.headers['Authorization'] = `Bearer ${token}`;
  }
  return config;
});

// ... rest of file ...
```

**File:** `app/src/socket.ts` (MODIFY)

```typescript
import { io } from "socket.io-client";
import { getToken } from "./utils/auth";  // NEW

const URL = process.env.NODE_ENV === 'production' ? undefined : 'http://localhost:5000';

export const socket = io(URL, {
  autoConnect: false,
  auth: (cb) => {
    // Send JWT token in auth on connect
    const token = getToken();
    cb({ token });
  }
});
```

**File:** `app/src/hooks/useAuth.ts` (NEW)

```typescript
/**
 * React hook for authentication.
 */
import { useEffect, useState } from 'react';
import { login, isAuthenticated, getClientId } from '../utils/auth';

export function useAuth() {
  const [isLoggedIn, setIsLoggedIn] = useState(false);
  const [clientId, setClientId] = useState<string | null>(null);

  useEffect(() => {
    // Check if already authenticated
    if (isAuthenticated()) {
      setIsLoggedIn(true);
      setClientId(getClientId());
    }
  }, []);

  const loginUser = async (userName: string) => {
    const response = await login(userName);
    setIsLoggedIn(true);
    setClientId(response.clientId);
  };

  return {
    isLoggedIn,
    clientId,
    loginUser,
  };
}
```

**Browser Client Auto-Login Flow:**

1. **First Visit:**
   - No token exists in localStorage
   - Generate UUID-based username (e.g., "user-abc12345")
   - Call `/api/login` with generated username
   - Store token, clientId, and username in localStorage

2. **Subsequent Visits:**
   - Token exists → use it for authentication
   - No manual login required

3. **Optional Custom Username:**
   - User can change username via settings UI
   - Calls `/api/login` again with new username
   - Updates localStorage with new credentials

**Benefits:**
- ✅ No username collisions (each browser gets unique UUID-based name)
- ✅ Seamless UX (no login screen on first visit)
- ✅ User can still customize username later
- ✅ Server can track individual clients properly

---

## 4. Implementation Strategy

This is a clean implementation with no backwards compatibility concerns.
All components will be updated simultaneously to use JWT authentication.

### Backend Components
- ❌ Create `src/zndraw/auth.py` with JWT utilities
- ❌ Add `/api/login` endpoint
- ❌ Update SECRET_KEY to use environment variables
- ❌ Update `/api/rooms/{room}/join` to require JWT
- ❌ Update WebSocket `connect` handler to require JWT
- ❌ Update frame upload/modify endpoints to use JWT for client identification
- ❌ Remove all X-Client-ID header handling

### Frontend Components
- ❌ Create `app/src/utils/auth.ts` with auto-login logic
- ❌ Update `app/src/socket.ts` to send JWT in auth
- ❌ Update `app/src/myapi/client.ts` to use JWT interceptor
- ❌ Remove `app/src/utils/clientId.ts` (no longer needed)
- ❌ Add username management UI (optional)

### Python Client Components
- ❌ Add `login()` method to APIManager
- ❌ Update ZnDraw.__post_init__ to call login before join_room
- ❌ Update SocketManager to use JWT for authentication
- ❌ Remove client-side UUID generation

### Testing & Documentation
- ❌ Write unit tests for auth.py
- ❌ Write integration tests for login flow
- ❌ Test WebSocket authentication
- ❌ Update API documentation
- ❌ Security audit

**Current Status:** No implementation started. Clean slate implementation.

---

## 5. Testing Strategy

### Unit Tests

**File:** `tests/test_auth.py` (NEW)

```python
import pytest
import time
from zndraw.auth import create_jwt_token, decode_jwt_token, AuthError


def test_create_and_decode_jwt(app):
    """Test JWT creation and decoding."""
    with app.app_context():
        token = create_jwt_token("client-123", "TestUser")
        assert token is not None

        payload = decode_jwt_token(token)
        assert payload["sub"] == "client-123"
        assert payload["userName"] == "TestUser"


def test_expired_token_raises_error(app):
    """Test that expired tokens raise AuthError."""
    with app.app_context():
        # Set expiration to past
        app.config["JWT_EXPIRATION_HOURS"] = -1
        token = create_jwt_token("client-123", "TestUser")

        # Reset to normal
        app.config["JWT_EXPIRATION_HOURS"] = 168

        with pytest.raises(AuthError, match="expired"):
            decode_jwt_token(token)


def test_invalid_token_raises_error(app):
    """Test that invalid tokens raise AuthError."""
    with app.app_context():
        with pytest.raises(AuthError, match="Invalid token"):
            decode_jwt_token("invalid.token.here")
```

**File:** `tests/test_login.py` (NEW)

```python
def test_login_endpoint(server):
    """Test login endpoint returns JWT."""
    response = requests.post(
        f"{server}/api/login",
        json={"userName": "TestUser"}
    )

    assert response.status_code == 200
    data = response.json()

    assert data["status"] == "ok"
    assert "token" in data
    assert "clientId" in data
    assert "expiresIn" in data


def test_login_requires_username(server):
    """Test login fails without username."""
    response = requests.post(
        f"{server}/api/login",
        json={}
    )

    assert response.status_code == 400


def test_join_room_requires_auth(server):
    """Test joining room requires JWT."""
    response = requests.post(
        f"{server}/api/rooms/test/join",
        json={}
    )

    assert response.status_code == 401
    assert "token" in response.json()["error"].lower()


def test_join_room_with_valid_token(server):
    """Test joining room with valid JWT."""
    # Login first
    login_resp = requests.post(
        f"{server}/api/login",
        json={"userName": "TestUser"}
    )
    token = login_resp.json()["token"]

    # Join room with token
    response = requests.post(
        f"{server}/api/rooms/test/join",
        json={},
        headers={"Authorization": f"Bearer {token}"}
    )

    assert response.status_code == 200
    assert response.json()["status"] == "ok"
```