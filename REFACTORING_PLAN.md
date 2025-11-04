# ZnDraw Redis Architecture Refactoring Plan

**Version:** 1.0
**Date:** 2025-11-04
**Status:** Planning Phase
**Goal:** Incremental migration from direct Redis access to abstracted, future-proof hybrid storage architecture

## Executive Summary

This document outlines a phased refactoring plan to improve the ZnDraw codebase's Redis usage patterns. The refactoring addresses:

1. **Structure & Maintainability**: Eliminate 164+ f-string key constructions across 17 files
2. **Type Safety**: Remove manual decode/encode operations scattered throughout
3. **Future-Proofing**: Enable gradual migration of historical/persistent data to SQL without touching routes
4. **Proactive Architecture**: Prepare for hybrid storage (Redis + PostgreSQL) while maintaining current Redis-only functionality

### Key Metrics
- **Current State**: 164 f-string key constructions across 17 files
- **Service Layer Coverage**: ~20% (2 services: RoomService, UserService)
- **Redis Key Helpers**: 2 dataclasses (ExtensionKeys, FilesystemKeys)
- **Target**: 100% abstracted access through repositories

## Design Principles

Following CLAUDE.md requirements:
- **KISS, DRY, SOLID, YAGNI**: Keep refactoring simple, eliminate duplication
- **No Breaking Changes**: This is a new application, break backwards compatibility freely
- **Incremental Migration**: One domain at a time, not a big-bang rewrite
- **Test-First**: Leverage existing test coverage to prevent regressions
- **Collections.abc**: Implement standard interfaces where sensible

## Architecture Overview

### Current Architecture (Problematic)
```
┌─────────────────────────────────────┐
│   Routes/Events (17 files)          │
│   - Direct Redis access (80%)       │
│   - Manual key construction         │
│   - Type coercion everywhere        │
└──────────────┬──────────────────────┘
               │
               ↓
         ┌─────────┐
         │  Redis  │
         └─────────┘
```

### Target Architecture (Phase 3)
```
┌─────────────────────────────────────────┐
│         Routes/Events Layer             │  ← Flask/SocketIO
│         (No direct storage access)      │
├─────────────────────────────────────────┤
│         Service Layer                   │  ← Business Logic
│  - RoomService, UserService             │
│  - ExtensionService, FrameService       │
│  - GeometryService, LockService         │
├──────────────────┬─────────────────────┤
│   Repository     │   Models            │  ← Data Access
│   Layer          │   (Type-safe)       │
│  - RoomRepo      │   - Room            │
│  - UserRepo      │   - User            │
│  - ExtensionRepo │   - Extension       │
│  - FrameRepo     │   - Geometry        │
│  - GeometryRepo  │                     │
├──────────────────┴─────────────────────┤
│         RedisKeys (Centralized)         │  ← Key Construction
│  - RoomKeys, UserKeys, SessionKeys      │
│  - GeometryKeys, FrameKeys              │
├──────────────────┬──────────────────────┤
│  Storage Layer   │  Cache/State         │  ← Physical Storage
│  (PostgreSQL)    │  (Redis)             │
│  [Future]        │  [Current]           │
└──────────────────┴──────────────────────┘
```

## Phase 1: Foundation (Week 1-2)

### Goal
Establish clean abstractions without changing storage backend. Zero functional changes.

### Tasks

#### 1.1: Expand RedisKeys Pattern
**Priority:** HIGH
**Effort:** 2-3 days
**Dependencies:** None

Expand `src/zndraw/app/redis_keys.py` to cover all Redis key patterns:

**New Key Classes:**
```python
@dataclass(frozen=True)
class RoomKeys:
    """All room-related Redis keys."""
    room_id: str

    def description(self) -> str: ...
    def locked(self) -> str: ...
    def hidden(self) -> str: ...
    def current_frame(self) -> str: ...
    def presenter_lock(self) -> str: ...
    def bookmarks(self) -> str: ...
    def geometries(self) -> str: ...
    def figures(self) -> str: ...
    def selections(self) -> str: ...
    def selection_groups(self) -> str: ...
    def active_selection_group(self) -> str: ...
    def settings(self, username: str) -> str: ...
    def lock(self, target: str) -> str: ...
    def lock_metadata(self, target: str) -> str: ...
    def trajectory_indices(self) -> str: ...
    def chat_index(self) -> str: ...
    def chat_message(self, msg_id: str) -> str: ...

@dataclass(frozen=True)
class UserKeys:
    """User-related Redis keys."""
    username: str

    def hash_key(self) -> str: ...
    def password_hash(self) -> str: ...  # Within hash
    def current_room(self) -> str: ...   # Within hash

@dataclass(frozen=True)
class SessionKeys:
    """Session mapping keys."""
    sid: str

    def username(self) -> str: ...
    def role(self) -> str: ...

@dataclass(frozen=True)
class GeometryKeys:
    """Geometry-related keys."""
    room_id: str

    def geometries(self) -> str: ...
    def figures(self) -> str: ...
    def selections(self) -> str: ...
    def selection_groups(self) -> str: ...
```

**Files to Update:**
- 17 files with f-string key construction (164 occurrences)
- Priority order:
  1. `room_routes.py` (23 occurrences)
  2. `events.py` (40 occurrences)
  3. `geometry_routes.py` (32 occurrences)
  4. Remaining 14 files

**Testing:**
- Unit tests for key construction
- Integration tests should pass unchanged (key strings identical)

**Success Criteria:**
- Zero f-string key construction in app code
- All tests pass
- No functional changes

---

#### 1.2: Create Type-Safe Models
**Priority:** HIGH
**Effort:** 2-3 days
**Dependencies:** None (can parallelize with 1.1)

Create `src/zndraw/models/` directory with domain models:

**File Structure:**
```
src/zndraw/models/
├── __init__.py
├── user.py         # User, SessionInfo
├── room.py         # Room, RoomMetadata
├── extension.py    # Extension, ExtensionSchema, WorkerState
├── geometry.py     # Geometry, Selection, SelectionGroup
├── frame.py        # FrameIndex, TrajectoryMapping
└── bookmark.py     # Bookmark
```

**Example: User Model**
```python
from dataclasses import dataclass
from typing import Optional
import time

@dataclass
class User:
    """User entity with type-safe fields."""
    username: str
    current_sid: Optional[str] = None
    current_room: Optional[str] = None
    last_activity: float = 0.0
    created_at: float = 0.0
    password_hash: Optional[str] = None
    password_salt: Optional[str] = None

    @classmethod
    def from_redis(cls, data: dict[bytes, bytes]) -> 'User':
        """Parse Redis hash data into User object."""
        return cls(
            username=data[b'userName'].decode('utf-8'),
            current_sid=data.get(b'currentSid', b'').decode('utf-8') or None,
            current_room=data.get(b'currentRoom', b'').decode('utf-8') or None,
            last_activity=float(data.get(b'lastActivity', b'0')),
            created_at=float(data.get(b'createdAt', b'0')),
            password_hash=data.get(b'passwordHash', b'').decode('utf-8') or None,
            password_salt=data.get(b'passwordSalt', b'').decode('utf-8') or None,
        )

    def to_redis(self) -> dict[str, str]:
        """Serialize to Redis hash."""
        result = {
            'userName': self.username,
            'lastActivity': str(self.last_activity),
            'createdAt': str(self.created_at),
        }
        if self.current_sid:
            result['currentSid'] = self.current_sid
        if self.current_room:
            result['currentRoom'] = self.current_room
        if self.password_hash:
            result['passwordHash'] = self.password_hash
        if self.password_salt:
            result['passwordSalt'] = self.password_salt
        return result

    @staticmethod
    def create_new(username: str) -> 'User':
        """Factory for new user."""
        return User(
            username=username,
            created_at=time.time(),
            last_activity=time.time(),
        )
```

**Example: Extension Model**
```python
from dataclasses import dataclass
from datetime import datetime
from enum import Enum

class WorkerState(Enum):
    IDLE = "idle"
    PROGRESSING = "progressing"

@dataclass
class ExtensionSchema:
    """Extension schema definition."""
    name: str
    category: str
    schema: dict
    registered_at: datetime
    worker_id: str
    scope: str  # "room:{room_id}" or "global"

    def to_json(self) -> dict:
        """Serialize for Redis storage."""
        return {
            'name': self.name,
            'category': self.category,
            'schema': self.schema,
            'registered_at': self.registered_at.isoformat(),
            'worker_id': self.worker_id,
            'scope': self.scope,
        }

    @classmethod
    def from_json(cls, data: dict) -> 'ExtensionSchema':
        """Deserialize from Redis."""
        return cls(
            name=data['name'],
            category=data['category'],
            schema=data['schema'],
            registered_at=datetime.fromisoformat(data['registered_at']),
            worker_id=data['worker_id'],
            scope=data['scope'],
        )
```

**Testing:**
- Unit tests for serialization/deserialization
- Round-trip tests (to_redis → from_redis)
- Validation edge cases (missing fields, invalid types)

**Success Criteria:**
- All models with type hints
- Serialization/deserialization tested
- No decode() calls in model usage

---

#### 1.3: Build Repository Layer
**Priority:** HIGH
**Effort:** 1-2 weeks (incremental)
**Dependencies:** 1.1 (RedisKeys), 1.2 (Models)

Create `src/zndraw/repositories/` directory with storage abstraction:

**File Structure:**
```
src/zndraw/repositories/
├── __init__.py
├── base.py             # Abstract interfaces
├── room.py             # RoomRepository
├── user.py             # UserRepository
├── extension.py        # ExtensionRepository
├── frame.py            # FrameRepository
├── geometry.py         # GeometryRepository
└── bookmark.py         # BookmarkRepository
```

**Base Repository (Storage-Agnostic Interface):**
```python
# src/zndraw/repositories/base.py
from abc import ABC, abstractmethod
from typing import Generic, TypeVar, Optional

T = TypeVar('T')

class Repository(ABC, Generic[T]):
    """Base repository interface.

    Defines storage-agnostic operations. Implementations can use
    Redis, SQL, or hybrid storage.
    """

    @abstractmethod
    def get_by_id(self, id: str) -> Optional[T]:
        """Retrieve entity by ID."""
        pass

    @abstractmethod
    def save(self, entity: T) -> None:
        """Persist entity."""
        pass

    @abstractmethod
    def delete(self, id: str) -> bool:
        """Remove entity."""
        pass

    @abstractmethod
    def exists(self, id: str) -> bool:
        """Check entity existence."""
        pass
```

**Example: RoomRepository (Redis Implementation):**
```python
# src/zndraw/repositories/room.py
from typing import Optional
from redis import Redis
from ..models.room import Room, RoomMetadata
from ..app.redis_keys import RoomKeys
from .base import Repository

class RoomRepository(Repository[Room]):
    """Repository for Room entities.

    Current implementation: Redis
    Future: Can add SQL for historical data
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def get_by_id(self, room_id: str) -> Optional[Room]:
        """Get room metadata."""
        keys = RoomKeys(room_id)

        # Check existence
        if not self.r.exists(keys.current_frame()):
            return None

        # Gather metadata
        pipe = self.r.pipeline()
        pipe.get(keys.description())
        pipe.get(keys.locked())
        pipe.get(keys.hidden())
        pipe.get(keys.current_frame())
        results = pipe.execute()

        return Room(
            room_id=room_id,
            description=results[0].decode('utf-8') if results[0] else None,
            locked=results[1].decode('utf-8') == '1' if results[1] else False,
            hidden=results[2].decode('utf-8') == '1' if results[2] else False,
            current_frame=int(results[3]) if results[3] else 0,
        )

    def save(self, room: Room) -> None:
        """Save room metadata."""
        keys = RoomKeys(room.room_id)

        pipe = self.r.pipeline()
        if room.description:
            pipe.set(keys.description(), room.description)
        pipe.set(keys.locked(), '1' if room.locked else '0')
        pipe.set(keys.hidden(), '1' if room.hidden else '0')
        pipe.set(keys.current_frame(), room.current_frame)
        pipe.execute()

    def exists(self, room_id: str) -> bool:
        """Check if room exists."""
        return self.r.exists(RoomKeys(room_id).current_frame()) > 0

    def delete(self, room_id: str) -> bool:
        """Delete room (would need to delete all keys)."""
        # Implementation omitted for brevity
        pass

    def get_frame_count(self, room_id: str) -> int:
        """Get number of frames in room trajectory."""
        return self.r.zcard(RoomKeys(room_id).trajectory_indices())

    def set_current_frame(self, room_id: str, frame: int) -> None:
        """Update current frame position."""
        self.r.set(RoomKeys(room_id).current_frame(), frame)

    def set_locked(self, room_id: str, locked: bool) -> None:
        """Update room locked status."""
        self.r.set(RoomKeys(room_id).locked(), '1' if locked else '0')
```

**Example: ExtensionRepository (Your Use Case):**
```python
# src/zndraw/repositories/extension.py
from typing import Optional
from redis import Redis
from ..models.extension import ExtensionSchema, WorkerState
from ..app.redis_keys import ExtensionKeys
from .base import Repository

class ExtensionRepository(Repository[ExtensionSchema]):
    """Repository for Extension schemas and workers.

    Current: Redis-only (active state)
    Future: Add SQL backend for historical queries
    """

    def __init__(self, redis_client: Redis, db_session=None):
        self.r = redis_client
        self.db = db_session  # Optional SQL backend

    def register_extension(
        self,
        scope: str,  # "room:{room_id}" or "global"
        extension: ExtensionSchema
    ) -> None:
        """Register extension and track in history.

        Active state → Redis (fast queries)
        Historical record → SQL (if configured)
        """
        # Parse scope
        if scope.startswith("room:"):
            room_id = scope.split(":")[1]
            keys = ExtensionKeys.for_extension(room_id, extension.category, extension.name)
            schema_key = ExtensionKeys.schema_key(room_id, extension.category)
        else:
            keys = ExtensionKeys.for_global_extension(extension.category, extension.name)
            schema_key = ExtensionKeys.global_schema_key(extension.category)

        # Store active schema in Redis
        self.r.hset(schema_key, extension.name, extension.to_json())

        # Add worker to idle pool
        self.r.sadd(keys.idle_workers, extension.worker_id)

        # Historical tracking (if SQL enabled)
        if self.db:
            from ..db.models import ExtensionRegistration
            self.db.add(ExtensionRegistration(
                scope=scope,
                category=extension.category,
                name=extension.name,
                schema_json=extension.schema,
                worker_id=extension.worker_id,
                registered_at=extension.registered_at,
            ))
            self.db.commit()

    def get_history(self, name: str) -> list[ExtensionSchema]:
        """Get historical registrations for an extension.

        Requires SQL backend to be configured.
        """
        if not self.db:
            raise NotImplementedError(
                "Extension history requires SQL backend. "
                "Configure database session to enable."
            )

        from ..db.models import ExtensionRegistration
        records = self.db.query(ExtensionRegistration)\
            .filter_by(name=name)\
            .order_by(ExtensionRegistration.registered_at.desc())\
            .all()

        return [ExtensionSchema.from_db_record(r) for r in records]

    def acquire_worker(
        self,
        scope: str,
        category: str,
        name: str
    ) -> Optional[str]:
        """Atomically move worker from idle to progressing."""
        if scope.startswith("room:"):
            room_id = scope.split(":")[1]
            keys = ExtensionKeys.for_extension(room_id, category, name)
        else:
            keys = ExtensionKeys.for_global_extension(category, name)

        # Atomic move with pipeline
        worker_id = self.r.spop(keys.idle_workers)
        if worker_id:
            self.r.sadd(keys.progressing_workers, worker_id)
        return worker_id.decode('utf-8') if worker_id else None

    def release_worker(
        self,
        scope: str,
        category: str,
        name: str,
        worker_id: str
    ) -> None:
        """Return worker to idle pool."""
        if scope.startswith("room:"):
            room_id = scope.split(":")[1]
            keys = ExtensionKeys.for_extension(room_id, category, name)
        else:
            keys = ExtensionKeys.for_global_extension(category, name)

        self.r.srem(keys.progressing_workers, worker_id)
        self.r.sadd(keys.idle_workers, worker_id)
```

**Repositories to Implement (Priority Order):**

1. **RoomRepository** (Week 1)
   - Straightforward, clear boundaries
   - Good starting point to establish pattern
   - ~300 LOC

2. **UserRepository** (Week 1)
   - Simple hash operations
   - Session management
   - ~200 LOC

3. **ExtensionRepository** (Week 2)
   - Your historical data use case
   - Most complex (worker tracking, queues)
   - SQL integration point
   - ~400 LOC

4. **GeometryRepository** (Week 2)
   - Selections, figures, geometries
   - ~250 LOC

5. **FrameRepository** (Week 2)
   - Trajectory operations (sorted sets)
   - ~300 LOC

6. **BookmarkRepository** (Week 2)
   - Simple hash operations
   - ~150 LOC

**Testing:**
- Unit tests for each repository
- Mock Redis for fast tests
- Integration tests with real Redis
- Ensure RoomService/UserService pass with new repositories

**Success Criteria:**
- All repositories implemented and tested
- Existing services refactored to use repositories
- All integration tests pass

---

#### 1.4: Update Routes to Use Repositories
**Priority:** HIGH
**Effort:** 1 week
**Dependencies:** 1.3 (Repositories)

Refactor routes to inject and use repositories instead of direct Redis access.

**Dependency Injection Pattern:**
```python
# In app/__init__.py
from flask import Flask
from .repositories import (
    RoomRepository,
    UserRepository,
    ExtensionRepository,
)

def create_app():
    app = Flask(__name__)

    # Initialize Redis
    redis_client = Redis(...)
    app.extensions["redis"] = redis_client

    # Initialize repositories
    app.extensions["room_repo"] = RoomRepository(redis_client)
    app.extensions["user_repo"] = UserRepository(redis_client)
    app.extensions["extension_repo"] = ExtensionRepository(redis_client)
    # ... etc

    return app
```

**Before (Direct Redis):**
```python
# room_routes.py
@room_bp.route("/api/room/<room_id>/description", methods=["GET"])
def get_room_description(room_id: str):
    r = current_app.extensions["redis"]
    desc = r.get(f"room:{room_id}:description")
    return {"description": desc.decode('utf-8') if desc else None}
```

**After (Repository):**
```python
# room_routes.py
@room_bp.route("/api/room/<room_id>/description", methods=["GET"])
def get_room_description(room_id: str):
    room_repo = current_app.extensions["room_repo"]
    room = room_repo.get_by_id(room_id)
    return {"description": room.description if room else None}
```

**Files to Update (Priority Order):**
1. `room_routes.py` - Use RoomRepository
2. `utility_routes.py` - Use UserRepository
3. `extension_routes.py` - Use ExtensionRepository
4. `geometry_routes.py` - Use GeometryRepository
5. `frame_routes.py` - Use FrameRepository
6. `bookmark_routes.py` - Use BookmarkRepository
7. `events.py` - Use all repositories (last, most complex)

**Testing Strategy:**
- Update one route file at a time
- Run full test suite after each file
- Use feature flags to toggle between old/new implementations if needed

**Success Criteria:**
- Zero direct Redis access in routes (no `r = current_app.extensions["redis"]`)
- All integration tests pass
- No functional changes
- Improved code readability

---

## Phase 2: Historical Data Support (Week 3-4)

### Goal
Add SQLAlchemy for historical queries without disrupting Redis fast path.

### Tasks

#### 2.1: Add SQLAlchemy Models
**Priority:** MEDIUM
**Effort:** 2-3 days
**Dependencies:** Phase 1 complete

Create `src/zndraw/db/` directory with SQLAlchemy models:

**File Structure:**
```
src/zndraw/db/
├── __init__.py
├── base.py             # Declarative base
├── models.py           # SQLAlchemy models
└── session.py          # Session management
```

**Example: Extension History Table:**
```python
# src/zndraw/db/models.py
from sqlalchemy import Column, Integer, String, DateTime, JSON, Index
from sqlalchemy.ext.declarative import declarative_base
from datetime import datetime

Base = declarative_base()

class ExtensionRegistration(Base):
    """Historical record of extension registrations.

    Enables queries like:
    - Show all versions of extension X
    - When was extension Y last registered?
    - Which worker registered extension Z?
    """
    __tablename__ = 'extension_registrations'

    id = Column(Integer, primary_key=True)
    scope = Column(String(255), nullable=False, index=True)  # "room:xyz" or "global"
    category = Column(String(100), nullable=False, index=True)
    name = Column(String(255), nullable=False, index=True)
    schema_json = Column(JSON, nullable=False)
    worker_id = Column(String(255), nullable=False)
    registered_at = Column(DateTime, nullable=False, default=datetime.utcnow, index=True)

    __table_args__ = (
        Index('idx_extension_history', 'name', 'registered_at'),
        Index('idx_scope_category', 'scope', 'category'),
    )

    def __repr__(self):
        return f"<ExtensionRegistration(name={self.name}, registered_at={self.registered_at})>"

class RoomEvent(Base):
    """Historical record of room lifecycle events.

    Tracks: creation, deletion, lock changes, etc.
    """
    __tablename__ = 'room_events'

    id = Column(Integer, primary_key=True)
    room_id = Column(String(255), nullable=False, index=True)
    event_type = Column(String(50), nullable=False)  # "created", "deleted", "locked", etc.
    user = Column(String(255), nullable=True)
    metadata = Column(JSON, nullable=True)
    timestamp = Column(DateTime, nullable=False, default=datetime.utcnow, index=True)

    __table_args__ = (
        Index('idx_room_events', 'room_id', 'timestamp'),
    )

class UserLoginEvent(Base):
    """Track user authentication events."""
    __tablename__ = 'user_login_events'

    id = Column(Integer, primary_key=True)
    username = Column(String(255), nullable=False, index=True)
    session_id = Column(String(255), nullable=False)
    login_at = Column(DateTime, nullable=False, default=datetime.utcnow, index=True)
    logout_at = Column(DateTime, nullable=True)
    ip_address = Column(String(45), nullable=True)  # IPv6 compatible
```

**Database Configuration:**
```python
# src/zndraw/db/session.py
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from .models import Base

class DatabaseManager:
    """Manages SQLAlchemy session lifecycle."""

    def __init__(self, database_url: str, echo: bool = False):
        self.engine = create_engine(database_url, echo=echo)
        self.session_factory = sessionmaker(bind=self.engine)
        self.Session = scoped_session(self.session_factory)

    def create_tables(self):
        """Create all tables."""
        Base.metadata.create_all(self.engine)

    def drop_tables(self):
        """Drop all tables (testing only)."""
        Base.metadata.drop_all(self.engine)

    def get_session(self):
        """Get thread-local session."""
        return self.Session()
```

**Dependencies:**
```toml
# pyproject.toml
[project]
dependencies = [
    "sqlalchemy>=2.0.0",
    "psycopg2-binary>=2.9.0",  # PostgreSQL driver
    # ... existing
]
```

**Testing:**
- Unit tests for model creation
- Migration scripts
- Integration tests with in-memory SQLite

**Success Criteria:**
- SQLAlchemy models created
- Migrations work
- Optional: Can run app with or without SQL backend

---

#### 2.2: Update Repositories for Hybrid Storage
**Priority:** MEDIUM
**Effort:** 3-4 days
**Dependencies:** 2.1 (SQL Models)

Update repositories to dual-write to Redis (active state) and SQL (history).

**Example: ExtensionRepository Hybrid Mode:**
```python
class ExtensionRepository:
    """Hybrid storage: Redis for active, SQL for history."""

    def __init__(self, redis_client: Redis, db_session=None):
        self.r = redis_client
        self.db = db_session
        self._history_enabled = db_session is not None

    def register_extension(self, scope: str, extension: ExtensionSchema) -> None:
        """Register extension with dual-write."""
        # Always write to Redis (active state)
        self._write_to_redis(scope, extension)

        # Optionally write to SQL (history)
        if self._history_enabled:
            self._write_to_sql(scope, extension)

    def _write_to_redis(self, scope: str, extension: ExtensionSchema) -> None:
        """Active state in Redis."""
        # ... existing Redis logic

    def _write_to_sql(self, scope: str, extension: ExtensionSchema) -> None:
        """Historical record in SQL."""
        from ..db.models import ExtensionRegistration
        self.db.add(ExtensionRegistration(
            scope=scope,
            category=extension.category,
            name=extension.name,
            schema_json=extension.schema,
            worker_id=extension.worker_id,
            registered_at=extension.registered_at,
        ))
        self.db.commit()

    def get_history(
        self,
        name: str,
        limit: int = 100,
        since: Optional[datetime] = None
    ) -> list[ExtensionSchema]:
        """Query historical registrations."""
        if not self._history_enabled:
            raise NotImplementedError("History requires SQL backend")

        from ..db.models import ExtensionRegistration
        query = self.db.query(ExtensionRegistration).filter_by(name=name)

        if since:
            query = query.filter(ExtensionRegistration.registered_at >= since)

        records = query.order_by(
            ExtensionRegistration.registered_at.desc()
        ).limit(limit).all()

        return [ExtensionSchema.from_db_record(r) for r in records]
```

**Configuration:**
```python
# app/__init__.py
def create_app(enable_sql: bool = False):
    app = Flask(__name__)

    redis_client = Redis(...)
    app.extensions["redis"] = redis_client

    # Optional SQL backend
    db_session = None
    if enable_sql:
        from .db.session import DatabaseManager
        db_manager = DatabaseManager(app.config['DATABASE_URL'])
        db_session = db_manager.get_session()
        app.extensions["db"] = db_session

    # Repositories with optional SQL
    app.extensions["extension_repo"] = ExtensionRepository(
        redis_client,
        db_session
    )
```

**Testing:**
- Test Redis-only mode (existing functionality)
- Test hybrid mode (dual-write)
- Test SQL queries for history
- Verify Redis performance not impacted

**Success Criteria:**
- Dual-write working
- History queries functional
- App works with SQL disabled (backwards compatible)

---

#### 2.3: Add History Query Endpoints
**Priority:** LOW
**Effort:** 1-2 days
**Dependencies:** 2.2 (Hybrid Repositories)

Add new API endpoints for historical queries:

```python
# extension_routes.py
@extension_bp.route("/api/extensions/<name>/history", methods=["GET"])
def get_extension_history(name: str):
    """Get historical registrations for an extension."""
    extension_repo = current_app.extensions["extension_repo"]

    # Parse query params
    limit = request.args.get("limit", 100, type=int)
    since = request.args.get("since", None)  # ISO format

    try:
        history = extension_repo.get_history(
            name=name,
            limit=limit,
            since=datetime.fromisoformat(since) if since else None
        )
        return {
            "name": name,
            "registrations": [e.to_json() for e in history],
            "count": len(history),
        }
    except NotImplementedError:
        return {"error": "History feature requires SQL backend"}, 501
```

**New Endpoints:**
- `GET /api/extensions/<name>/history` - Extension registration history
- `GET /api/rooms/<room_id>/events` - Room lifecycle events
- `GET /api/users/<username>/sessions` - User login history

---

## Phase 3: Performance & Polish (Week 5+)

### Goal
Optimize patterns and add convenience features.

### Tasks

#### 3.1: Maintain Explicit Indices
**Priority:** MEDIUM
**Effort:** 2-3 days

Replace `scan_iter()` with maintained sets:

```python
class RoomRepository:
    def create(self, room: Room) -> None:
        """Create room and add to index."""
        pipe = self.r.pipeline()
        # ... create room keys
        pipe.sadd("all_rooms", room.room_id)  # Maintain index
        pipe.execute()

    def list_all(self) -> list[str]:
        """O(N) in rooms, not keyspace."""
        return [r.decode('utf-8') for r in self.r.smembers("all_rooms")]

    def delete(self, room_id: str) -> None:
        """Delete room and remove from index."""
        # ... delete room keys
        self.r.srem("all_rooms", room_id)
```

**Apply to:**
- Rooms (`all_rooms` set)
- Users (`all_users` set)
- Global extensions
- Global filesystems

---

#### 3.2: Create Service Layer for Complex Operations
**Priority:** MEDIUM
**Effort:** 1 week

Build high-level services that orchestrate repositories:

```python
# src/zndraw/services/extension_service.py
class ExtensionService:
    """High-level extension management service."""

    def __init__(
        self,
        extension_repo: ExtensionRepository,
        user_repo: UserRepository
    ):
        self.extensions = extension_repo
        self.users = user_repo

    def register_extension_for_user(
        self,
        username: str,
        scope: str,
        category: str,
        name: str,
        schema: dict
    ) -> ExtensionSchema:
        """Register extension and link to user."""
        user = self.users.get_by_name(username)
        if not user:
            raise ValueError(f"User {username} not found")

        extension = ExtensionSchema(
            name=name,
            category=category,
            schema=schema,
            registered_at=datetime.now(),
            worker_id=user.current_sid,
            scope=scope,
        )

        self.extensions.register_extension(scope, extension)
        return extension

    def cleanup_user_extensions(self, username: str) -> int:
        """Remove all extensions for disconnected user."""
        # Orchestrate cleanup across multiple repositories
        count = 0
        # ... implementation
        return count
```

**Services to Create:**
- `ExtensionService` - Orchestrate extension lifecycle
- `LockService` - Centralized lock management
- `SessionService` - User session lifecycle

---

#### 3.3: Add Comprehensive Logging
**Priority:** LOW
**Effort:** 2 days

Add structured logging for Redis operations:

```python
import logging
import time

logger = logging.getLogger(__name__)

class RoomRepository:
    def get_by_id(self, room_id: str) -> Optional[Room]:
        start = time.time()
        result = self._fetch_from_redis(room_id)
        elapsed = (time.time() - start) * 1000

        logger.info(
            "redis.get_room",
            extra={
                "room_id": room_id,
                "hit": result is not None,
                "latency_ms": elapsed,
            }
        )
        return result
```

---

## Migration Strategy

### Rollout Plan

**Week 1:** Phase 1.1 + 1.2 (RedisKeys + Models)
- Low risk, zero functional changes
- Can deploy incrementally

**Week 2:** Phase 1.3 (Repositories)
- Build all repositories
- Keep existing routes unchanged
- Parallel testing

**Week 3:** Phase 1.4 (Update Routes)
- Update routes one file at a time
- Deploy and test after each file
- Use feature flags if needed

**Week 4:** Phase 2.1 + 2.2 (SQL Integration)
- Add SQL models
- Enable hybrid mode (optional)
- Test in staging

**Week 5+:** Phase 3 (Polish)
- Performance optimizations
- Service layer
- Cleanup

### Rollback Strategy

Each phase is independently reversible:

1. **Phase 1.1 (RedisKeys):** Revert key construction, use f-strings
2. **Phase 1.2 (Models):** Revert to decode() calls
3. **Phase 1.3 (Repositories):** Repositories wrap existing logic, no risk
4. **Phase 1.4 (Routes):** Feature flag to toggle repository usage
5. **Phase 2 (SQL):** SQL is optional, disable via config

### Testing Strategy

**Unit Tests:**
- Test each component in isolation
- Mock Redis for fast tests
- 100% coverage for new code

**Integration Tests:**
- Test with real Redis
- Test SQL migrations
- Test hybrid mode

**Regression Tests:**
- Existing test suite must pass
- No functional changes until Phase 2

**Performance Tests:**
- Benchmark key Redis operations
- Ensure no slowdown from abstraction
- Target: <5% overhead

---

## Success Metrics

### Code Quality
- [ ] Zero direct Redis access in routes/events
- [ ] Zero f-string key construction
- [ ] Zero manual decode() calls in business logic
- [ ] 100% type-annotated models
- [ ] 90%+ test coverage

### Performance
- [ ] No measurable latency increase (<5%)
- [ ] Room discovery via index (not scan)
- [ ] All Redis operations in O(1) or O(log N)

### Maintainability
- [ ] New feature can add repository in <1 day
- [ ] Storage backend swappable without touching routes
- [ ] Clear separation of concerns (Route → Service → Repository → Storage)

### Historical Data (Phase 2)
- [ ] Extension history queryable
- [ ] Room event tracking
- [ ] User session history

---

## Open Questions

1. **Database Choice:** PostgreSQL vs SQLite vs MongoDB?
   - Recommendation: PostgreSQL (production), SQLite (dev/testing)

2. **Caching Strategy:** Add in-memory cache above repositories?
   - Recommendation: Phase 3 optimization if needed

3. **Schema Migrations:** Use Alembic for SQL migrations?
   - Recommendation: Yes, add in Phase 2.1

4. **Backwards Compatibility:** Support downgrade to Redis-only?
   - Recommendation: Yes, SQL is optional via config

---

## Appendix

### File Inventory

**Current Files with Direct Redis Access (164 f-strings):**
1. events.py (40)
2. geometry_routes.py (32)
3. room_routes.py (23)
4. job_manager.py (12)
5. extension_routes.py (7)
6. route_utils.py (7)
7. frame_routes.py (6)
8. bookmark_routes.py (6)
9. file_browser.py (6)
10. room_manager.py (6)
11. chat_utils.py (5)
12. redis_keys.py (5)
13. screenshot_chat_routes.py (4)
14. utility_routes.py (2)
15. filesystem_routes.py (1)
16. queue_manager.py (1)
17. metadata_manager.py (1)

### Estimated Effort

| Phase | Tasks | Effort | Dependencies |
|-------|-------|--------|--------------|
| 1.1 | RedisKeys | 2-3 days | None |
| 1.2 | Models | 2-3 days | None |
| 1.3 | Repositories | 1-2 weeks | 1.1, 1.2 |
| 1.4 | Update Routes | 1 week | 1.3 |
| 2.1 | SQL Models | 2-3 days | Phase 1 |
| 2.2 | Hybrid Storage | 3-4 days | 2.1 |
| 2.3 | History Endpoints | 1-2 days | 2.2 |
| 3.1 | Indices | 2-3 days | Phase 1 |
| 3.2 | Services | 1 week | Phase 1 |
| 3.3 | Logging | 2 days | Phase 1 |
| **Total** | | **6-8 weeks** | |

### Glossary

- **Repository Pattern:** Data access abstraction that isolates storage implementation
- **Domain Model:** Type-safe representation of business entities
- **Service Layer:** Business logic coordination between repositories
- **Hybrid Storage:** Using multiple data stores (Redis + SQL) for different use cases
- **Dual-Write:** Writing to both Redis and SQL to maintain consistency

---

**End of Plan**
