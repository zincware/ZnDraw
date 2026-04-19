# src/zndraw_joblib/schemas.py
from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING, Any, Generic, TypeVar
from uuid import UUID

from pydantic import BaseModel, Field

from zndraw_joblib.models import TaskStatus

if TYPE_CHECKING:
    from zndraw_joblib.models import ProviderRecord


class JobRegisterRequest(BaseModel):
    category: str
    name: str
    schema_: dict[str, Any] = Field(default={}, alias="schema")
    worker_id: UUID | None = None  # Optional: use existing worker

    model_config = {"populate_by_name": True}


class JobResponse(BaseModel):
    id: UUID
    room_id: str
    category: str
    name: str
    full_name: str
    schema_: dict[str, Any] = Field(alias="schema")
    workers: list[UUID]
    worker_id: UUID | None = (
        None  # The worker that was created/used for this registration
    )

    model_config = {"populate_by_name": True}


class JobSummary(BaseModel):
    full_name: str
    category: str
    name: str
    workers: list[UUID]


class WorkerSummary(BaseModel):
    id: UUID
    last_heartbeat: datetime
    job_count: int


class WorkerResponse(BaseModel):
    id: UUID
    last_heartbeat: datetime


class TaskSubmitRequest(BaseModel):
    payload: dict[str, Any] = Field(default_factory=dict)


class TaskResponse(BaseModel):
    id: UUID
    job_name: str
    room_id: str
    status: TaskStatus
    created_at: datetime
    started_at: datetime | None = None
    completed_at: datetime | None = None
    worker_id: UUID | None = None
    error: str | None = None
    payload: dict[str, Any] = {}
    queue_position: int | None = None


class TaskClaimRequest(BaseModel):
    worker_id: UUID


class TaskClaimResponse(BaseModel):
    task: TaskResponse | None = None


class TaskUpdateRequest(BaseModel):
    status: TaskStatus
    error: str | None = None


T = TypeVar("T")


class PaginatedResponse(BaseModel, Generic[T]):
    items: list[T]
    total: int
    limit: int
    offset: int


class ProviderRegisterRequest(BaseModel):
    category: str
    name: str
    schema_: dict[str, Any] = Field(default={}, alias="schema")
    content_type: str = "application/json"
    worker_id: UUID | None = None

    model_config = {"populate_by_name": True}


class ProviderResponse(BaseModel):
    id: UUID
    room_id: str
    category: str
    name: str
    full_name: str
    schema_: dict[str, Any] = Field(alias="schema")
    worker_id: UUID | None = None
    created_at: datetime

    model_config = {"populate_by_name": True}

    @classmethod
    def from_record(cls, record: ProviderRecord) -> ProviderResponse:
        return cls(
            id=record.id,
            room_id=record.room_id,
            category=record.category,
            name=record.name,
            full_name=record.full_name,
            schema=record.schema_,
            worker_id=record.worker_id,
            created_at=record.created_at,
        )
