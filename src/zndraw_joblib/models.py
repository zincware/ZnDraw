# src/zndraw_joblib/models.py
from datetime import datetime, timedelta, timezone
from enum import Enum
from typing import Any, Optional
from uuid import UUID, uuid4

from sqlalchemy import (
    Boolean,
    DateTime,
    ForeignKey,
    Index,
    String,
    Text,
    UniqueConstraint,
    func,
)
from sqlalchemy.orm import Mapped, mapped_column, relationship
from sqlalchemy.types import JSON
from zndraw_auth import Base


class TaskStatus(str, Enum):
    PENDING = "pending"
    CLAIMED = "claimed"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


TERMINAL_STATUSES = {TaskStatus.COMPLETED, TaskStatus.FAILED, TaskStatus.CANCELLED}


class WorkerJobLink(Base):
    """Bare M:N link between Worker and Job."""

    __tablename__ = "worker_job_link"

    worker_id: Mapped[UUID] = mapped_column(
        ForeignKey("worker.id", ondelete="CASCADE"), primary_key=True
    )
    job_id: Mapped[UUID] = mapped_column(
        ForeignKey("job.id", ondelete="CASCADE"), primary_key=True
    )


class Job(Base):
    __tablename__ = "job"
    __table_args__ = (
        UniqueConstraint("room_id", "category", "name", name="unique_job"),
    )

    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)
    room_id: Mapped[str] = mapped_column(String, index=True)
    category: Mapped[str] = mapped_column(String, index=True)
    name: Mapped[str] = mapped_column(String, index=True)
    schema_: Mapped[dict[str, Any]] = mapped_column("schema", JSON, default=dict)
    deleted: Mapped[bool] = mapped_column(Boolean, default=False, index=True)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        server_default=func.now(),
        index=True,
    )

    # Relationships
    tasks: Mapped[list["Task"]] = relationship(back_populates="job")
    workers: Mapped[list["Worker"]] = relationship(
        back_populates="jobs", secondary="worker_job_link", passive_deletes=True
    )

    @property
    def full_name(self) -> str:
        return f"{self.room_id}:{self.category}:{self.name}"


class Worker(Base):
    __tablename__ = "worker"

    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)
    user_id: Mapped[UUID] = mapped_column(
        ForeignKey("user.id", ondelete="CASCADE"), index=True
    )
    last_heartbeat: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=lambda: datetime.now(timezone.utc), index=True
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        server_default=func.now(),
        index=True,
    )

    # Relationships
    jobs: Mapped[list[Job]] = relationship(
        back_populates="workers", secondary="worker_job_link", passive_deletes=True
    )
    tasks: Mapped[list["Task"]] = relationship(back_populates="worker")
    providers: Mapped[list["ProviderRecord"]] = relationship(back_populates="worker")

    def is_alive(self, threshold: timedelta) -> bool:
        return datetime.now(timezone.utc) - self.last_heartbeat < threshold


class ProviderRecord(Base):
    __tablename__ = "provider"
    __table_args__ = (
        UniqueConstraint("room_id", "category", "name", name="unique_provider"),
    )

    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)
    room_id: Mapped[str] = mapped_column(String, index=True)
    category: Mapped[str] = mapped_column(String, index=True)
    name: Mapped[str] = mapped_column(String, index=True)
    schema_: Mapped[dict[str, Any]] = mapped_column("schema", JSON, default=dict)
    content_type: Mapped[str] = mapped_column(String, default="application/json")
    user_id: Mapped[UUID] = mapped_column(
        ForeignKey("user.id", ondelete="CASCADE"), index=True
    )
    worker_id: Mapped[UUID] = mapped_column(
        ForeignKey("worker.id", ondelete="CASCADE"), index=True
    )
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
        server_default=func.now(),
        index=True,
    )

    # Relationships
    worker: Mapped[Optional["Worker"]] = relationship(back_populates="providers")

    @property
    def full_name(self) -> str:
        return f"{self.room_id}:{self.category}:{self.name}"


class Task(Base):
    __tablename__ = "task"
    __table_args__ = (
        Index("ix_task_job_status_created", "job_id", "status", "created_at"),
    )

    id: Mapped[UUID] = mapped_column(primary_key=True, default=uuid4)

    job_id: Mapped[UUID] = mapped_column(ForeignKey("job.id"), index=True)
    job: Mapped[Optional[Job]] = relationship(back_populates="tasks")

    worker_id: Mapped[Optional[UUID]] = mapped_column(
        ForeignKey("worker.id"), default=None, index=True, nullable=True
    )
    worker: Mapped[Optional[Worker]] = relationship(back_populates="tasks")

    room_id: Mapped[str] = mapped_column(String, index=True)
    created_by_id: Mapped[Optional[UUID]] = mapped_column(
        ForeignKey("user.id"), default=None, index=True, nullable=True
    )

    payload: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    status: Mapped[TaskStatus] = mapped_column(default=TaskStatus.PENDING, index=True)

    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=lambda: datetime.now(timezone.utc)
    )
    started_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), default=None, nullable=True
    )
    completed_at: Mapped[Optional[datetime]] = mapped_column(
        DateTime(timezone=True), default=None, nullable=True
    )
    error: Mapped[Optional[str]] = mapped_column(Text, default=None, nullable=True)
