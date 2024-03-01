from datetime import datetime

import znframe
from sqlalchemy import (
    JSON,
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    String,
)
from sqlalchemy.orm import declarative_base, relationship

Base = declarative_base()


class Room(Base):
    __tablename__ = "rooms"

    token = Column(String, primary_key=True)
    currentStep = Column(Integer)
    points = Column(JSON)
    selection = Column(JSON)
    camera = Column(JSON)

    web_clients = relationship("WebClient", back_populates="room")
    frames = relationship("Frame", back_populates="room", cascade="all,delete")
    modifiers = relationship("Modifier", back_populates="room")
    bookmarks = relationship("Bookmark", back_populates="room", cascade="all,delete")
    launched_jobs = relationship("QueueItem", back_populates="room")

    def __repr__(self):
        return f"<Room(token={self.token}, currentStep={self.currentStep})>"


class Bookmark(Base):
    __tablename__ = "bookmarks"

    id = Column(Integer, primary_key=True, autoincrement=True)
    step = Column(Integer)
    text = Column(String)
    room_token = Column(String, ForeignKey("rooms.token"))

    room = relationship("Room", back_populates="bookmarks")


class WebClient(Base):
    __tablename__ = "web_clients"

    sid = Column(String, primary_key=True)
    name = Column(String)
    host = Column(Boolean, default=False)
    connected_at = Column(DateTime, default=datetime.now)
    disconnected_at = Column(DateTime, nullable=True, default=None)

    room_token = Column(String, ForeignKey("rooms.token"))
    room = relationship("Room", back_populates="web_clients")

    camera_controller_sid = Column(String, ForeignKey("web_clients.sid"), nullable=True)
    camera_controller = relationship(
        "WebClient",
        remote_side=[sid],
        uselist=False,
        foreign_keys=[camera_controller_sid],
    )

    step_controller_sid = Column(String, ForeignKey("web_clients.sid"), nullable=True)
    step_controller = relationship(
        "WebClient",
        remote_side=[sid],
        uselist=False,
        foreign_keys=[step_controller_sid],
    )


class Frame(Base):
    __tablename__ = "frames"

    id = Column(Integer, primary_key=True, autoincrement=True)
    index = Column(Integer)
    data = Column(JSON)
    room_token = Column(String, ForeignKey("rooms.token"))

    room = relationship("Room", back_populates="frames")

    def to_frame(self):
        return znframe.Frame.from_dict(self.data)

    def __repr__(self):
        return (
            f"<Frame(id={self.id}, index={self.index}, room_token={self.room_token})>"
        )


class Modifier(Base):
    __tablename__ = "modifiers"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    schema = Column(JSON)
    room_token = Column(String, ForeignKey("rooms.token"))
    room = relationship("Room", back_populates="modifiers")

    modifier_clients = relationship("ModifierClient", back_populates="modifier")


class ModifierClient(Base):
    __tablename__ = "modifier_clients"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sid = Column(String)
    timeout = Column(Float)
    available = Column(Boolean)

    modifier_id = Column(Integer, ForeignKey("modifiers.id"))
    modifier = relationship("Modifier", back_populates="modifier_clients")


class Queue(Base):
    __tablename__ = "queues"

    name = Column(String, primary_key=True)

    jobs = relationship("QueueItem", back_populates="queue")


class QueueItem(Base):
    __tablename__ = "queue_items"

    id = Column(Integer, primary_key=True, autoincrement=True)
    job_name = Column(String)
    job_id = Column(String)
    datetime = Column(String)
    status = Column(String, default="queued")
    parameters = Column(JSON)

    queue_name = Column(String, ForeignKey("queues.name"))
    queue = relationship("Queue", back_populates="jobs")
    room_token = Column(String, ForeignKey("rooms.token"))
    room = relationship("Room", back_populates="launched_jobs")
