from sqlalchemy import (
    JSON,
    Boolean,
    Column,
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

    clients = relationship("Client", back_populates="room")
    frames = relationship("Frame", back_populates="room")
    room_modifiers = relationship("RoomModifier", back_populates="room")
    bookmarks = relationship("Bookmark", back_populates="room")

    def __repr__(self):
        return f"<Room(token={self.token}, currentStep={self.currentStep})>"


class Bookmark(Base):
    __tablename__ = "bookmarks"

    id = Column(Integer, primary_key=True, autoincrement=True)
    step = Column(Integer)
    text = Column(String)
    room_token = Column(String, ForeignKey("rooms.token"))

    room = relationship("Room", back_populates="bookmarks")


class Client(Base):
    __tablename__ = "clients"

    sid = Column(String, primary_key=True)
    cameras = Column(JSON)

    room_token = Column(String, ForeignKey("rooms.token"))
    room = relationship("Room", back_populates="clients")


class Frame(Base):
    __tablename__ = "frames"

    id = Column(Integer, primary_key=True, autoincrement=True)
    index = Column(Integer)
    data = Column(JSON)
    room_token = Column(String, ForeignKey("rooms.token"))

    room = relationship("Room", back_populates="frames")

    def __repr__(self):
        return (
            f"<Frame(id={self.id}, index={self.index}, room_token={self.room_token})>"
        )


class GlobalModifier(Base):
    __tablename__ = "global_modifiers"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    schema = Column(JSON)

    global_modifier_clients = relationship(
        "GlobalModifierClient", back_populates="global_modifier"
    )


class GlobalModifierClient(Base):
    __tablename__ = "global_modifier_clients"

    sid = Column(String, primary_key=True)
    timeout = Column(Float)
    available = Column(Boolean)

    modifier = Column(String, ForeignKey("global_modifiers.id"))
    global_modifier = relationship(
        "GlobalModifier", back_populates="global_modifier_clients"
    )


class RoomModifier(Base):
    __tablename__ = "room_modifiers"

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    schema = Column(JSON)
    room_token = Column(String, ForeignKey("rooms.token"))

    room = relationship("Room", back_populates="room_modifiers")
    room_modifier_clients = relationship(
        "RoomModifierClient", back_populates="room_modifier"
    )


class RoomModifierClient(Base):
    __tablename__ = "room_modifier_clients"

    sid = Column(String, primary_key=True)
    timeout = Column(Float)
    available = Column(Boolean)

    modifier = Column(String, ForeignKey("room_modifiers.id"))
    room_modifier = relationship("RoomModifier", back_populates="room_modifier_clients")

class Queue(Base):
    __tablename__ = "queues"
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String)
    jobs = Column(String)