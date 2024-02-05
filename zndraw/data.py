import dataclasses

import znframe


@dataclasses.dataclass
class RoomSetData:
    """Update the room with new data.

    Attributes
    ----------
    frames: dict[int, znframe.Frame | None]
        If the frame is None, it is deleted from the room.
    update_database: bool
        Whether to update the database with the new data.
    """

    points: list[list[float]] | None = None
    bookmarks: dict[int, str] | None = None
    step: int | None = None
    selection: list[int] | None = None
    frames: dict[int, znframe.Frame | None] | None = dataclasses.field(
        default=None, repr=False
    )

    update_database: bool = False

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


@dataclasses.dataclass
class RoomGetData:
    points: bool | list[list[float]] = False
    bookmarks: bool | dict[str, str] = False
    step: bool | int = False
    selection: bool | list[int] = False
    length: bool | int = False
    segments: bool | list[list[float]] = False
    frames: list[int] | None | list[dict] = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def get_current_state(cls):
        return cls(
            points=True,
            bookmarks=True,
            step=True,
            selection=True,
            length=True,
            segments=True,
            frames=["current"],
        )


@dataclasses.dataclass
class CeleryTaskData:
    """A message to emit a 'emit' or 'call' from the server.

    Attributes
    ----------
    target: str
        The target of the message, e.g. a room name or sid.
    event: str
        The event to emit or call, e.g. 'message:log'.
    data: dict
        The data to send with the message.
    disconnect: bool
        Whether to tell the server to  disconnect this client,
        after it has received the message. Using 'disconnect' after
        'emit' on the client might loose the message.
    timeout: int
        The timeout in seconds, when using 'call'.
    authentication: str
        Authentication token, used for ensuring that not every client
        can send arbitrary messages through the server.
    """

    target: str
    event: str
    data: dict | None | str
    disconnect: bool = False
    timeout: int = 60
    authentication: str = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


@dataclasses.dataclass
class FrameData:
    """
    index: int
        The index of the frame in the trajectory.
    data: dict
        The data of the frame.
    update: bool
        Whether the UI should be updated.
    update_database: bool
        Whether the configuration should be saved in the database
    """

    index: int
    update: bool
    data: dict
    update_database: bool


@dataclasses.dataclass
class JoinData:
    token: str
    auth_token: str


@dataclasses.dataclass
class ModifierRunData:
    params: dict
    url: dict
    sid: str = None
    target: str = None

    @property
    def name(self) -> str:
        return self.params["method"]["discriminator"]


@dataclasses.dataclass
class AnalysisRunData:
    params: dict
    target: str = None


@dataclasses.dataclass
class AnalysisFigureData:
    figure: dict
    token: str


@dataclasses.dataclass
class SceneSetData:
    index: int
    token: str


@dataclasses.dataclass
class SceneStepData:
    token: str


@dataclasses.dataclass
class AtomsDownloadData:
    token: str
    indices: list[int]


@dataclasses.dataclass
class DeleteAtomsData:
    index: list
    token: str


@dataclasses.dataclass
class AtomsLengthData:
    token: str


@dataclasses.dataclass
class SchemaData:
    schema: dict
    sid: str


@dataclasses.dataclass
class SelectionSetData:
    selection: list[int]
    token: str = None


@dataclasses.dataclass
class SelectionRunData:
    params: dict
    target: str = None


@dataclasses.dataclass
class MessageData:
    message: str
    token: str


@dataclasses.dataclass
class PlayData:
    token: str


@dataclasses.dataclass
class ModifierRegisterData:
    schema: dict
    name: str
    default: bool
    timeout: float = 60


@dataclasses.dataclass
class BookmarksSetData:
    bookmarks: dict
    token: str


@dataclasses.dataclass
class PointsSetData:
    value: list[list[float]]
    token: str = None


@dataclasses.dataclass
class ModifierRunRunningData:
    token: str
    name: str


@dataclasses.dataclass
class SubscribedUserData:
    user: str


@dataclasses.dataclass
class SceneUpdateData:
    """
    Attributes:
        camera: dict with {position: list[float], quaternion: list[float]}
        step: int
    """

    camera: dict = None
    step: int = None
