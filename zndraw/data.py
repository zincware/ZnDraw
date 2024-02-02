import dataclasses


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
