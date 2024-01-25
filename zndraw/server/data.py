import dataclasses

@dataclasses.dataclass
class CeleryTaskData:
    target: str
    event: str
    data: dict

@dataclasses.dataclass
class JoinData:
    token: str
    uuid: str
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
class AnalysisSchemaData:
    schema: dict
    sid: str


@dataclasses.dataclass
class ModifierSchemaData:
    schema: dict
    token: str


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
    uuid: str
    modifiers: list[dict]
    token: str = None

    @property
    def name(self) -> str:
        return self.modifiers[0]["name"]

    @property
    def is_default(self) -> bool:
        return self.modifiers[0]["default"]

    @property
    def schema(self) -> dict:
        return self.modifiers[0]["schema"]


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
