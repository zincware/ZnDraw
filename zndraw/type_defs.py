import typing as t

if t.TYPE_CHECKING:
    from zndraw.base import Extension


class RegisterModifier(t.TypedDict):
    cls: t.Type["Extension"]
    run_kwargs: dict
    public: bool


class TimeoutConfig(t.TypedDict):
    """Timeout configuration for the ZnDraw client."""

    connection: int
    modifier: float
    between_calls: float

    emit_retries: int
    call_retries: int
    connect_retries: int


class JupyterConfig(t.TypedDict):
    width: str | int
    height: str | int


class CameraData(t.TypedDict):
    position: list[float]
    target: list[float]


class ASEDict(t.TypedDict):
    numbers: list[int]
    positions: list[list[float]]
    connectivity: list[tuple[int, int, int]]
    arrays: dict[str, list[float | int | list[float | int]]]
    info: dict[str, float | int]
    # calc: dict[str, float|int|np.ndarray] # should this be split into arrays and info?
    pbc: list[bool]
    cell: list[list[float]]
    vectors: list[list[list[float]]]
    constraints: list[dict]


class ASEJson(t.TypedDict):
    _type: t.Literal["ase.Atoms"]
    value: ASEDict


# Type hint is string, but correctly it is 'json.dumps(ASEJson)'
ATOMS_LIKE = t.Union[ASEDict, str]
