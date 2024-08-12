import dataclasses
import typing as t

from zndraw.scene import Scene
from zndraw.utils import emit_with_retry

if t.TYPE_CHECKING:
    from zndraw import ZnDraw

HSLColor = t.Tuple[float, float, float]


@dataclasses.dataclass
class ArrowsConfig:
    colormap: list[HSLColor] = ((-0.5, 0.9, 0.5), (0.0, 0.9, 0.5))
    normalize: bool = True
    colorrange: tuple[float, float] = (0, 1.0)
    scale_vector_thickness: bool = False
    opacity: float = 1.0

    _vis = None  # not a dataclass field

    def __setattr__(self, name: str, value) -> None:
        super().__setattr__(name, value)
        if self._vis is not None:
            if name != "vis" and name in [x.name for x in dataclasses.fields(self)]:
                data = {
                    "arrows": dataclasses.asdict(self),
                }
                emit_with_retry(
                    self._vis.socket,
                    "room:config:set",
                    data,
                    retries=self._vis.timeout["emit_retries"],
                )
                self._vis.socket.sleep(0.1)  # maybe use call?

    def set_vis(self, vis: "ZnDraw") -> None:
        self._vis = vis


@dataclasses.dataclass
class ZnDrawConfig:
    vis: "ZnDraw|None" = dataclasses.field(repr=False)
    arrows: ArrowsConfig = dataclasses.field(default_factory=ArrowsConfig)
    scene: Scene = dataclasses.field(default_factory=Scene)

    def __post_init__(self) -> None:
        if self.vis is not None:
            self.arrows.set_vis(self.vis)
            self.scene.set_vis(self.vis)

    def to_dict(self) -> dict:
        return {
            "arrows": dataclasses.asdict(self.arrows),
            "scene": self.scene.model_dump(),
        }
