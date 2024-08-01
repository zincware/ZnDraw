import dataclasses
import typing as t

from zndraw.utils import emit_with_retry

if t.TYPE_CHECKING:
    from zndraw import ZnDraw

HSLColor = t.Tuple[float, float, float]


@dataclasses.dataclass
class ArrowsConfig:
    colormap: list[HSLColor]
    normalize: bool
    colorrange: tuple[float, float]
    scale_vector_thickness: bool
    opacity: float

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


def _create_arrows_config() -> ArrowsConfig:
    return ArrowsConfig(
        colormap=[[-0.5, 0.9, 0.5], [0.0, 0.9, 0.5]],
        normalize=True,
        colorrange=[0, 1],
        scale_vector_thickness=False,
        opacity=1.0,
    )


@dataclasses.dataclass
class ZnDrawConfig:
    vis: "ZnDraw" = dataclasses.field(repr=False)
    arrows: ArrowsConfig = dataclasses.field(default_factory=_create_arrows_config)

    def __post_init__(self) -> None:
        self.arrows.set_vis(self.vis)
