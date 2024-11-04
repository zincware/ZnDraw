import base64
import dataclasses

import matplotlib.pyplot as plt
import znjson


@dataclasses.dataclass(frozen=True)
class Figure:
    """Visualize a file or a matplotlib figure."""

    path: str | None = None
    figure: plt.Figure | None = None

    def __post_init__(self):
        if self.path is not None and self.figure is not None:
            raise ValueError("Figure can't have both path and figure")

    def to_base64(self) -> str:
        if self.path is not None:
            with open(self.path, "rb") as image_file:
                return base64.b64encode(image_file.read()).decode("utf-8")
        else:
            import io

            buf = io.BytesIO()
            self.figure.savefig(buf, format="png")
            buf.seek(0)
            return base64.b64encode(buf.read()).decode("utf-8")


class FigureConverter(znjson.ConverterBase):
    level = 100
    instance = Figure
    representation = "zndraw.Figure"

    def encode(self, obj: Figure) -> dict:
        return {"path": obj.path, "base64": obj.to_base64()}

    def decode(self, value: dict) -> Figure:
        return Figure(value["path"])
