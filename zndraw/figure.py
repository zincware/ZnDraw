import base64
import dataclasses
import io

import matplotlib.pyplot as plt
import znjson


@dataclasses.dataclass(frozen=True)
class Figure:
    """Visualize a file or a matplotlib figure.
    
    
    Parameters
    ----------
    path : str | None
        Path to the image file.
    figure : plt.Figure | None
        Matplotlib figure object.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from zndraw import ZnDraw, Figure
    >>> vis = ZnDraw(url="http://localhost:8000", token="test_token")
    >>> fig, ax = plt.subplots()
    >>> ax.plot([1, 2, 3], [1, 2, 3])
    >>> vis.figures["mtpl"] = Figure.from_matplotlib(fig)
    >>> vis.figures["png"] = Figure.from_path("path/to/image.png")
    """

    base64: str | None = None

    @classmethod
    def from_path(cls, path: str) -> "Figure":
        """Create a Figure from a file path."""
        with open(path, "rb") as image_file:
            content = base64.b64encode(image_file.read()).decode("utf-8")
        return cls(base64=content)

    @classmethod
    def from_matplotlib(cls, figure: plt.Figure) -> "Figure":
        """Create a Figure from a matplotlib figure."""
        buf = io.BytesIO()
        figure.savefig(buf, format="png")
        buf.seek(0)
        content = base64.b64encode(buf.read()).decode("utf-8")
        return cls(base64=content)


class FigureConverter(znjson.ConverterBase):
    level = 100
    instance = Figure
    representation = "zndraw.Figure"

    def encode(self, obj: Figure) -> dict:
        return {"base64": obj.base64}

    def decode(self, value: dict) -> Figure:
        return Figure(base64=value["base64"])
