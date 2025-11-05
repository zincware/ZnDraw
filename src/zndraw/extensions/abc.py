import enum
import typing as t

from pydantic import BaseModel

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw


class Category(str, enum.Enum):
    MODIFIER = "modifiers"
    SELECTION = "selections"
    ANALYSIS = "analysis"


class Extension(BaseModel):
    """The base class for all ZnDraw extensions."""

    # This is a class attribute, not an instance attribute.
    category: t.ClassVar[Category]

    def run(self, vis: "ZnDraw", **kwargs):
        """Run the extension.

        This method should be overridden by subclasses to implement the extension's functionality.

        Args:
            vis (ZnDraw): The ZnDraw instance.
            **kwargs: Additional keyword arguments specific to the extension.
        """
        raise NotImplementedError("Subclasses must implement the run method.")
