import enum
import typing as t
from pydantic import BaseModel

class ExtensionType(str, enum.Enum):
    MODIFIER = "modifiers"
    SELECTION = "selections"
    ANALYSIS = "analyses"

class Extension(BaseModel):
    """The base class for all ZnDraw extensions."""
    # This is a class attribute, not an instance attribute.
    category: t.ClassVar[ExtensionType]