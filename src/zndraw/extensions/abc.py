"""Base classes for ZnDraw extensions."""

import enum
import typing as t

from pydantic import BaseModel


class Category(str, enum.Enum):
    """Extension category types."""

    MODIFIER = "modifiers"
    SELECTION = "selections"
    ANALYSIS = "analysis"


class Extension(BaseModel):
    """Base class for all ZnDraw extensions.

    Extensions are Pydantic models that define their parameters as fields.
    The JSON schema is generated from the model for the frontend forms.
    """

    # Class attribute (not instance attribute)
    category: t.ClassVar[Category]

    def run(self, vis: t.Any, **kwargs: t.Any) -> t.Any:
        """Run the extension.

        Parameters
        ----------
        vis
            The ZnDraw instance.
        **kwargs
            Additional keyword arguments specific to the extension.

        Returns
        -------
        Any
            Extension-specific return value.

        Raises
        ------
        NotImplementedError
            If the subclass does not implement this method.
        """
        raise NotImplementedError("Subclasses must implement the run method.")
