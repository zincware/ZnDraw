"""Arrow geometry for ZnDraw."""

import typing as t
from pydantic import Field

from .base import BaseGeometry, DataProp


class Arrow(BaseGeometry):
    """Arrow geometry with direction vector.

    Arrows are defined by a start position and a direction vector.
    The direction vector determines both the arrow's orientation and its length.
    """

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        
        # Position field: dropdown of atom props only
        schema["properties"]["position"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["position"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["position"]["type"] = "string"
        schema["properties"]["position"].pop("anyOf", None)

        # Direction field: dropdown of atom props (vectors)
        schema["properties"]["direction"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["direction"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["direction"]["type"] = "string"
        schema["properties"]["direction"].pop("anyOf", None)
        
        # Color field: dropdown + free text + color picker
        schema["properties"]["color"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["color"]["x-features"] = [
            "color-picker",
            "dynamic-atom-props",
            "free-solo",
        ]
        schema["properties"]["color"]["type"] = "string"
        schema["properties"]["color"].pop("anyOf", None)
        
        return schema

    # Use 'start' as the primary field with 'position' alias for compatibility
    position: DataProp = Field(
        default="arrays.positions",
        description="Arrow start position [x,y,z]. String for dynamic data key, tuple for static value.",
    )

    direction: DataProp = Field(
        default="calc.forces",
        description="Direction vector [x,y,z]. Defines arrow orientation and base length. String for dynamic data key, tuple for static value.",
    )

    radius: float = Field(
        default=1,
        description="Arrow shaft radius. String for dynamic data key, float for static value.",
    )

    scale: float = Field(
        default=1.0,
        description="Length scale multiplier applied to direction vector length. String for dynamic data key, float for static value.",
    )
