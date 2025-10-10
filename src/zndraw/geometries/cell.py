from pydantic import Field
import typing as t
from .base import BaseGeometry, DataProp


class Cell(BaseGeometry):
    """A geometry to draw the simulation cell."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        schema["properties"]["position"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["position"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["position"]["type"] = "string"
        schema["properties"]["position"].pop("anyOf", None)

        schema["properties"]["color"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["color"]["x-features"] = ["color-picker"]
        # TODO: add dropdown with "default" option, where default means follow color scheme
        schema["properties"]["color"]["type"] = "string"
        schema["properties"]["color"].pop("anyOf", None)
        
        return schema
    
    color: str = Field(default="default", description="Color of the cell")

    position: str = Field(
        default="cell",
        description="Cell vectors. Should be a 3x3 matrix.",
    )

    thickness: float = Field(
        default=2.0,
        ge=0.0,
        description="Line thickness for the cell vectors.",
    )

    material: t.Literal[
        "LineBasicMaterial",
        "LineDashedMaterial",
    ] = Field(
        default="LineBasicMaterial",
        description="Material type (static config, not fetched from server)",
    )
