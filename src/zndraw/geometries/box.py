"""Box geometry for ZnDraw."""

import typing as t
from pydantic import Field

from .base import BaseGeometry, DataProp, InteractionSettings


class Box(BaseGeometry):
    """Box/cuboid geometry with customizable dimensions.

    Boxes are defined by a center position and size dimensions [width, height, depth].
    Optional rotation can be applied as Euler angles [x, y, z] in radians.
    """

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Position field: dropdown of atom props
        schema["properties"]["position"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["position"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["position"]["type"] = "string"
        schema["properties"]["position"].pop("anyOf", None)

        # Size field: dropdown of atom props
        schema["properties"]["size"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["size"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["size"]["type"] = "string"
        schema["properties"]["size"].pop("anyOf", None)

        # Color field: dropdown + free text + color picker
        schema["properties"]["color"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["color"]["x-features"] = [
            "color-picker",
            "dynamic-atom-props",
            "free-solo",
        ]
        schema["properties"]["color"]["type"] = "string"
        schema["properties"]["color"].pop("anyOf", None)

        # Rotation field: dropdown of atom props
        schema["properties"]["rotation"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["rotation"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["rotation"]["type"] = "string"
        schema["properties"]["rotation"].pop("anyOf", None)

        # Selection/hover color pickers
        schema["$defs"]["InteractionSettings"]["properties"]["color"][
            "x-custom-type"
        ] = "dynamic-enum"
        schema["$defs"]["InteractionSettings"]["properties"]["color"][
            "x-features"
        ] = ["color-picker", "dynamic-atom-props", "free-solo"]

        return schema

    size: DataProp = Field(
        default=(1.0, 1.0, 1.0),
        description="Box dimensions [width, height, depth]. String for dynamic data key, tuple for static value.",
    )

    rotation: DataProp = Field(
        default=(0.0, 0.0, 0.0),
        description="Rotation as Euler angles [x, y, z] in radians. String for dynamic data key, tuple for static value.",
    )

    scale: float = Field(
        default=1.0,
        ge=0.0,
        description="Uniform scale factor applied to box size.",
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Box opacity, between 0 (transparent) and 1 (opaque).",
    )

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )

    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )
