"""Property inspector configuration as a scene object."""

from pydantic import BaseModel, ConfigDict, Field


class PropertyInspector(BaseModel):
    """Per-particle property display configuration.

    Controls which atom/particle properties are shown in the inspector panel.
    """

    model_config = ConfigDict(frozen=True)

    owner: str | None = Field(
        default=None,
        description="User ID of the property inspector owner.",
        json_schema_extra={"x-custom-type": "ownership-toggle"},
    )
    active: bool = Field(
        default=False,
        description="Show property inspector panel.",
    )
    enabled_properties: list[str] = Field(
        default_factory=list,
        description="Selected property keys to display in the inspector table.",
        json_schema_extra={"x-custom-type": "property-inspector"},
    )
