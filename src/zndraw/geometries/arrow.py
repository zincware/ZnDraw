from pydantic import BaseModel, ConfigDict, Field


class Arrow(BaseModel):
    """A arrow geometry."""

    model_config = ConfigDict(frozen=True)
    start: str | tuple[float, float, float] | list[tuple[float, float, float]] = (
        Field(
            default="arrays.positions",
            description="The start position of the arrow. Dynamic via e.g. `arrays.positions`, or static for one or more arrows via `(x, y, z)`.",
            # examples=[
            #     "0,0,0",
            #     (0.0, 0.0, 0.0),
            #     [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)],
            # ],
        )
    )
    direction: str | tuple[float, float, float] | list[tuple[float, float, float]] = (
        Field(
            default="calc.forces",
            description="The direction of the arrow. Dynamic via e.g. `arrays.directions`, or static for one or more arrows via `(x, y, z)`.",
        )
    )
    radius: str | float | int = Field(
        default=1,
        description="The radius of the arrow. Dynamic via e.g. `arrays.radii`, or static via a float.",
    )
    color: str | tuple[float, float, float] | list[tuple[float, float, float]] = Field(
        default="arrays.colors",
        description="The color of the arrow. Dynamic via e.g. `arrays.colors`, or static via an RGB tuple or list of RGB tuples.",
    )
    scale: str | float | int = Field(
        default=1,
        description="A global scale factor for the arrow size. Dynamic via e.g. `arrays.scales`, or static via a float.",
    )
