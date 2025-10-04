from pydantic import BaseModel, ConfigDict, Field


class Sphere(BaseModel):
    """A sphere geometry."""

    model_config = ConfigDict(frozen=True)
    position: str | tuple[float, float, float] | list[tuple[float, float, float]] = (
        Field(
            default="arrays.positions",
            description="The position of the sphere. Dynamic via e.g. `arrays.positions`, or static for one or more spheres via `(x, y, z)`.",
            # examples=[
            #     "0,0,0",
            #     (0.0, 0.0, 0.0),
            #     [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)],
            # ],
        )
    )
    radius: str | float = Field(
        default="arrays.radii",
        description="The radius of the sphere. Dynamic via e.g. `arrays.radii`, or static via a float.",
    )
    color: str | tuple[float, float, float] | list[tuple[float, float, float]] = Field(
        default="arrays.colors",
        description="The color of the sphere. Dynamic via e.g. `arrays.colors`, or static via an RGB tuple or list of RGB tuples.",
    )
