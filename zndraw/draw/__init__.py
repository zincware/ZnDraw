import typing as t

from pydantic import BaseModel, ConfigDict, Field

from zndraw.base import Extension, MethodsCollection


def _update_material_schema(schema: dict) -> dict:
    schema["properties"]["wireframe"]["format"] = "checkbox"
    schema["properties"]["color"]["format"] = "color"
    schema["properties"]["opacity"]["format"] = "range"
    schema["properties"]["opacity"]["step"] = 0.01
    schema["properties"]["outlines"]["format"] = "checkbox"
    return schema


def _update_object3d_schema(schema: dict) -> dict:
    """Remote position, rotation, and scale from the schema."""
    schema["properties"].pop("position", None)
    schema["properties"].pop("rotation", None)
    schema["properties"].pop("scale", None)
    return schema


class Material(BaseModel):
    """Base class for the Object3D material."""

    # TODO, reuse / combine with scene materials
    color: str = "#62929E"
    opacity: float = Field(0.2, ge=0.0, le=1.0)
    wireframe: bool = False
    outlines: bool = False

    model_config = ConfigDict(json_schema_extra=_update_material_schema)


class Object3D(Extension):
    """Base class for all 3D objects."""

    material: Material = Field(default_factory=Material)

    position: t.Tuple[float, float, float] | list[float] = (0, 0, 0)
    rotation: t.Tuple[float, float, float] | list[float] = (0, 0, 0)
    scale: t.Tuple[float, float, float] | list[float] = (1, 1, 1)

    model_config = ConfigDict(json_schema_extra=_update_object3d_schema)

    def run(self, vis, **kwargs) -> None:
        # get the selected particles and compute the COM
        if len(vis.selection) > 0:
            selected = vis.atoms[vis.selection]
            self.position = selected.get_center_of_mass().tolist()
        # self.position = vis.points[0].tolist()
        print(f"Running {self.__class__.__name__} at {self.position}")

        vis.geometries = vis.geometries + [self]


class Plane(Object3D):
    width: float = 10.0
    height: float = 10.0


class Box(Object3D):
    width: float = 10.0
    height: float = 10.0
    depth: float = 10.0


class Circle(Object3D):
    radius: float = 5.0


class Cone(Object3D):
    radius: float = 5.0
    height: float = 10.0


class Cylinder(Object3D):
    radius_top: float = 5.0
    radius_bottom: float = 5.0
    height: float = 10.0


class Dodecahedron(Object3D):
    radius: float = 5.0


class Icosahedron(Object3D):
    radius: float = 5.0


class Octahedron(Object3D):
    radius: float = 5.0


class Ring(Object3D):
    inner_radius: float = 1
    outer_radius: float = 4.0


class Sphere(Object3D):
    radius: float = 4.0


class Tetrahedron(Object3D):
    radius: float = 5.0


class Torus(Object3D):
    radius: float = 3.0
    tube: float = 1.0


class TorusKnot(Object3D):
    radius: float = 3.0
    tube: float = 1.0


methods = t.Union[
    Plane,
    Sphere,
    Box,
    Circle,
    Cone,
    Cylinder,
    Dodecahedron,
    Icosahedron,
    Octahedron,
    Ring,
    Tetrahedron,
    Torus,
    TorusKnot,
]


class Geometry(MethodsCollection):
    method: methods = Field(default_factory=Box, description="Select a geometry method.")

    # @classmethod
    # def updated_schema(cls) -> dict:
    #     schema = super().updated_schema()

    #     schema["properties"]["wireframe"]["format"] = "checkbox"
    #     schema["properties"]["color"]["format"] = "color"
    #     schema["properties"]["opacity"]["format"] = "range"
    #     schema["properties"]["opacity"]["step"] = 0.01

    #     return schema
