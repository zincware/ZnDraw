import typing as t

from pydantic import BaseModel, ConfigDict, Field


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
    color: str = "#62929E"
    opacity: float = Field(default=0.2, ge=0.0, le=1.0)
    wireframe: bool = False
    outlines: bool = False

    model_config = ConfigDict(json_schema_extra=_update_material_schema)


class Object3D(BaseModel):
    """Base class for all 3D objects."""

    material: Material = Material()

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

        vis.geometries.append(self)  # TODO: dump / load without pickle


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


class Rhomboid(Object3D):
    vectorA: t.Tuple[float, float, float] | list[float] = (10, 0, 0)
    vectorB: t.Tuple[float, float, float] | list[float] = (0, 10, 0)
    vectorC: t.Tuple[float, float, float] | list[float] = (0, 0, 10)


class Custom2DShape(Object3D):
    points: list[tuple[float, float]]


class Ellipsoid(Object3D):
    a: float = 10.0
    b: float = 5.0
    c: float = 5.0


geometries: dict[str, t.Type[Object3D]] = {
    Plane.__name__: Plane,
    Box.__name__: Box,
    Circle.__name__: Circle,
    Cone.__name__: Cone,
    Cylinder.__name__: Cylinder,
    Dodecahedron.__name__: Dodecahedron,
    Icosahedron.__name__: Icosahedron,
    Octahedron.__name__: Octahedron,
    Ring.__name__: Ring,
    Sphere.__name__: Sphere,
    Tetrahedron.__name__: Tetrahedron,
    Torus.__name__: Torus,
    TorusKnot.__name__: TorusKnot,
    Rhomboid.__name__: Rhomboid,
    Ellipsoid.__name__: Ellipsoid,
}
