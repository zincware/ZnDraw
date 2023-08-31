import typing as t

from pydantic import BaseModel, Field


class PlaneGeometry(BaseModel):
    method: t.Literal["PlaneGeometry"] = Field("PlaneGeometry")
    width: float = 10.0
    height: float = 10.0


class BoxGeometry(BaseModel):
    method: t.Literal["BoxGeometry"] = Field("BoxGeometry")
    width: float = 10.0
    height: float = 10.0
    depth: float = 10.0


class CircleGeometry(BaseModel):
    method: t.Literal["CircleGeometry"] = Field("CircleGeometry")
    radius: float = 5.0


class ConeGeometry(BaseModel):
    method: t.Literal["ConeGeometry"] = Field("ConeGeometry")
    radius: float = 5.0
    height: float = 10.0


class CylinderGeometry(BaseModel):
    method: t.Literal["CylinderGeometry"] = Field("CylinderGeometry")
    radius_top: float = 5.0
    radius_bottom: float = 5.0
    height: float = 10.0


class DodecahedronGeometry(BaseModel):
    method: t.Literal["DodecahedronGeometry"] = Field("DodecahedronGeometry")
    radius: float = 5.0


class IcosahedronGeometry(BaseModel):
    method: t.Literal["IcosahedronGeometry"] = Field("IcosahedronGeometry")
    radius: float = 5.0


class OctahedronGeometry(BaseModel):
    method: t.Literal["OctahedronGeometry"] = Field("OctahedronGeometry")
    radius: float = 5.0


class RingGeometry(BaseModel):
    method: t.Literal["RingGeometry"] = Field("RingGeometry")
    inner_radius: float = 1
    outer_radius: float = 4.0


class SphereGeometry(BaseModel):
    method: t.Literal["SphereGeometry"] = Field("SphereGeometry")
    radius: float = 4.0


class TetrahedronGeometry(BaseModel):
    method: t.Literal["TetrahedronGeometry"] = Field("TetrahedronGeometry")
    radius: float = 5.0


class TorusGeometry(BaseModel):
    method: t.Literal["TorusGeometry"] = Field("TorusGeometry")
    radius: float = 3.0
    tube: float = 1.0


class TorusKnotGeometry(BaseModel):
    method: t.Literal["TorusKnotGeometry"] = Field("TorusKnotGeometry")
    radius: float = 3.0
    tube: float = 1.0


methods = t.Union[
    SphereGeometry,
    PlaneGeometry,
    BoxGeometry,
    CircleGeometry,
    ConeGeometry,
    CylinderGeometry,
    DodecahedronGeometry,
    IcosahedronGeometry,
    OctahedronGeometry,
    RingGeometry,
    TetrahedronGeometry,
    TorusGeometry,
    TorusKnotGeometry,
]


class Geometry(BaseModel):
    geometry: methods = Field(discriminator="method")
    wireframe: bool = True
    color: str = "#62929E"
    opacity: float = Field(0.2, ge=0.0, le=1.0)

    @classmethod
    def updated_schema(cls):
        schema = cls.model_json_schema()
        for prop in [x.__name__ for x in t.get_args(methods)]:
            schema["$defs"][prop]["properties"]["method"]["options"] = {"hidden": True}
            schema["$defs"][prop]["properties"]["method"]["type"] = "string"

        schema["properties"]["wireframe"]["format"] = "checkbox"
        schema["properties"]["color"]["format"] = "color"
        schema["properties"]["opacity"]["format"] = "range"
        schema["properties"]["opacity"]["step"] = 0.01

        return schema
