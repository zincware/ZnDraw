import typing as t
from collections.abc import MutableMapping

from pydantic import BaseModel

from zndraw.geometries import Arrow, Sphere

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Geometries(MutableMapping):
    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def __repr__(self) -> str:
        return f"Scene({self.vis._geometries!r})"

    def __str__(self) -> str:
        return f"Scene({self.vis._geometries!r})"

    def _geometry_to_model(self, geometry_type, geometry_data) -> BaseModel:
        if geometry_type == "Sphere":
            return Sphere(**geometry_data)
        elif geometry_type == "Arrow":
            return Arrow(**geometry_data)
        else:
            raise ValueError(f"Unknown geometry type: {geometry_type}")

    def __getitem__(self, key: str) -> BaseModel:
        if key not in self.vis._geometries:
            response = self.vis.api.get_geometries()
            if response is None:
                raise KeyError(f"Geometry with key '{key}' does not exist")
            self.vis._geometries = response
        if key not in self.vis._geometries:
            raise KeyError(f"Geometry with key '{key}' does not exist")
        return self._geometry_to_model(
            geometry_type=self.vis._geometries[key]["type"],
            geometry_data=self.vis._geometries[key]["data"],
        )

    def __setitem__(self, key: str, value: BaseModel) -> None:
        from zndraw.geometries import geometries

        geometry_type = type(value).__name__
        if geometry_type not in geometries:
            raise ValueError(f"Unknown geometry type: {geometry_type}")
        self.vis._geometries[key] = {"type": geometry_type, "data": value.model_dump()}
        self.vis.api.set_geometry(
            data=value.model_dump(), key=key, geometry_type=geometry_type
        )

    def __delitem__(self, key: str) -> None:
        if key not in self.vis._geometries:
            raise KeyError(f"Geometry with key '{key}' does not exist")
        del self.vis._geometries[key]
        self.vis.api.del_geometry(key=key)

    def __iter__(self):
        return iter(self.vis._geometries)

    def __len__(self) -> int:
        return len(self.vis._geometries)
