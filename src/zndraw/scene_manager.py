import typing as t
from collections.abc import MutableMapping

from pydantic import BaseModel

from zndraw.geometries import (
    Arrow,
    Bond,
    Box,
    Camera,
    Cell,
    Curve,
    Floor,
    Plane,
    Sphere,
)

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Geometries(MutableMapping):
    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def __repr__(self) -> str:
        keys = list(self.vis._geometries.keys())
        return f"Geometries(keys={keys})"

    def __str__(self) -> str:
        keys = list(self.vis._geometries.keys())
        return f"Geometries(keys={keys})"

    def _geometry_to_model(self, geometry_type, geometry_data) -> BaseModel:
        if geometry_type == "Sphere":
            return Sphere(**geometry_data)
        elif geometry_type == "Arrow":
            return Arrow(**geometry_data)
        elif geometry_type == "Bond":
            return Bond(**geometry_data)
        elif geometry_type == "Curve":
            return Curve(**geometry_data)
        elif geometry_type == "Camera":
            return Camera(**geometry_data)
        elif geometry_type == "Cell":
            return Cell(**geometry_data)
        elif geometry_type == "Floor":
            return Floor(**geometry_data)
        elif geometry_type == "Box":
            return Box(**geometry_data)
        elif geometry_type == "Plane":
            return Plane(**geometry_data)
        else:
            raise ValueError(f"Unknown geometry type: {geometry_type}")

    def __getitem__(self, key: str) -> BaseModel:
        """Get a geometry by key.

        First checks local cache, then fetches from server if needed.
        """
        if key not in self.vis._geometries:
            # Try to fetch this specific geometry from server
            response = self.vis.api.get_geometry(key)
            if response is None:
                raise KeyError(f"Geometry with key '{key}' does not exist")
            # Response has correct structure with 'type' and 'data' keys
            from zndraw.zndraw import _GeometryStore

            self.vis._geometries[key] = t.cast(_GeometryStore, response)

        geometry_data = self.vis._geometries[key]
        return self._geometry_to_model(
            geometry_type=geometry_data["type"],
            geometry_data=geometry_data["data"],
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
        """Iterate over geometry keys.

        If local cache is empty, fetches keys from server first.
        """
        if not self.vis._geometries:
            # Fetch keys from server
            keys = self.vis.api.list_geometries()
            # Populate cache with None to track keys exist
            for key in keys:
                if key not in self.vis._geometries:
                    # Mark as known but not loaded
                    pass
            return iter(keys)
        return iter(self.vis._geometries)

    def __len__(self) -> int:
        """Return the number of geometries.

        If local cache is empty, fetches keys from server first.
        """
        if not self.vis._geometries:
            keys = self.vis.api.list_geometries()
            return len(keys)
        return len(self.vis._geometries)
