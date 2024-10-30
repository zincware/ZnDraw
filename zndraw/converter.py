# TODO: move ASEConverter to here
import znjson

from zndraw.draw import Object3D


class Object3DConverter(znjson.ConverterBase):
    instance: type = Object3D
    representation: str = "zndraw.Object3D"
    level: int = 100

    def encode(self, obj: Object3D) -> dict:
        return {"class": obj.__class__.__name__, "data": obj.model_dump()}

    def decode(self, value: str) -> Object3D:
        import zndraw

        cls = getattr(zndraw, value["class"])

        return cls(**value["data"])
