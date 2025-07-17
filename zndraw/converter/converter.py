import ase
import numpy as np
import znjson
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms

from zndraw.draw import Object3D
from zndraw.figure import Figure
from zndraw.type_defs import ASEDict

class ASEConverter(znjson.ConverterBase):
    level = 100
    representation = "ase.Atoms"
    instance = ase.Atoms

    def encode(self, obj: ase.Atoms) -> dict:
        def recursive_encode(val):
            if isinstance(val, (str, int, float, bool, type(None))):
                return val
            elif isinstance(val, np.ndarray):
                return {"type": "ndarray", "value": val.tolist()}
            elif isinstance(val, Figure):
                return {"type": "zndraw.Figure", "base64": val.base64}
            elif isinstance(val, list):
                return [recursive_encode(v) for v in val]
            elif isinstance(val, tuple):
                return {"type": "tuple", "value": [recursive_encode(v) for v in val]}
            elif isinstance(val, dict):
                return {k: recursive_encode(v) for k, v in val.items()}
            else:
                raise TypeError(f"Unsupported type during encoding: {type(val)}")

        data = {
            "numbers": obj.numbers.tolist(),
            "positions": obj.positions.tolist(),
            "pbc": obj.pbc.tolist(),
            "cell": obj.cell.tolist(),
            "info": recursive_encode(obj.info),
            "arrays": {k: recursive_encode(v) for k, v in obj.arrays.items()},
            "constraints": [],
            "calc": {},
        }

        if obj.calc is not None:
            data["calc"] = recursive_encode(obj.calc.results)

        for constraint in obj.constraints:
            if isinstance(constraint, FixAtoms):
                data["constraints"].append(
                    {"type": "FixAtoms", "indices": constraint.index.tolist()}
                )
            else:
                raise TypeError(f"Unsupported constraint type: {type(constraint)}")

        return data

    def decode(self, value: dict) -> ase.Atoms:
        def recursive_decode(val):
            if isinstance(val, (str, int, float, bool, type(None))):
                return val
            elif isinstance(val, list):
                return [recursive_decode(v) for v in val]
            elif isinstance(val, dict):
                if "type" in val:
                    if val["type"] == "zndraw.Figure":
                        return Figure(base64=val["base64"])
                    elif val["type"] == "ndarray":
                        return np.array(val["value"])
                    elif val["type"] == "tuple":
                        return tuple(recursive_decode(v) for v in val["value"])
                    else:
                        raise TypeError(f"Unsupported type during decoding: {val['type']}")
                else:
                    return {k: recursive_decode(v) for k, v in val.items()}
            else:
                raise TypeError(f"Unsupported type during decoding: {type(val)}")

        atoms = ase.Atoms(
            numbers=value["numbers"],
            positions=value["positions"],
            pbc=value["pbc"],
            cell=value["cell"],
        )
        atoms.info = recursive_decode(value.get("info", {}))

        for key, val in value.get("arrays", {}).items():
            atoms.arrays[key] = recursive_decode(val)

        if "calc" in value and value["calc"]:
            atoms.calc = SinglePointCalculator(atoms)
            atoms.calc.results = recursive_decode(value["calc"])

        if "constraints" in value:
            for constraint in value["constraints"]:
                if constraint["type"] == "FixAtoms":
                    atoms.set_constraint(FixAtoms(constraint["indices"]))
                else:
                    raise TypeError(f"Unsupported constraint type: {constraint['type']}")

        return atoms



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
