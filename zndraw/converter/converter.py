import ase
import numpy as np
import znjson
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms

from zndraw.draw import Object3D
from zndraw.figure import Figure
from zndraw.type_defs import ASEDict
import znsocket
import plotly.graph_objects as go
import plotly.io as pio

# TODO: there is an issue with using `_type` with znjson and the numpy array type,
TYPE_KEY = "type"

class ASEConverter(znjson.ConverterBase):
    level = 100
    representation = "ase.Atoms"
    instance = ase.Atoms

    def encode(self, obj: ase.Atoms) -> ASEDict:
        def recursive_encode(val):
            if isinstance(val, (str, int, float, bool, type(None))):
                return val
            elif isinstance(val, np.ndarray):
                return {TYPE_KEY: "ndarray", "value": val.tolist()}
            elif isinstance(val, np.generic):
                return val.item()  # Convert numpy generic types to native Python types
            elif isinstance(val, Figure):
                return {TYPE_KEY: "zndraw.Figure", "base64": val.base64}
            elif isinstance(val, (list, tuple)):
                return [recursive_encode(v) for v in val]
            elif isinstance(val, dict):
                return {k: recursive_encode(v) for k, v in val.items()}
            elif isinstance(val, go.Figure):
                return {TYPE_KEY: "plotly.graph_objs.Figure", "value": val.to_json()}
            # elif isinstance(val, znsocket.Dict):
            #     return {k: recursive_encode(v) for k, v in val.items()}
            # elif isinstance(val, (znsocket.List, znsocket.Segments)):
            #     return [recursive_encode(v) for v in val]
            else:
                raise TypeError(f"Unsupported type during encoding: {type(val)} / {val}")

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
                    {TYPE_KEY: "FixAtoms", "indices": constraint.index.tolist()}
                )
            else:
                raise TypeError(f"Unsupported constraint type: {type(constraint)}")

        return data

    def decode(self, value: ASEDict) -> ase.Atoms:
        def recursive_decode(val):
            if isinstance(val, (str, int, float, bool, type(None))):
                return val
            elif isinstance(val, list):
                return [recursive_decode(v) for v in val]
            elif isinstance(val, dict):
                if TYPE_KEY in val:
                    if val[TYPE_KEY] == "zndraw.Figure":
                        return Figure(base64=val["base64"])
                    elif val[TYPE_KEY] == "ndarray":
                        return np.array(val["value"])
                    elif val[TYPE_KEY] == "plotly.graph_objs.Figure":
                        return pio.from_json(val["value"])
                    else:
                        raise TypeError(f"Unsupported type during decoding: {val[TYPE_KEY]}")
                else:
                    return {k: recursive_decode(v) for k, v in val.items()}
            # special znsocket cases, they need to be resolved
            elif isinstance(val, znsocket.Dict):
                return recursive_decode(dict(val))
            elif isinstance(val, (znsocket.List, znsocket.Segments)):
                return recursive_decode(val[:])
            else:
                raise TypeError(f"Unsupported type during decoding: {type(val)} / {val}")

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
                if constraint[TYPE_KEY] == "FixAtoms":
                    atoms.set_constraint(FixAtoms(constraint["indices"]))
                else:
                    raise TypeError(f"Unsupported constraint type: {constraint[TYPE_KEY]}")

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
