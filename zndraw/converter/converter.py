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
        data = {
            "numbers": obj.numbers.tolist(),
            "positions": obj.positions.tolist(),
            "pbc": obj.pbc.tolist(),
            "cell": obj.cell.tolist(),
        }
        info = {}
        for key, value in obj.info.items():
            if isinstance(value, (str, int, float, bool)):
                info[key] = value
            elif isinstance(value, np.ndarray):
                info[key] = value.tolist()
            elif isinstance(value, Figure):
                info[key] = {"_type": "zndraw.Figure", "base64": value.base64}
            else:
                raise TypeError(f"Unsupported type in info: {type(value)}")
        data["info"] = info

        arrays = {}
        for key, value in obj.arrays.items():
            if isinstance(value, np.ndarray):
                arrays[key] = value.tolist()
            elif isinstance(value, Figure):
                arrays[key] = {"_type": "zndraw.Figure", "base64": value.base64}
            else:
                raise TypeError(f"Unsupported type in arrays: {type(value)}")
        data["arrays"] = arrays

        if obj.calc is not None:
            calc_results = {}
            for key, value in obj.calc.results.items():
                if isinstance(value, (str, int, float, bool)):
                    calc_results[key] = value
                elif isinstance(value, np.ndarray):
                    calc_results[key] = value.tolist()
                else:
                    raise TypeError(f"Unsupported type in calc results: {type(value)}")
            data["calc"] = calc_results
        else:
            data["calc"] = {}

        constraints = []
        for constraint in obj.constraints:
            if isinstance(constraint, FixAtoms):
                constraints.append(
                    {"_type": "FixAtoms", "indices": constraint.index.tolist()}
                )
            else:
                raise TypeError(
                    f"Unsupported constraint type: {type(constraint)}"
                )
        data["constraints"] = constraints

        return data

    def decode(self, value: dict) -> ase.Atoms:
        atoms = ase.Atoms(
            numbers=value["numbers"],
            positions=value["positions"],
            pbc=value["pbc"],
            cell=value["cell"],
        )
        for key, val in value["info"].items():
            if isinstance(val, dict):
                if val.keys() == {"_type", "base64"} and val["_type"] == "zndraw.Figure":
                    atoms.info[key] = Figure(base64=val["base64"])
                else:
                    atoms.info[key] = val
            else:
                atoms.info[key] = val
        for key, val in value["arrays"].items():
            atoms.arrays[key] = np.array(val)
        if "calc" in value:
            atoms.calc = SinglePointCalculator(atoms)
            for key, val in value["calc"].items():
                if isinstance(val, list):
                    atoms.calc.results[key] = np.array(val)
                else:
                    atoms.calc.results[key] = val
        if "constraints" in value:
            for constraint in value["constraints"]:
                if constraint["_type"] == "FixAtoms":
                    atoms.set_constraint(FixAtoms(constraint["indices"]))
                else:
                    raise TypeError(
                        f"Unsupported constraint type: {constraint['_type']}"
                    )

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
