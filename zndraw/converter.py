import ase
import numpy as np
import znjson
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms
from ase.data.colors import jmol_colors

from zndraw.draw import Object3D
from zndraw.type_defs import ASEDict
from zndraw.utils import get_scaled_radii, rgb2hex


class ASEConverter(znjson.ConverterBase):
    """Encode/Decode datetime objects

    Attributes
    ----------
    level: int
        Priority of this converter over others.
        A higher level will be used first, if there
        are multiple converters available
    representation: str
        An unique identifier for this converter.
    instance:
        Used to select the correct converter.
        This should fulfill isinstance(other, self.instance)
        or __eq__ should be overwritten.
    """

    level = 100
    representation = "ase.Atoms"
    instance = ase.Atoms

    def encode(self, obj: ase.Atoms) -> ASEDict:
        """Convert the datetime object to str / isoformat"""

        numbers = obj.numbers.tolist()
        positions = obj.positions.tolist()
        pbc = obj.pbc.tolist()
        cell = obj.cell.tolist()

        info = {
            k: v
            for k, v in obj.info.items()
            if isinstance(v, (float, int, str, bool, list))
        }
        info |= {k: v.tolist() for k, v in obj.info.items() if isinstance(v, np.ndarray)}
        vectors = info.pop("vectors", [])
        if isinstance(vectors, np.ndarray):
            vectors = vectors.tolist()
        for idx, vector in enumerate(vectors):
            if isinstance(vector, np.ndarray):
                vectors[idx] = vector.tolist()

        if len(vectors) != 0:
            vectors = np.array(vectors)
            if vectors.ndim != 3:
                raise ValueError(
                    f"Vectors must be of shape (n, 2, 3), found '{vectors.shape}'"
                )
            if vectors.shape[1] != 2:
                raise ValueError(
                    f"Vectors must be of shape (n, 2, 3), found '{vectors.shape}'"
                )
            if vectors.shape[2] != 3:
                raise ValueError(
                    f"Vectors must be of shape (n, 2, 3), found '{vectors.shape}'"
                )

            vectors = vectors.tolist()

        if obj.calc is not None:
            calc = {
                k: v
                for k, v in obj.calc.results.items()
                if isinstance(v, (float, int, str, bool, list))
            }
            calc |= {
                k: v.tolist()
                for k, v in obj.calc.results.items()
                if isinstance(v, np.ndarray)
            }
        else:
            calc = {}

        # All additional information should be stored in calc.results
        # and not in calc.arrays, thus we will not convert it here!
        arrays = {}
        if ("colors" not in obj.arrays) or ("" in obj.arrays["colors"]):
            arrays["colors"] = [rgb2hex(jmol_colors[number]) for number in numbers]
        else:
            arrays["colors"] = (
                obj.arrays["colors"].tolist()
                if isinstance(obj.arrays["colors"], np.ndarray)
                else obj.arrays["colors"]
            )

        if ("radii" not in obj.arrays) or (0 in obj.arrays["radii"]):
            # arrays["radii"] = [covalent_radii[number] for number in numbers]
            arrays["radii"] = [get_scaled_radii()[number] for number in numbers]
        else:
            arrays["radii"] = (
                obj.arrays["radii"].tolist()
                if isinstance(obj.arrays["radii"], np.ndarray)
                else obj.arrays["radii"]
            )

        for key in obj.arrays:
            if key in ["colors", "radii"]:
                continue
            if isinstance(obj.arrays[key], np.ndarray):
                arrays[key] = obj.arrays[key].tolist()
            else:
                arrays[key] = obj.arrays[key]

        if hasattr(obj, "connectivity") and obj.connectivity is not None:
            connectivity = (
                obj.connectivity.tolist()
                if isinstance(obj.connectivity, np.ndarray)
                else obj.connectivity
            )
        else:
            connectivity = []

        constraints = []
        if len(obj.constraints) > 0:
            for constraint in obj.constraints:
                if isinstance(constraint, FixAtoms):
                    constraints.append(
                        {"type": "FixAtoms", "indices": constraint.index.tolist()}
                    )
                else:
                    # Can not serialize other constraints
                    pass

        # We don't want to send positions twice
        arrays.pop("positions", None)
        arrays.pop("numbers", None)

        return ASEDict(
            numbers=numbers,
            positions=positions,
            connectivity=connectivity,
            arrays=arrays,
            info=info,
            calc=calc,
            pbc=pbc,
            cell=cell,
            vectors=vectors,
            constraints=constraints,
        )

    def decode(self, value: ASEDict) -> ase.Atoms:
        """Create datetime object from str / isoformat"""
        atoms = ase.Atoms(
            numbers=value["numbers"],
            positions=value["positions"],
            info=value["info"],
            pbc=value["pbc"],
            cell=value["cell"],
        )
        if connectivity := value.get("connectivity"):
            # or do we want this to be nx.Graph?
            atoms.connectivity = np.array(connectivity)
        for key, val in value["arrays"].items():
            atoms.arrays[key] = np.array(val)
        if calc := value.get("calc"):
            atoms.calc = SinglePointCalculator(atoms)
            atoms.calc.results.update(calc)
        if vectors := value.get("vectors"):
            atoms.info["vectors"] = vectors
        if constraints := value.get("constraints"):
            for constraint in constraints:
                if constraint["type"] == "FixAtoms":
                    atoms.set_constraint(FixAtoms(constraint["indices"]))
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
