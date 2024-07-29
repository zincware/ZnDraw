import functools
import importlib.util
import json
import logging
import pathlib
import socket
import sys
import tempfile
import typing as t
import uuid

import ase
import datamodel_code_generator
import numpy as np
import socketio.exceptions
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from znjson import ConverterBase

log = logging.getLogger(__name__)


class ASEDict(t.TypedDict):
    numbers: list[int]
    positions: list[list[float]]
    connectivity: list[tuple[int, int, int]]
    arrays: dict[str, list[float | int | list[float | int]]]
    info: dict[str, float | int]
    # calc: dict[str, float|int|np.ndarray] # should this be split into arrays and info?
    pbc: list[bool]
    cell: list[list[float]]
    vectors: list[list[list[float]]]


def rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


@functools.lru_cache(maxsize=128)
def get_scaled_radii() -> np.ndarray:
    """Scale down the covalent radii to visualize bonds better."""
    radii = covalent_radii
    # shift the values such that they are in [0.3, 1.3]
    radii = radii - np.min(radii)
    radii = radii / np.max(radii)
    radii = radii + 0.3
    return radii


class ASEConverter(ConverterBase):
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

        if hasattr(obj, "connectivity") and obj.connectivity is not None:
            connectivity = (
                obj.connectivity.tolist()
                if isinstance(obj.connectivity, np.ndarray)
                else obj.connectivity
            )
        else:
            connectivity = []

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
        if "colors" in value["arrays"]:
            atoms.arrays["colors"] = np.array(value["arrays"]["colors"])
        if "radii" in value["arrays"]:
            atoms.arrays["radii"] = np.array(value["arrays"]["radii"])
        if calc := value.get("calc"):
            atoms.calc = SinglePointCalculator(atoms)
            atoms.calc.results.update(calc)
        if vectors := value.get("vectors"):
            atoms.info["vectors"] = vectors
        return atoms


def get_port(default: int) -> int:
    """Get an open port."""
    try:
        sock = socket.socket()
        sock.bind(("", default))
        port = default
    except OSError:
        sock = socket.socket()
        sock.bind(("", 0))
        port = sock.getsockname()[1]
    finally:
        sock.close()
    return port


class ZnDrawLoggingHandler(logging.Handler):
    """Logging handler which emits log messages to the ZnDraw server."""

    def __init__(self, zndraw):
        super().__init__()
        self.zndraw = zndraw

    def emit(self, record):
        try:
            msg = self.format(record)
            self.zndraw.log(msg)
        except RecursionError:  # See StreamHandler
            raise
        except Exception:
            print("Something went wrong")
            self.handleError(record)


def get_cls_from_json_schema(schema: dict, name: str, **kwargs):
    """Get a python class from a json schema."""

    # TODO: needs tests
    # TODO: do not write file but use in-memory

    kwargs["strict_nullable"] = True

    with tempfile.TemporaryDirectory() as temporary_directory_name:
        temporary_directory = pathlib.Path(temporary_directory_name)
        output = temporary_directory / "model.py"
        datamodel_code_generator.generate(
            json.dumps(schema),
            input_file_type=datamodel_code_generator.InputFileType.JsonSchema,
            input_filename="example.json",
            output=output,
            # set up the output model types
            output_model_type=datamodel_code_generator.DataModelType.PydanticV2BaseModel,
            **kwargs,
        )

        ref_module = uuid.uuid4().hex

        spec = importlib.util.spec_from_file_location(ref_module, output)
        module = importlib.util.module_from_spec(spec)
        sys.modules[ref_module] = module
        spec.loader.exec_module(module)

        return getattr(module, name)


def emit_with_retry(
    socket: socketio.Client,
    event,
    data=None,
    namespace=None,
    callback=None,
    retries: int = 1,
    delay: float = 0.1,
    increase_delay: float = 0.1,
    reconnect: bool = False,
) -> None:
    """Emit data to a socket with retries."""
    for idx in range(retries):
        try:
            socket.emit(event=event, data=data, namespace=namespace, callback=callback)
            break
        except socketio.exceptions.BadNamespaceError as err:
            log.error(f"Retrying {event} due to {err}")
            if idx == retries - 1:
                raise err
            if reconnect:
                raise ValueError("Reconnect not implemented yet")
        except Exception as err:
            if idx == retries - 1:
                raise err
        socket.sleep(delay)
        delay += increase_delay


def call_with_retry(
    socket: socketio.Client,
    event,
    data=None,
    namespace=None,
    timeout=60,
    retries: int = 1,
    delay: float = 0.1,
    increase_delay: float = 0.1,
    reconnect: bool = False,
) -> t.Any:
    """Call a function with retries."""
    for idx in range(retries):
        try:
            return socket.call(
                event=event, data=data, namespace=namespace, timeout=timeout
            )
        except (
            socketio.exceptions.TimeoutError,
            socketio.exceptions.BadNamespaceError,
        ) as err:
            log.error(f"Retrying {event} due to {err} ({idx} / {retries})")
            if idx == retries - 1:
                raise err
            if reconnect:
                raise ValueError("Reconnect not implemented yet")
        except Exception as err:
            if idx == retries - 1:
                raise err
        socket.sleep(delay)
        delay += increase_delay
    return None


def direction_to_euler(direction, roll=0):
    """
    Convert a direction vector to euler angles.

    You get an increased degree of freedom by setting the roll angle.
    """
    direction = np.array(direction)
    direction = direction / np.linalg.norm(direction)

    x, y, z = direction

    if x == 0 and y == 0:
        if z > 0:
            yaw = 0
            pitch = 0
        else:
            yaw = 0
            pitch = 180
    else:
        yaw = np.arctan2(y, x)
        pitch = np.arctan2(z, np.sqrt(x**2 + y**2))

    return np.array([yaw, pitch, roll])


def euler_to_direction(angles):
    """
    Convert euler angles to a direction vector.

    Roll gets discarded as this is a reduction in degrees of freedom.
    """
    yaw, pitch, roll = angles

    x = np.cos(yaw) * np.cos(pitch)
    y = np.sin(yaw) * np.cos(pitch)
    z = np.sin(pitch)

    return np.array([x, y, z])
