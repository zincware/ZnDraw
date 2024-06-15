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
import socketio.exceptions
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
        return ASEDict(
            numbers=obj.numbers.tolist(),
            positions=obj.positions.tolist(),
            connectivity=None,
            arrays={},
            info={},
            # calc=obj.calc.results,
            pbc=obj.pbc.tolist(),
            cell=obj.cell.tolist(),
        )

    def decode(self, value: ASEDict) -> ase.Atoms:
        """Create datetime object from str / isoformat"""
        return ase.Atoms(
            numbers=value["numbers"],
            positions=value["positions"],
            # connectivity=value["connectivity"],
            # arrays=value["arrays"],
            info=value["info"],
            # calc=value["calc"],
            pbc=value["pbc"],
            cell=value["cell"],
        )


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
