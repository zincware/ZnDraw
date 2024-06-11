import importlib.util
import json
import logging
import pathlib
import socket
import sys
import tempfile
import typing as t
import uuid

import datamodel_code_generator
import socketio.exceptions

log = logging.getLogger(__name__)


def get_port(default: int = 1234) -> int:
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


def ensure_path(path: str):
    """Ensure that a path exists."""
    p = pathlib.Path(path).expanduser()
    p.mkdir(parents=True, exist_ok=True)
    return p.as_posix()


def wrap_and_check_index(index: int | slice | list[int], length: int) -> list[int]:
    is_slice = isinstance(index, slice)
    if is_slice:
        index = list(range(*index.indices(length)))
    index = [index] if isinstance(index, int) else index
    index = [i if i >= 0 else length + i for i in index]
    # check if index is out of range
    for i in index:
        if i >= length:
            raise IndexError(f"Index {i} out of range for length {length}")
        if i < 0:
            raise IndexError(f"Index {i-length} out of range for length {length}")
    return index


def check_selection(value: list[int], maximum: int):
    """Check if the selection is valid

    Attributes
    ----------
        value: list[int]
            the selected indices
        maximum: int
            len(vis.step), will be incremented by one, to account for
    """
    if not isinstance(value, list):
        raise ValueError("Selection must be a list")
    if any(not isinstance(i, int) for i in value):
        raise ValueError("Selection must be a list of integers")
    if len(value) != len(set(value)):
        raise ValueError("Selection must be unique")
    if any(i < 0 for i in value):
        raise ValueError("Selection must be positive integers")
    if any(i >= maximum for i in value):
        raise ValueError(
            f"Can not select particles indices larger than size of the scene: {maximum }. Got {value}"
        )


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
        except socketio.exceptions.TimeoutError as err:
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
    return None
