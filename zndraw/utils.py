import contextlib
import importlib.util
import json
import logging
import pathlib
import socket
import sys
import tempfile
import uuid
from typing import get_type_hints

import ase
import datamodel_code_generator
import numpy as np
import znframe
from decorator import decorator

SHARED = {"atoms": None}


def split_list_into_chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def estimate_max_batch_size_for_socket(frames: list[znframe.Frame]):
    from .settings import GlobalConfig

    max_size = GlobalConfig().max_socket_data_size
    sizes = [
        sys.getsizeof(json.dumps(frame.to_dict(built_in_types=False)))
        for frame in frames
    ]
    largest_frame = max(sizes)
    chunk_size = int(max_size / largest_frame * 0.9)
    return max(chunk_size, 1)


def rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


def get_radius(value):
    return (0.25 * (2 - np.exp(-0.2 * value)),)


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


@contextlib.contextmanager
def set_global_atoms(atoms: ase.Atoms):
    """Temporarily create a 'SHARED["atoms"]' variable."""
    # TODO: I know this is bad, but I don't know how to do it better - send help!
    SHARED["atoms"] = atoms
    yield
    SHARED["atoms"] = None


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


def hide_discriminator_field(d: dict) -> None:
    """Inplace method to set discriminator field to hidden."""
    raise ValueError("This method is not used anymore!")
    for v in d.values():
        if isinstance(v, dict):
            if v.get("title") == "Discriminator":
                v["options"] = {"hidden": True}
                v["type"] = "string"
            hide_discriminator_field(v)


@decorator
def typecast(func, *args, **kwargs):
    annotations = get_type_hints(func)
    updated_args = []
    updated_kwargs = {}
    for arg, arg_type in zip(args, annotations.values()):
        if not isinstance(arg, arg_type):
            if isinstance(arg, dict):
                updated_args.append(arg_type(**arg))
            else:
                updated_args.append(arg_type(arg))
        else:
            updated_args.append(arg)

    for kwarg, kwarg_type in annotations.items():
        if kwarg in kwargs:
            if not isinstance(kwargs[kwarg], kwarg_type):
                if isinstance(kwargs[kwarg], dict):
                    updated_kwargs[kwarg] = kwarg_type(**kwargs[kwarg])
                else:
                    updated_kwargs[kwarg] = kwarg_type(kwargs[kwarg])
            else:
                updated_kwargs[kwarg] = kwargs[kwarg]
    return func(*updated_args, **updated_kwargs)


@decorator
def typecast_kwargs(func, **kwargs):
    annotations = get_type_hints(func)
    arg_type = list(annotations.values())[0]
    return arg_type(**kwargs)


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
