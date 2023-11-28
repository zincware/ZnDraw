import contextlib
import importlib.util
import json
import logging
import pathlib
import socket
import sys
import tempfile
import uuid

import ase
import datamodel_code_generator
import numpy as np

SHARED = {"atoms": None}


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
        port = 1234
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
            self.handleError(record)


def get_cls_from_json_schema(schema: dict, name: str, **kwargs):
    """Get a python class from a json schema."""

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
    for v in d.values():
        if isinstance(v, dict):
            if v.get("title") == "Discriminator":
                v["options"] = {"hidden": True}
                v["type"] = "string"
            hide_discriminator_field(v)
