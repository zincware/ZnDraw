"""Output helpers for zndraw-cli.

All JSON output goes to stdout. Errors go to stderr as RFC 9457 problem JSON.
"""

from __future__ import annotations

import json
import sys
from typing import Any

from pydantic import BaseModel


def json_print(data: Any) -> None:
    """Print JSON to stdout.

    Parameters
    ----------
    data
        Data to serialize. Pydantic models are serialized via
        ``model_dump()``, other types via ``json.dumps()``.
    """
    if isinstance(data, BaseModel):
        sys.stdout.write(data.model_dump_json(indent=2) + "\n")
    elif isinstance(data, str):
        sys.stdout.write(data + "\n")
    else:
        sys.stdout.write(json.dumps(data, indent=2, default=str) + "\n")


def text_print(data: str) -> None:
    """Print raw text to stdout."""
    sys.stdout.write(data)
    if not data.endswith("\n"):
        sys.stdout.write("\n")
