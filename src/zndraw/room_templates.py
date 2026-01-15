"""Room templates for initializing rooms with predefined content.

Templates are functions that take a writer object and populate it with atoms.
Template names are reserved and cannot be used as room IDs.
"""

import logging
from collections.abc import Callable
from typing import Any

import ase

log = logging.getLogger(__name__)


# Registry of all templates
TEMPLATES: dict[str, Callable[[Any], None]] = {}


def register_template(name: str):
    """Decorator to register a template function.

    Parameters
    ----------
    name : str
        Template name (will be reserved and cannot be used as room ID)
    """

    def decorator(func: Callable[[Any], None]):
        if name in TEMPLATES:
            log.warning(f"Template '{name}' is being overwritten")
        TEMPLATES[name] = func
        return func

    return decorator


# --- Built-in Templates ---


@register_template("none")
def none_template(vis: Any) -> None:
    """Truly empty room with 0 frames.

    Used by Python clients that will upload their own data.
    """
    pass  # No frames added


@register_template("empty")
def empty_template(vis: Any) -> None:
    """Empty room with single empty atoms frame."""
    vis.append(ase.Atoms())
