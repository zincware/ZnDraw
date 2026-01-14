"""Room templates for initializing rooms with predefined content.

Templates are functions that take a ZnDraw-like writer and populate it.
Template names are reserved and cannot be used as room IDs.
"""

from collections.abc import Callable
from typing import Protocol

import ase


class RoomWriter(Protocol):
    """Protocol for objects that can receive atoms frames."""

    def append(self, atoms: ase.Atoms) -> None: ...

    def extend(self, frames: list[ase.Atoms]) -> None: ...


# Registry of all templates
TEMPLATES: dict[str, Callable[[RoomWriter], None]] = {}


def register_template(name: str):
    """Decorator to register a template function.

    Parameters
    ----------
    name : str
        Template name (will be reserved and cannot be used as room ID)
    """

    def decorator(func: Callable[[RoomWriter], None]):
        TEMPLATES[name] = func
        return func

    return decorator


def get_template_names() -> list[str]:
    """Return list of all reserved template names."""
    return list(TEMPLATES.keys())


# --- Built-in Templates ---


@register_template("empty")
def empty(vis: RoomWriter) -> None:
    """Empty room with single empty atoms frame."""
    vis.append(ase.Atoms())


@register_template("water")
def water(vis: RoomWriter) -> None:
    """Room with a water molecule."""
    from ase.build import molecule

    vis.append(molecule("H2O"))


@register_template("ethanol")
def ethanol(vis: RoomWriter) -> None:
    """Room with an ethanol molecule."""
    from ase.build import molecule

    vis.append(molecule("CH3CH2OH"))


@register_template("benzene")
def benzene(vis: RoomWriter) -> None:
    """Room with a benzene molecule."""
    from ase.build import molecule

    vis.append(molecule("C6H6"))
