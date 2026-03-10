"""ZnDraw extensions module.

Provides extension classes for modifiers, analysis, and selections.
"""

from zndraw.extensions.abc import Category, Extension
from zndraw.extensions.analysis import (
    DihedralAngle,
    Distance,
    Properties1D,
    Properties2D,
)
from zndraw.extensions.modifiers import (
    Center,
    ChangeType,
    Delete,
    Duplicate,
    Empty,
    FixAtoms,
    NewCanvas,
    RemoveAtoms,
    Replicate,
    Wrap,
)
from zndraw.extensions.molecule_building import AddFromSMILES, PackBox
from zndraw.extensions.selections import (
    All,
    ConnectedParticles,
    IdenticalSpecies,
    Invert,
    Neighbour,
    NoneSelection,
    Random,
    Range,
    UpdateSelection,
)

__all__ = [
    # ABC
    "Category",
    "Extension",
    # Analysis
    "DihedralAngle",
    "Distance",
    "Properties1D",
    "Properties2D",
    # Modifiers
    "Center",
    "ChangeType",
    "Delete",
    "Duplicate",
    "Empty",
    "FixAtoms",
    "NewCanvas",
    "RemoveAtoms",
    "Replicate",
    "Wrap",
    # Molecule Building
    "AddFromSMILES",
    "PackBox",
    # Selections
    "All",
    "ConnectedParticles",
    "IdenticalSpecies",
    "Invert",
    "Neighbour",
    "NoneSelection",
    "Random",
    "Range",
    "UpdateSelection",
]
