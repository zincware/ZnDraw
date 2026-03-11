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
    # Molecule Building
    "AddFromSMILES",
    # Selections
    "All",
    # ABC
    "Category",
    # Modifiers
    "Center",
    "ChangeType",
    "ConnectedParticles",
    "Delete",
    # Analysis
    "DihedralAngle",
    "Distance",
    "Duplicate",
    "Empty",
    "Extension",
    "FixAtoms",
    "IdenticalSpecies",
    "Invert",
    "Neighbour",
    "NewCanvas",
    "NoneSelection",
    "PackBox",
    "Properties1D",
    "Properties2D",
    "Random",
    "Range",
    "RemoveAtoms",
    "Replicate",
    "UpdateSelection",
    "Wrap",
]
