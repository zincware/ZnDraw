import dataclasses
import typing as t
from collections.abc import MutableSequence

import ase
import numpy as np

from zndraw.client import Client, _TemplateValue
from zndraw.extensions import Extension, ExtensionType
from zndraw.utils import atoms_from_dict, atoms_to_dict, update_colors_and_radii


@dataclasses.dataclass
class ZnDraw(MutableSequence):
    url: str
    room: str
    user: str
    auto_pickup_jobs: bool = True
    template: str | None | t.Type[_TemplateValue] = _TemplateValue

    def __post_init__(self):
        self.client = Client(
            url=self.url,
            room=self.room,
            user=self.user,
            auto_pickup_jobs=self.auto_pickup_jobs,
            template=self.template,
        )
        self.client.connect()

    @property
    def step(self) -> int:
        return self.client.step

    @property
    def lock(self):
        return self.client.lock

    @step.setter
    def step(self, value: int):
        self.client.step = value

    def __len__(self) -> int:
        return len(self.client)

    @property
    def settings(self):
        return self.client.settings

    def register_extension(self, *args, **kwargs):
        self.client.register_extension(*args, **kwargs)

    def append(self, atoms: ase.Atoms) -> None:
        """Append an Atoms object to the visualization."""
        if not isinstance(atoms, ase.Atoms):
            raise TypeError("Only ase.Atoms objects are supported")
        update_colors_and_radii(atoms)
        self.client.append(atoms_to_dict(atoms))

    def extend(self, atoms_list: t.Iterable[ase.Atoms]) -> None:
        """Extend the visualization with a list of Atoms objects."""
        dicts = []
        for atoms in atoms_list:
            if not isinstance(atoms, ase.Atoms):
                raise TypeError("Only ase.Atoms objects are supported")
            update_colors_and_radii(atoms)
            dicts.append(atoms_to_dict(atoms))
        self.client.extend(dicts)

    @t.overload
    def __getitem__(self, index: int) -> ase.Atoms: ...
    @t.overload
    def __getitem__(self, index: slice) -> list[ase.Atoms]: ...
    @t.overload
    def __getitem__(self, index: list[int]) -> list[ase.Atoms]: ...
    @t.overload
    def __getitem__(self, index: np.ndarray) -> ase.Atoms | list[ase.Atoms]: ...
    def __getitem__(
        self, index: int | slice | list[int] | np.ndarray
    ) -> ase.Atoms | list[ase.Atoms]:
        """Get an Atoms object by index."""
        data = self.client[index]
        if isinstance(data, list):
            return [atoms_from_dict(d) for d in data]
        return atoms_from_dict(data)


    @t.overload
    def __setitem__(self, index: int, atoms: ase.Atoms) -> None: ...
    @t.overload
    def __setitem__(self, index: slice, atoms: list[ase.Atoms]) -> None: ...
    @t.overload
    def __setitem__(self, index: list[int], atoms: list[ase.Atoms]) -> None: ...
    @t.overload
    def __setitem__(self, index: np.ndarray, atoms: list[ase.Atoms] | ase.Atoms) -> None: ...
    def __setitem__(self, index: int | slice | list[int] | np.ndarray, atoms: list[ase.Atoms] | ase.Atoms) -> None:
        """Set an Atoms object at a specific index."""
        # Handle single atom vs list of atoms
        if isinstance(atoms, list):
            # List of atoms - validate each one
            if not all(isinstance(a, ase.Atoms) for a in atoms):
                raise TypeError("All elements must be ase.Atoms objects")
            # Convert each atom to dict with colors and radii
            dicts = []
            for atom in atoms:
                update_colors_and_radii(atom)
                dicts.append(atoms_to_dict(atom))
            self.client[index] = dicts
        elif isinstance(atoms, ase.Atoms):
            # Single atom
            update_colors_and_radii(atoms)
            self.client[index] = atoms_to_dict(atoms)
        else:
            raise TypeError("Only ase.Atoms objects or lists of ase.Atoms are supported")


    @t.overload
    def __delitem__(self, index: int) -> None: ...
    @t.overload
    def __delitem__(self, index: slice) -> None: ...
    @t.overload
    def __delitem__(self, index: list[int]) -> None: ...
    @t.overload
    def __delitem__(self, index: np.ndarray) -> None: ...
    def __delitem__(self, index: int | slice | list[int] | np.ndarray) -> None:
        """Delete an Atoms object at a specific index."""
        del self.client[index]

    def insert(self, index: int, atoms: ase.Atoms) -> None:
        """Insert an Atoms object at a specific index."""
        if not isinstance(atoms, ase.Atoms):
            raise TypeError("Only ase.Atoms objects are supported")
        update_colors_and_radii(atoms)
        self.client.insert(index, atoms_to_dict(atoms))

    @property
    def selection(self) -> frozenset[int]:
        """Get the current selection of frame indices."""
        return self.client.selection

    @selection.setter
    def selection(self, value: t.Iterable[int] | None):
        """Set the current selection of frame indices."""
        self.client.selection = value

    @property
    def frame_selection(self) -> frozenset[int]:
        """Get the current selection of frame indices."""
        return self.client.frame_selection

    @frame_selection.setter
    def frame_selection(self, value: t.Iterable[int] | None):
        """Set the current selection of frame indices."""
        self.client.frame_selection = value

    @property
    def bookmarks(self) -> dict[int, str]:
        """Get the current bookmarks of frame indices."""
        return self.client.bookmarks
    
    @bookmarks.setter
    def bookmarks(self, value: dict[int, str] | None):
        """Set the current bookmarks of frame indices."""
        self.client.bookmarks = value

    def run(self, *args, **kwargs) -> None:
        """Run an extension on the current selection."""
        self.client.run(*args, **kwargs)

    def log(self, message: str) -> None:
        """Log a message to the client."""
        self.client.log(message)