import dataclasses
from zndraw.client import Client
from zndraw.extensions import Extension, ExtensionType
import typing as t
import ase
from zndraw.utils import update_colors_and_radii, atoms_to_dict, atoms_from_dict
from collections.abc import MutableSequence



@dataclasses.dataclass
class ZnDraw(MutableSequence):
    url: str
    room: str
    user: str
    auto_pickup_jobs: bool = True

    def __post_init__(self):
        self.client = Client(url=self.url, room=self.room, user=self.user, auto_pickup_jobs=self.auto_pickup_jobs)
        self.client.connect()

    @property
    def step(self) -> int:
        return self.client.step

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
    def __getitem__(self, index: int|slice|list[int]) -> ase.Atoms|list[ase.Atoms]:
        """Get an Atoms object by index."""
        data = self.client[index]
        return atoms_from_dict(data)
    
    def __setitem__(self, index: int, atoms: ase.Atoms) -> None:
        """Set an Atoms object at a specific index."""
        if not isinstance(atoms, ase.Atoms):
            raise TypeError("Only ase.Atoms objects are supported")
        update_colors_and_radii(atoms)
        self.client[index] = atoms_to_dict(atoms)

    def __delitem__(self, index: int) -> None:
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
