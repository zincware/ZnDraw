import abc
import enum
import logging
import time
import typing as t

import ase
import numpy as np
from ase.data import chemical_symbols
from pydantic import Field

from zndraw.base import Extension, MethodsCollection

try:
    from zndraw.modify import extras  # noqa: F401
except ImportError:
    # mdanalysis is not installed
    pass

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw


log = logging.getLogger("zndraw")

Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})


class UpdateScene(Extension, abc.ABC):
    @abc.abstractmethod
    def run(self, vis: "ZnDraw", timeout: float, **kwargs) -> None:
        """Method called when running the modifier."""
        pass

    def apply_selection(
        self, atom_ids: list[int], atoms: ase.Atoms
    ) -> t.Tuple[ase.Atoms, ase.Atoms]:
        """Split the atoms object into the selected and remaining atoms."""
        atoms_selected = atoms[atom_ids]
        atoms_remaining_ids = [x for x in range(len(atoms)) if x not in atom_ids]
        if len(atoms_remaining_ids) > 0:
            atoms_remaining = atoms[atoms_remaining_ids]
        else:
            atoms_remaining = ase.Atoms()
        return atoms_selected, atoms_remaining


class Connect(UpdateScene):
    """Create guiding curve between selected atoms."""

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        atoms = vis.atoms
        atom_ids = vis.selection
        atom_positions = vis.atoms.get_positions()
        camera_position = np.array(vis.camera["position"])[None, :]  # 1,3

        new_points = atom_positions[atom_ids]  # N, 3
        radii: np.ndarray = atoms.arrays["radii"][atom_ids][:, None]
        direction = camera_position - new_points
        direction /= np.linalg.norm(direction, axis=1, keepdims=True)
        new_points += direction * radii

        vis.points = new_points
        vis.selection = []


class Rotate(UpdateScene):
    """Rotate the selected atoms around a the line (2 points only)."""

    angle: float = Field(90, le=360, ge=0, description="Angle in degrees")
    direction: t.Literal["left", "right"] = Field(
        "left", description="Direction of rotation"
    )
    steps: int = Field(
        30, ge=1, description="Number of steps to take to complete the rotation"
    )
    sleep: float = Field(0.1, ge=0, description="Sleep time between steps")

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        # split atoms object into the selected from atoms_ids and the remaining
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        points = vis.points
        atom_ids = vis.selection
        atoms = vis.atoms
        if len(points) != 2:
            raise ValueError("Please draw exactly 2 points to rotate around.")

        angle = self.angle if self.direction == "left" else -self.angle
        angle = angle / self.steps

        atoms_selected, atoms_remaining = self.apply_selection(atom_ids, atoms)
        # create a vector from the two points
        vector = points[1] - points[0]
        for _ in range(self.steps):
            # rotate the selected atoms around the vector
            atoms_selected.rotate(angle, vector, center=points[0])
            # update the positions of the selected atoms
            atoms.positions[atom_ids] = atoms_selected.positions
            vis.append(atoms)
            time.sleep(self.sleep)


class Delete(UpdateScene):
    """Delete the selected atoms."""

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        atom_ids = vis.selection
        atoms = vis.atoms

        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]
        vis.log(f"Deleting atoms {atom_ids}")
        if len(atom_ids) == len(atoms):
            vis.append(ase.Atoms())
        else:
            for idx, atom_id in enumerate(sorted(atom_ids)):
                atoms.pop(atom_id - idx)  # we remove the atom and shift the index
            if hasattr(atoms, "connectivity"):
                del atoms.connectivity
        vis.append(atoms)
        vis.selection = []
        vis.step += 1


class Move(UpdateScene):
    """Move the selected atoms along the line."""

    steps: int = Field(10, ge=1)

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        atoms_ids = vis.selection
        atoms_selected, atoms_remaining = self.apply_selection(atoms_ids, atoms)
        if self.steps > len(vis.segments):
            raise ValueError(
                "The number of steps must be less than the number of segments. You can add more points to increase the number of segments."
            )

        segments = vis.segments
        steps = self.steps

        for idx in range(1, steps):
            # get the vector between the two points
            start_idx = int((idx - 1) * len(segments) / steps)
            end_idx = int(idx * len(segments) / steps)

            vector = segments[end_idx] - segments[start_idx]
            # move the selected atoms along the vector
            atoms_selected.positions += vector
            # merge the selected and remaining atoms
            atoms.positions[atoms_ids] = atoms_selected.positions
            vis.append(atoms)
            vis.step += 1


class Duplicate(UpdateScene):
    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)
    symbol: Symbols

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        atoms = vis.atoms
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        for atom_id in vis.selection:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom

        vis.append(atoms)
        vis.selection = []


class ChangeType(UpdateScene):
    symbol: Symbols

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        for atom_id in vis.selection:
            atoms[atom_id].symbol = self.symbol.name

        del atoms.arrays["colors"]
        del atoms.arrays["radii"]

        vis.append(atoms)
        vis.selection = []


class AddLineParticles(UpdateScene):
    symbol: Symbols
    steps: int = Field(10, le=100, ge=1)

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        for point in vis.points:
            atoms += ase.Atom(self.symbol.name, position=point)

        for _ in range(self.steps):
            vis.append(atoms)


class Wrap(UpdateScene):
    """Wrap the atoms to the cell."""

    recompute_bonds: bool = True

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        vis.log("Downloading atoms...")
        atoms_list = list(vis)
        vis.step = 0

        del vis[1:]

        for idx, atoms in enumerate(atoms_list):
            atoms.wrap()
            if self.recompute_bonds:
                delattr(atoms, "connectivity")
            vis[idx] = atoms


class Center(UpdateScene):
    """Move the atoms, such that the selected atom is in the center of the cell."""

    recompute_bonds: bool = True
    dynamic: bool = Field(
        False, description="Move the atoms to the center of the cell at each step"
    )
    wrap: bool = Field(True, description="Wrap the atoms to the cell")

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        selection = vis.selection
        if len(selection) < 1:
            vis.log("Please select at least one atom.")
            return

        vis.log("Downloading atoms...")
        atoms_list = list(vis)

        if not self.dynamic:
            center = atoms_list[vis.step][selection].get_center_of_mass()
        else:
            center = None

        vis.step = 0

        for idx, atoms in enumerate(atoms_list):
            if self.dynamic:
                center = atoms[selection].get_center_of_mass()
            atoms.positions -= center
            atoms.positions += np.diag(atoms.cell) / 2
            if self.wrap:
                atoms.wrap()
            if self.recompute_bonds:
                delattr(atoms, "connectivity")

            vis[idx] = atoms


class Replicate(UpdateScene):
    x: int = Field(2, ge=1)
    y: int = Field(2, ge=1)
    z: int = Field(2, ge=1)

    keep_box: bool = Field(False, description="Keep the original box size")

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        vis.log("Downloading atoms...")
        atoms_list = list(vis)
        vis.step = 0

        del vis[1:]

        for idx, atoms in enumerate(atoms_list):
            atoms = atoms.repeat((self.x, self.y, self.z))
            if self.keep_box:
                atoms.cell = atoms_list[idx].cell
            vis[idx] = atoms


class NewCanvas(UpdateScene):
    """Clear the scene, deleting all atoms and points."""

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        from zndraw.draw import Plane

        del vis[vis.step + 1 :]
        vis.points = []
        vis.append(ase.Atoms())
        vis.selection = []
        step = len(vis) - 1
        vis.step = step
        vis.bookmarks = vis.bookmarks | {step: "New Scene"}
        vis.camera = {"position": [0, 0, -15], "target": [0, 0, 0]}
        vis.geometries = [
            Plane(
                position=[0, 0, 0],
                rotation=[0, 0, 0],
                scale=[1, 1, 1],
                width=10,
                height=10,
            )
        ]


methods = t.Union[
    Delete,
    Rotate,
    Move,
    Duplicate,
    ChangeType,
    AddLineParticles,
    Wrap,
    Center,
    Replicate,
    Connect,
    NewCanvas,
]


class Modifier(MethodsCollection):
    """Run modifications on the scene"""

    method: methods = Field(
        ..., description="Modify method", discriminator="discriminator"
    )
