import enum
import logging
import typing as t

import ase
import ase.constraints
import numpy as np
from ase.data import chemical_symbols
from pydantic import Field

from zndraw.base import Extension

try:
    from zndraw.modify import extras  # noqa: F401
except ImportError:
    # mdanalysis is not installed
    pass

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw


log = logging.getLogger("zndraw")

Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})


class UpdateScene(Extension):
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
        frames = []
        for _ in range(self.steps):
            # rotate the selected atoms around the vector
            atoms_selected.rotate(angle, vector, center=points[0])
            # update the positions of the selected atoms
            atoms.positions[atom_ids] = atoms_selected.positions
            frames.append(atoms.copy())

        vis.extend(frames)


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


class Translate(UpdateScene):
    """Move the selected atoms along the line."""

    steps: int = Field(10, ge=1)

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        if self.steps > len(vis.segments):
            raise ValueError(
                "The number of steps must be less than the number of segments. You can add more points to increase the number of segments."
            )

        segments = vis.segments
        atoms = vis.atoms
        selection = np.array(vis.selection)

        frames = []

        for idx in range(self.steps):
            end_idx = int((idx + 1) * (len(segments) - 1) / self.steps)
            tmp_atoms = atoms.copy()
            vector = segments[end_idx] - segments[0]
            positions = tmp_atoms.positions
            positions[selection] += vector
            tmp_atoms.positions = positions
            frames.append(tmp_atoms)

        vis.extend(frames)


class Duplicate(UpdateScene):
    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)
    symbol: Symbols = Field(Symbols.X, description="Symbol of the new atoms")

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        atoms = vis.atoms
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        for atom_id in vis.selection:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom
            del atoms.arrays["colors"]
            del atoms.arrays["radii"]
            if hasattr(atoms, "connectivity"):
                del atoms.connectivity

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
        if hasattr(atoms, "connectivity"):
            # vdW radii might change
            del atoms.connectivity

        vis.append(atoms)
        vis.selection = []


class AddLineParticles(UpdateScene):
    symbol: Symbols

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        for point in vis.points:
            atoms += ase.Atom(self.symbol.name, position=point)

        del atoms.arrays["colors"]
        del atoms.arrays["radii"]
        if hasattr(atoms, "connectivity"):
            del atoms.connectivity

        vis.append(atoms)


class Wrap(UpdateScene):
    """Wrap the atoms to the cell."""

    recompute_bonds: bool = True
    all: bool = Field(
        False,
        description="Apply to the full trajectory",
    )

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if self.all:
            for idx, atoms in enumerate(vis):
                atoms.wrap()
                if self.recompute_bonds:
                    delattr(atoms, "connectivity")
                vis[idx] = atoms
        else:
            atoms = vis.atoms
            atoms.wrap()
            if self.recompute_bonds:
                delattr(atoms, "connectivity")
            vis[vis.step] = atoms


class Center(UpdateScene):
    """Move the atoms, such that the selected atom is in the center of the cell."""

    recompute_bonds: bool = True
    dynamic: bool = Field(
        False, description="Move the atoms to the center of the cell at each step"
    )
    wrap: bool = Field(True, description="Wrap the atoms to the cell")
    all: bool = Field(
        False,
        description="Apply to the full trajectory",
    )

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        selection = vis.selection
        if len(selection) < 1:
            vis.log("Please select at least one atom.")
            return

        if not self.dynamic:
            center = vis.atoms[selection].get_center_of_mass()
        else:
            center = None

        if self.all:
            for idx, atoms in enumerate(vis):
                if self.dynamic:
                    center = atoms[selection].get_center_of_mass()
                atoms.positions -= center
                atoms.positions += np.diag(atoms.cell) / 2
                if self.wrap:
                    atoms.wrap()
                if self.recompute_bonds:
                    delattr(atoms, "connectivity")

                vis[idx] = atoms
        else:
            atoms = vis.atoms
            center = atoms[selection].get_center_of_mass()
            atoms.positions -= center
            atoms.positions += np.diag(atoms.cell) / 2
            if self.wrap:
                atoms.wrap()
            if self.recompute_bonds:
                delattr(atoms, "connectivity")

            vis[vis.step] = atoms


class Replicate(UpdateScene):
    """Replicate the atoms in the cell."""

    x: int = Field(2, ge=1, le=10)
    y: int = Field(2, ge=1, le=10)
    z: int = Field(2, ge=1, le=10)

    keep_box: bool = Field(False, description="Keep the original box size")
    all: bool = Field(
        False,
        description="Apply to the full trajectory",
    )

    @classmethod
    def model_json_schema(cls):
        schema = super().model_json_schema()
        # make it sliders
        schema["properties"]["x"]["format"] = "range"
        schema["properties"]["y"]["format"] = "range"
        schema["properties"]["z"]["format"] = "range"
        # and checkboxes
        schema["properties"]["keep_box"]["format"] = "checkbox"
        schema["properties"]["all"]["format"] = "checkbox"
        return schema

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if self.all:
            for idx, atoms in enumerate(vis):
                atoms = atoms.repeat((self.x, self.y, self.z))
                if self.keep_box:
                    atoms.cell = vis[idx].cell
                vis[idx] = atoms
        else:
            atoms = vis.atoms
            atoms = atoms.repeat((self.x, self.y, self.z))
            if self.keep_box:
                atoms.cell = vis.atoms.cell
            vis[vis.step] = atoms


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


class RemoveAtoms(UpdateScene):
    """Remove the current scene."""

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        del vis[vis.step]


class FixAtoms(UpdateScene):
    """Fix the selected atoms."""

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        selection = vis.selection
        atoms = vis.atoms
        if len(selection) == 0:
            atoms.set_constraint()
        else:
            constraint = ase.constraints.FixAtoms(indices=selection)
            atoms.set_constraint(constraint)
        vis.atoms = atoms


modifier: dict[str, t.Type[UpdateScene]] = {
    Delete.__name__: Delete,
    Rotate.__name__: Rotate,
    Translate.__name__: Translate,
    Duplicate.__name__: Duplicate,
    ChangeType.__name__: ChangeType,
    AddLineParticles.__name__: AddLineParticles,
    Wrap.__name__: Wrap,
    Center.__name__: Center,
    Replicate.__name__: Replicate,
    Connect.__name__: Connect,
    NewCanvas.__name__: NewCanvas,
    RemoveAtoms.__name__: RemoveAtoms,
    FixAtoms.__name__: FixAtoms,
}
