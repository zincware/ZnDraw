"""Tests for room templates."""

import pytest

from zndraw.room_templates import TEMPLATES, get_template_names


def test_empty_template_exists():
    """The 'empty' template must exist."""
    assert "empty" in TEMPLATES


def test_get_template_names():
    """get_template_names returns all registered templates."""
    names = get_template_names()
    assert "empty" in names
    assert isinstance(names, list)


@pytest.mark.parametrize("name", list(TEMPLATES.keys()))
def test_template_is_callable(name):
    """Template functions are callable."""
    assert callable(TEMPLATES[name])


class MockVis:
    """Mock ZnDraw-like object for testing templates."""

    def __init__(self):
        self.frames = []

    def append(self, atoms):
        self.frames.append(atoms)

    def extend(self, frames):
        self.frames.extend(frames)


def test_none_template_creates_zero_frames():
    """None template creates no frames (for Python clients uploading own data)."""
    vis = MockVis()
    TEMPLATES["none"](vis)

    assert len(vis.frames) == 0


def test_empty_template_creates_empty_atoms():
    """Empty template creates one empty ase.Atoms."""
    import ase

    vis = MockVis()
    TEMPLATES["empty"](vis)

    assert len(vis.frames) == 1
    assert isinstance(vis.frames[0], ase.Atoms)
    assert len(vis.frames[0]) == 0  # Empty atoms


def test_water_template_creates_molecule():
    """Water template creates H2O molecule."""
    vis = MockVis()
    TEMPLATES["water"](vis)

    assert len(vis.frames) == 1
    # Water has 3 atoms (H, H, O)
    assert len(vis.frames[0]) == 3


def test_ethanol_template_creates_molecule():
    """Ethanol template creates C2H5OH molecule."""
    vis = MockVis()
    TEMPLATES["ethanol"](vis)

    assert len(vis.frames) == 1
    # Ethanol has 9 atoms (C2H5OH)
    assert len(vis.frames[0]) == 9


def test_benzene_template_creates_molecule():
    """Benzene template creates C6H6 molecule."""
    vis = MockVis()
    TEMPLATES["benzene"](vis)

    assert len(vis.frames) == 1
    # Benzene has 12 atoms (C6H6)
    assert len(vis.frames[0]) == 12
