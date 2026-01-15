"""Tests for room templates."""

import pytest

from zndraw.room_templates import TEMPLATES


def test_empty_template_exists():
    """The 'empty' template must exist."""
    assert "empty" in TEMPLATES


def test_none_template_exists():
    """The 'none' template must exist."""
    assert "none" in TEMPLATES


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
