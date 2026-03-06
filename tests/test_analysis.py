"""TDD tests for analysis extensions.

Tests for Distance, DihedralAngle, Properties1D, Properties2D.
Calls extension.run(vis) directly against a real server to verify
that analysis extensions create correct Plotly figures.
"""

import uuid

import ase
import numpy as np
import plotly.graph_objects as go
import pytest

from zndraw import ZnDraw
from zndraw.extensions.analysis import (
    DihedralAngle,
    Distance,
    Properties1D,
    Properties2D,
)
from zndraw.schemas import FigureData


class TestFigureDataSerialization:
    """Unit tests for FigureData.from_figure / .to_figure Pydantic model."""

    def test_roundtrip_preserves_traces(self):
        """Figure roundtrips with correct trace count and type."""
        fig = go.Figure(data=[go.Scatter(x=[1, 2, 3], y=[4, 5, 6])])
        fd = FigureData.from_figure(fig)
        restored = fd.to_figure()
        assert isinstance(restored, go.Figure)
        assert len(restored.data) == 1
        assert restored.data[0].type == "scatter"

    def test_roundtrip_preserves_layout(self):
        """Figure roundtrip preserves layout properties."""
        fig = go.Figure(
            data=[go.Scatter(x=[1], y=[2])],
            layout=go.Layout(title="Test Title", dragmode="lasso"),
        )
        fd = FigureData.from_figure(fig)
        restored = fd.to_figure()
        assert restored.layout.title.text == "Test Title"
        assert restored.layout.dragmode == "lasso"

    def test_roundtrip_preserves_meta(self):
        """Figure roundtrip preserves trace meta (interactions)."""
        fig = go.Figure(data=[go.Scatter(x=[1], y=[2])])
        fig.update_traces(meta={"interactions": [{"click": "step"}]})
        fd = FigureData.from_figure(fig)
        restored = fd.to_figure()
        assert restored.data[0].meta == {"interactions": [{"click": "step"}]}

    def test_roundtrip_with_numpy_data(self):
        """Figure with numpy data roundtrips to a valid figure."""
        fig = go.Figure(
            data=[
                go.Scatter(
                    x=np.array([1.0, 2.0, 3.0]),
                    y=np.array([4.0, 5.0, 6.0]),
                )
            ]
        )
        fd = FigureData.from_figure(fig)
        restored = fd.to_figure()
        assert isinstance(restored, go.Figure)
        assert len(restored.data) == 1

    def test_type_defaults_to_plotly(self):
        """Default type is 'plotly'."""
        fd = FigureData.from_figure(go.Figure())
        assert fd.type == "plotly"

    def test_model_dump_structure(self):
        """model_dump produces expected wire format."""
        fig = go.Figure(data=[go.Scatter(x=[1], y=[2])])
        fd = FigureData.from_figure(fig)
        dumped = fd.model_dump()
        assert dumped["type"] == "plotly"
        assert isinstance(dumped["data"], str)


def _make_trajectory(n_frames: int = 5) -> list[ase.Atoms]:
    """Create a trajectory of H2 molecules with increasing bond length.

    Frame i has atom 0 at origin, atom 1 at (1+i, 0, 0).
    Distance increases from 1.0 to 1+n_frames-1.
    """
    frames = []
    for i in range(n_frames):
        atoms = ase.Atoms("HH", positions=[[0, 0, 0], [1.0 + i, 0, 0]])
        atoms.info["energy"] = -10.0 + i * 0.5
        atoms.info["temperature"] = 300 + i * 10
        frames.append(atoms)
    return frames


def _make_dihedral_trajectory(n_frames: int = 5) -> list[ase.Atoms]:
    """Create a trajectory with 4 atoms suitable for dihedral angle measurement.

    Atoms 0-1-2-3 form a dihedral that rotates with frame index.
    """
    frames = []
    for i in range(n_frames):
        angle = np.radians(i * 30)  # 0°, 30°, 60°, 90°, 120°
        positions = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.5, 1.0, 0.0],
            [2.0, 1.5, np.sin(angle)],
        ]
        atoms = ase.Atoms("HHHH", positions=positions)
        frames.append(atoms)
    return frames


class TestDistance:
    """Tests for Distance analysis extension."""

    def test_distance_creates_figure(self, server: str):
        """Distance.run() creates a Plotly figure named 'Distance' in vis.figures."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(5))
            vis.selection = [0, 1]

            Distance().run(vis)

            assert "Distance" in vis.figures
            fig = vis.figures["Distance"]
            assert isinstance(fig, go.Figure)

    def test_distance_figure_has_correct_traces(self, server: str):
        """Distance figure has scatter + line traces for 5 frames."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(5))
            vis.selection = [0, 1]

            Distance().run(vis)

            fig = vis.figures["Distance"]
            # First trace is scatter (markers), second is line (trend)
            assert len(fig.data) == 2
            assert fig.data[0].type == "scatter"

    def test_distance_figure_has_interaction_metadata(self, server: str):
        """Distance figure has interactions schema for frame synchronization."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(3))
            vis.selection = [0, 1]

            Distance().run(vis)

            fig = vis.figures["Distance"]
            scatter = fig.data[0]
            assert scatter.meta is not None
            interactions = scatter.meta["interactions"]
            assert len(interactions) == 1
            assert interactions[0]["click"] == "step"

    def test_distance_wrong_selection_count(self, server: str):
        """Distance raises ValueError when selection doesn't have exactly 2 atoms."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(3))
            vis.selection = [0, 1, 2]  # 3 atoms, need 2

            with pytest.raises(ValueError, match="exactly 2"):
                Distance().run(vis)


class TestDihedralAngle:
    """Tests for DihedralAngle analysis extension."""

    def test_dihedral_creates_figure(self, server: str):
        """DihedralAngle.run() creates a Plotly figure named 'DihedralAngle'."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_dihedral_trajectory(5))
            vis.selection = [0, 1, 2, 3]

            DihedralAngle().run(vis)

            assert "DihedralAngle" in vis.figures
            fig = vis.figures["DihedralAngle"]
            assert isinstance(fig, go.Figure)

    def test_dihedral_figure_has_interaction_metadata(self, server: str):
        """DihedralAngle figure has interactions schema for frame sync."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_dihedral_trajectory(3))
            vis.selection = [0, 1, 2, 3]

            DihedralAngle().run(vis)

            fig = vis.figures["DihedralAngle"]
            scatter = fig.data[0]
            assert scatter.meta is not None
            assert scatter.meta["interactions"][0]["click"] == "step"

    def test_dihedral_wrong_selection_count(self, server: str):
        """DihedralAngle raises ValueError with wrong selection."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_dihedral_trajectory(3))
            vis.selection = [0, 1]  # 2 atoms, need 4

            with pytest.raises(ValueError, match="exactly 4"):
                DihedralAngle().run(vis)


class TestProperties1D:
    """Tests for Properties1D analysis extension."""

    def test_properties1d_uses_metadata_keys(self, server: str):
        """Properties1D uses keys from the /metadata endpoint (e.g. 'info.energy')."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(5))

            # Discover keys from metadata endpoint (same as frontend dropdown)
            resp = vis.api.http.get(
                f"/v1/rooms/{vis.room}/frames/0/metadata",
                headers=vis.api._headers(),
            )
            metadata_keys = list(resp.json()["metadata"].keys())
            assert "info.energy" in metadata_keys

            Properties1D(value="info.energy").run(vis)

            assert "Properties1D-info.energy" in vis.figures
            fig = vis.figures["Properties1D-info.energy"]
            assert isinstance(fig, go.Figure)

    def test_properties1d_figure_has_correct_traces(self, server: str):
        """Properties1D figure has scatter + line traces."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(5))

            Properties1D(value="info.energy").run(vis)

            fig = vis.figures["Properties1D-info.energy"]
            assert len(fig.data) == 2
            assert fig.data[0].type == "scatter"

    def test_properties1d_figure_has_interaction_metadata(self, server: str):
        """Properties1D figure has interactions schema for frame sync."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(3))

            Properties1D(value="info.energy").run(vis)

            fig = vis.figures["Properties1D-info.energy"]
            scatter = fig.data[0]
            assert scatter.meta is not None
            assert scatter.meta["interactions"][0]["click"] == "step"


class TestProperties2D:
    """Tests for Properties2D analysis extension."""

    def test_properties2d_uses_metadata_keys(self, server: str):
        """Properties2D uses keys from the /metadata endpoint."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(5))

            resp = vis.api.http.get(
                f"/v1/rooms/{vis.room}/frames/0/metadata",
                headers=vis.api._headers(),
            )
            metadata_keys = list(resp.json()["metadata"].keys())
            assert "info.energy" in metadata_keys
            assert "info.temperature" in metadata_keys

            Properties2D(
                x_data="info.energy",
                y_data="info.temperature",
                color="step",
                fix_aspect_ratio=False,
            ).run(vis)

            assert "Properties2D" in vis.figures
            fig = vis.figures["Properties2D"]
            assert isinstance(fig, go.Figure)

    def test_properties2d_with_step_as_axis(self, server: str):
        """Properties2D supports 'step' as a special axis value for frame indices."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(5))

            Properties2D(
                x_data="step",
                y_data="info.energy",
                color="step",
                fix_aspect_ratio=False,
            ).run(vis)

            fig = vis.figures["Properties2D"]
            assert isinstance(fig, go.Figure)
            assert len(fig.data) >= 1

    def test_properties2d_figure_has_interaction_metadata(self, server: str):
        """Properties2D figure has interactions schema for frame sync."""
        with ZnDraw(url=server, room=uuid.uuid4().hex) as vis:
            vis.extend(_make_trajectory(3))

            Properties2D(
                x_data="info.energy",
                y_data="info.temperature",
                color="step",
                fix_aspect_ratio=False,
            ).run(vis)

            fig = vis.figures["Properties2D"]
            scatter = fig.data[0]
            assert scatter.meta is not None
            assert scatter.meta["interactions"][0]["click"] == "step"
