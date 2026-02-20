"""Analysis extensions for creating plots and visualizations.

These extensions analyze the atomic trajectory and create interactive plots.
"""

import typing as t

import numpy as np
from pydantic import Field

from zndraw.extensions.abc import Category, Extension


class Analysis(Extension):
    """Base class for all analysis extensions."""

    category: t.ClassVar[Category] = Category.ANALYSIS


class Distance(Analysis):
    """Calculate distances between 2 selected atoms across all frames.

    Creates an interactive plot with bidirectional synchronization:
    - Clicking points in the plot sets the frame
    - Frame changes highlight corresponding points
    """

    def run(self, vis: t.Any) -> None:
        import pandas as pd
        import plotly.express as px

        atoms_lst = vis[:]
        distances = []

        selection = list(vis.selection)

        if len(selection) != 2:
            raise ValueError("Please select exactly 2 atoms")

        for atoms in atoms_lst:
            distances.append(atoms.get_distance(selection[0], selection[1], mic=True))

        df = pd.DataFrame({"step": list(range(len(atoms_lst))), "distance": distances})

        fig = px.scatter(
            df,
            x="step",
            y="distance",
            labels={"step": "Frame", "distance": "Distance (Angstrom)"},
            title="Distance Over Time",
        )

        # Add line trace for visual continuity
        fig.add_scatter(
            x=df["step"],
            y=df["distance"],
            mode="lines",
            name="trend",
            line=dict(color="rgba(0, 0, 0, 0.1)"),
            hoverinfo="skip",
            showlegend=False,
        )

        meta_step = np.arange(len(atoms_lst))

        # Set up customdata and interactions schema
        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
            selector=dict(mode="markers"),
            meta={
                "interactions": [
                    {
                        "click": "step",
                        "select": "step",
                        "hover": "step",
                    }
                ]
            },
        )

        fig.update_layout(
            dragmode="lasso",
            hovermode="closest",
        )

        vis.figures["Distance"] = fig


class DihedralAngle(Analysis):
    """Calculate dihedral angles for 4 selected atoms across all frames.

    Creates an interactive plot with bidirectional synchronization.
    """

    def run(self, vis: t.Any) -> None:
        import pandas as pd
        import plotly.express as px

        atoms_lst = vis[:]
        dihedral_angles = []

        selection = list(vis.selection)

        if len(selection) != 4:
            raise ValueError("Please select exactly 4 atoms")
        for atoms in atoms_lst:
            dihedral_angles.append(
                atoms.get_dihedrals(indices=[selection], mic=True)[0]
            )
        df = pd.DataFrame(
            {"step": list(range(len(atoms_lst))), "dihedral": dihedral_angles}
        )

        fig = px.scatter(
            df,
            x="step",
            y="dihedral",
            labels={"step": "Frame", "dihedral": "Dihedral Angle (degrees)"},
            title="Dihedral Angle Over Time",
        )

        # Add line trace for visual continuity
        fig.add_scatter(
            x=df["step"],
            y=df["dihedral"],
            mode="lines",
            name="trend",
            line=dict(color="rgba(0, 0, 0, 0.1)"),
            hoverinfo="skip",
            showlegend=False,
        )

        meta_step = np.arange(len(atoms_lst))

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
            selector=dict(mode="markers"),
            meta={
                "interactions": [
                    {
                        "click": "step",
                        "select": "step",
                        "hover": "step",
                    }
                ]
            },
        )

        fig.update_layout(
            dragmode="lasso",
            hovermode="closest",
        )

        vis.figures["DihedralAngle"] = fig


class Properties1D(Analysis):
    """Create scatter plot of a 1D property over frames.

    Supports interactive frame selection via plot interactions.
    """

    value: str = Field(
        ...,
        description="The property value to plot",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props"],
        },
    )

    def run(self, vis: t.Any) -> None:
        import pandas as pd
        import plotly.express as px

        frames = vis.get(slice(None), keys=[self.value])
        values = [frame[self.value] for frame in frames]
        frame_indices = list(range(len(frames)))

        df = pd.DataFrame({"frame": frame_indices, self.value: values})

        fig = px.scatter(
            df,
            x="frame",
            y=self.value,
            labels={
                "frame": "Frame",
                self.value: self.value,
            },
            title=f"Property: {self.value}",
        )

        fig.add_scatter(
            x=df["frame"],
            y=df[self.value],
            mode="lines",
            name="trend",
            line={"color": "rgba(0, 0, 0, 0.1)"},
            hoverinfo="skip",
            showlegend=False,
        )

        meta_step = np.arange(len(frames))

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
            selector={"mode": "markers"},
            meta={
                "interactions": [
                    {
                        "click": "step",
                        "select": "step",
                        "hover": "step",
                    }
                ]
            },
        )

        fig.update_layout(
            dragmode="lasso",
            hovermode="closest",
        )

        vis.figures[f"Properties1D-{self.value}"] = fig


class Properties2D(Analysis):
    """Create 2D scatter plot of two properties with color mapping.

    Supports interactive frame selection via plot interactions.
    """

    x_data: str = Field(
        ...,
        description="Property for x-axis",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "step"],
        },
    )
    y_data: str = Field(
        ...,
        description="Property for y-axis",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "step"],
        },
    )
    color: str = Field(
        ...,
        description="Property for color mapping",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "step"],
        },
    )
    fix_aspect_ratio: bool = Field(
        True,
        description="Fix aspect ratio to 1:1 for equal scaling",
        json_schema_extra={"format": "checkbox"},
    )

    def run(self, vis: t.Any) -> None:
        import pandas as pd
        import plotly.express as px

        # Only fetch the keys we actually need (skip "step" — it's computed)
        fetch_keys = [k for k in {self.x_data, self.y_data, self.color} if k != "step"]
        frames = vis.get(slice(None), keys=fetch_keys) if fetch_keys else []
        num_frames = len(frames) if frames else len(vis)

        def get_property_data(prop_name: str) -> list:
            if prop_name == "step":
                return list(range(num_frames))
            return [frame[prop_name] for frame in frames]

        df = pd.DataFrame(
            {
                self.x_data: get_property_data(self.x_data),
                self.y_data: get_property_data(self.y_data),
                self.color: get_property_data(self.color),
            }
        )

        fig = px.scatter(
            df,
            x=self.x_data,
            y=self.y_data,
            color=self.color,
            labels={
                self.x_data: self.x_data,
                self.y_data: self.y_data,
                self.color: self.color,
            },
            title=f"2D Properties: {self.x_data} vs {self.y_data}",
        )

        if self.fix_aspect_ratio:
            fig.update_yaxes(scaleanchor="x", scaleratio=1)

        meta_step = np.arange(num_frames)

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
            selector={"mode": "markers"},
            meta={
                "interactions": [
                    {
                        "click": "step",
                        "select": "step",
                        "hover": "step",
                    }
                ]
            },
        )

        fig.update_layout(
            dragmode="lasso",
            hovermode="closest",
        )

        vis.figures["Properties2D"] = fig


# Registry of all analysis extensions
analysis: dict[str, type[Analysis]] = {
    "Distance": Distance,
    "DihedralAngle": DihedralAngle,
    "Properties1D": Properties1D,
    "Properties2D": Properties2D,
}
