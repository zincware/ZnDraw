import typing as t

import numpy as np
import pandas as pd
import plotly.express as px
from pydantic import Field

from zndraw.extensions.abc import Extension, Category

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Analysis(Extension):
    """The base class for all analysis extensions."""

    category: t.ClassVar[Category] = Category.ANALYSIS


class DihedralAngle(Analysis):
    """Dihedral angle analysis with interactive frame selection.

    Calculates dihedral angles for 4 selected atoms across all frames.
    Creates an interactive plot with bidirectional synchronization:
    - Clicking/selecting points in the plot sets the frame
    - Frame changes in the 3D view highlight corresponding points in the plot
    """

    def run(self, vis: "ZnDraw") -> None:
        """Create interactive dihedral angle plot.

        Requires exactly 4 atoms to be selected in the 3D view.
        """
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

        # Create scatter plot with line for visual continuity
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

        # Set up customdata and interactions schema
        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
            selector=dict(mode="markers"),  # Only update scatter points
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

        # Enable lasso and box selection
        fig.update_layout(
            dragmode="lasso",
            hovermode="closest",
        )

        vis.figures["DihedralAngle"] = fig


class Properties1D(Analysis):
    """1D property analysis with interactive frame selection.

    Creates an interactive scatter plot of 1D properties over frames with
    bidirectional synchronization:
    - Clicking/selecting points in the plot sets the frame
    - Frame changes in the 3D view highlight corresponding points in the plot

    Supports lasso and box selection for selecting multiple frames at once.
    """

    value: str = Field(..., description="The property value")

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        schema["properties"]["value"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["value"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["value"]["type"] = "string"
        schema["properties"]["value"].pop("anyOf", None)
        return schema

    def run(self, vis: "ZnDraw") -> None:
        """Create interactive scatter plot with frame interactions.

        The plot uses:
        - customdata to store frame indices (step)
        - meta.interactions schema to enable frame selection/hovering
        - lasso dragmode for flexible selection
        """
        properties = vis.get(slice(None, None, None), keys=[self.value])
        df = pd.DataFrame(properties)

        # Add frame index for x-axis (continuous positioning)
        frame_indices = np.arange(len(properties))
        df.insert(0, "frame", frame_indices)

        # Create scatter plot with line connecting points for visual continuity
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

        # Add line traces for visual continuity (lower opacity)
        fig.add_scatter(
            x=df["frame"],
            y=df[self.value],
            mode="lines",
            name="trend",
            line=dict(color="rgba(0, 0, 0, 0.1)"),
            hoverinfo="skip",
            showlegend=False,
        )

        # Set up customdata for interactions (frame indices only)
        meta_step = np.arange(len(properties))

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
            selector=dict(mode="markers"),  # Only update scatter points, not line
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

        # Enable lasso and box selection
        fig.update_layout(
            dragmode="lasso",  # Default to lasso, users can toggle to box with toolbar
            hovermode="closest",
        )

        vis.figures[f"Properties1D-{self.value}"] = fig


class Properties2D(Analysis):
    """2D property scatter plot with interactive frame selection.

    Creates an interactive scatter plot showing the relationship between two
    properties across all frames with color-coded data points.
    Supports bidirectional synchronization with frame selection.
    """

    x_data: str = Field(..., description="Property for x-axis")
    y_data: str = Field(..., description="Property for y-axis")
    color: str = Field(..., description="Property for color mapping")
    fix_aspect_ratio: bool = Field(
        True, description="Fix aspect ratio to 1:1 for equal scaling"
    )

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        # Mark property fields as dynamic enums
        for prop in ["x_data", "y_data", "color"]:
            schema["properties"][prop]["x-custom-type"] = "dynamic-enum"
            schema["properties"][prop]["x-features"] = ["dynamic-atom-props", "step"]
            schema["properties"][prop]["type"] = "string"
            schema["properties"][prop].pop("anyOf", None)
        # Mark fix_aspect_ratio as checkbox
        schema["properties"]["fix_aspect_ratio"]["format"] = "checkbox"
        return schema

    def run(self, vis: "ZnDraw") -> None:
        """Create interactive 2D property scatter plot.

        Supports special 'step' property to use frame indices.
        All other properties are retrieved from frame data.
        """
        # Determine which properties to fetch (exclude 'step')
        keys_to_fetch = []
        for prop in [self.x_data, self.y_data, self.color]:
            if prop != "step" and prop not in keys_to_fetch:
                keys_to_fetch.append(prop)

        # Get all frames with required properties
        if keys_to_fetch:
            properties = vis.get(slice(None, None, None), keys=keys_to_fetch)
        else:
            # All properties are 'step', just need frame count
            properties = [{} for _ in range(len(vis))]

        num_frames = len(properties)

        # Extract or generate data for each axis
        def get_property_data(prop_name: str) -> np.ndarray:
            if prop_name == "step":
                return np.arange(num_frames)
            else:
                data = np.array([frame[prop_name] for frame in properties])
                return data.reshape(-1)

        x_data = get_property_data(self.x_data)
        y_data = get_property_data(self.y_data)
        color_data = get_property_data(self.color)

        # Create dataframe for plotting
        df = pd.DataFrame(
            {
                self.x_data: x_data,
                self.y_data: y_data,
                self.color: color_data,
            }
        )

        # Create scatter plot with color mapping
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

        # Fix aspect ratio if requested
        if self.fix_aspect_ratio:
            fig.update_yaxes(scaleanchor="x", scaleratio=1)

        # Set up customdata and interactions for frame selection
        meta_step = np.arange(num_frames)

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

        # Enable lasso and box selection
        fig.update_layout(
            dragmode="lasso",
            hovermode="closest",
        )

        vis.figures["Properties2D"] = fig


class ForceCorrelation(Analysis):
    """Analyze correlation between two properties for the current frame.

    Creates an interactive scatter plot showing the correlation between two
    properties (e.g., force magnitudes) across particles in the current frame.
    Clicking/selecting particles in the plot highlights them in the 3D view.
    """

    x_data: str = Field(..., description="Property for x-axis")
    y_data: str = Field(..., description="Property for y-axis")

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        # Mark these as dynamic enums for the UI to populate dynamically
        for prop in ["x_data", "y_data"]:
            schema["properties"][prop]["x-custom-type"] = "dynamic-enum"
            schema["properties"][prop]["x-features"] = ["dynamic-atom-props"]
            schema["properties"][prop]["type"] = "string"
            schema["properties"][prop].pop("anyOf", None)
        return schema

    def run(self, vis: "ZnDraw") -> None:
        """Create interactive correlation plot for the current frame.

        Maps particle indices to the 'particles' geometry for interactive
        particle selection in the 3D view.
        """
        step = vis.step

        # Get property data for current frame only
        frame_slice = slice(step, step + 1)
        frame_list = vis.get(frame_slice, keys=[self.x_data, self.y_data])

        if not frame_list:
            raise ValueError("No data available for current frame")

        # Extract values from the single frame dict
        frame_data = frame_list[0]
        if self.x_data not in frame_data or self.y_data not in frame_data:
            raise ValueError(
                f"Properties '{self.x_data}' or '{self.y_data}' not found in frame data"
            )

        x_values = frame_data[self.x_data]
        y_values = frame_data[self.y_data]

        # Normalize to get magnitudes if vectors
        if x_values.ndim > 1:
            x_values = np.linalg.norm(x_values, axis=-1)
        if y_values.ndim > 1:
            y_values = np.linalg.norm(y_values, axis=-1)

        x_values = x_values.reshape(-1)
        y_values = y_values.reshape(-1)

        # Create dataframe for plotting
        df = pd.DataFrame(
            {
                self.x_data: x_values,
                self.y_data: y_values,
            }
        )

        # Create scatter plot
        fig = px.scatter(
            df,
            x=self.x_data,
            y=self.y_data,
            labels={
                self.x_data: self.x_data,
                self.y_data: self.y_data,
            },
            title=f"Correlation: {self.x_data} vs {self.y_data} (Step {step})",
            render_mode="svg",
        )

        # Set up interactions: customdata is particle indices, mapped to "particles" geometry
        particle_indices = np.arange(len(x_values))

        fig.update_traces(
            customdata=particle_indices,
            selector=dict(mode="markers"),
            meta={
                "interactions": [
                    {
                        "click": "particles",
                        "select": "particles",
                    }
                ]
            },
        )

        # Enable lasso and box selection for flexible particle selection
        fig.update_layout(
            dragmode="lasso",
            hovermode="closest",
        )

        vis.figures[f"ForceCorrelation-{step}"] = fig


analysis: dict[str, t.Type[Analysis]] = {
    "DihedralAngle": DihedralAngle,
    "Properties1D": Properties1D,
    "Properties2D": Properties2D,
    "ForceCorrelation": ForceCorrelation,
}
