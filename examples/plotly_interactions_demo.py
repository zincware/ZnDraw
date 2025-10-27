"""
Example: Plotly Figure Interactions

This example demonstrates how to create interactive Plotly figures that synchronize
with ZnDraw's 3D visualization. Three common patterns are shown:

1. Step mapping only - clicking points jumps to specific frames
2. Geometry mapping - selecting points highlights 3D objects
3. Multi-dimensional - combining step and geometry interactions
"""

import numpy as np
import plotly.express as px
import plotly.graph_objs as go


def example_1_step_mapping():
    """Example 1: Energy vs Accuracy plot with frame selection.

    Clicking a point jumps to that frame.
    Selecting a region adds frames to the frame selection.
    """
    print("Example 1: Step mapping only")

    # Generate sample data
    n_frames = 10
    pred_energy = np.random.rand(n_frames)
    true_energy = np.random.rand(n_frames)
    frame_indices = np.arange(n_frames)

    # Create scatter plot
    fig = px.scatter(
        x=pred_energy,
        y=true_energy,
        labels={"x": "Predicted Energy", "y": "True Energy"},
        title="Energy vs Accuracy",
    )

    # Add interaction schema
    fig.update_traces(
        customdata=frame_indices,
        meta={
            "interactions": [
                {"click": "step", "select": "step"}  # dimension 0: customdata[0]
            ]
        },
    )

    print("  - Clicking a point sets current frame")
    print("  - Selecting a region adds frames to frame_selection")
    print()
    return fig


def example_2_geometry_mapping():
    """Example 2: Force distribution histogram.

    Each bin represents a range of forces. Clicking selects all forces in that range.
    This maps to the 'forces' geometry in vis.geometries.
    """
    print("Example 2: Geometry mapping (forces)")

    # Generate sample data
    n_forces = 100
    forces = np.random.exponential(scale=0.5, size=n_forces)
    force_ids = np.arange(
        n_forces
    )  # Must be valid indices into vis.geometries["forces"]

    # Create histogram
    fig = px.histogram(
        x=forces, nbins=10, labels={"x": "Force Magnitude"}, title="Force Distribution"
    )

    # Extract the bin values for customdata (bin index)
    # For histograms, we need to map points to their bin index
    customdata = []
    for i, y_value in enumerate(fig.data[0].y):
        # Each y value represents the count in that bin
        # We'll use the bin index as the force ID (this is a simplification)
        force_bin_ids = (
            force_ids[i * 10 : (i + 1) * 10]
            if i * 10 < len(force_ids)
            else force_ids[i * 10 :]
        )
        customdata.extend(force_bin_ids)

    fig.update_traces(
        customdata=force_ids[: len(fig.data[0].x)],
        meta={"interactions": [{"click": "forces", "select": "forces"}]},
    )

    print("  - Clicking a bin selects all forces in that range")
    print("  - Maps to vis.geometries['forces'] selection")
    print()
    return fig


def example_3_multi_dimensional():
    """Example 3: Particle forces scatter plot.

    A more complex example combining frame selection and particle geometry selection.
    This demonstrates the full power of the interaction schema.
    """
    print("Example 3: Multi-dimensional (step + particles)")

    # Generate sample data
    n_steps = 5
    n_particles = 10
    pred_forces = np.random.rand(n_steps * n_particles)
    particle_ids = np.tile(np.arange(n_particles), n_steps)
    steps = np.repeat(np.arange(n_steps), n_particles)

    # Create scatter plot
    fig = px.scatter(
        x=particle_ids,
        y=pred_forces,
        labels={"x": "Particle ID", "y": "Force Magnitude"},
        title="Particle Forces Over Time",
    )

    # Add 2D customdata: [step, particle_id]
    customdata = np.column_stack([steps, particle_ids])

    # Define interactions for each dimension
    fig.update_traces(
        customdata=customdata,
        meta={
            "interactions": [
                {"click": "step", "select": "step"},  # dimension 0
                {"click": "particles", "select": "particles"},  # dimension 1
            ]
        },
    )

    print("  - Clicking a point:")
    print("    - Sets current frame (from particle_ids dimension)")
    print("    - Selects the particle (from steps dimension)")
    print("  - Selecting a region does both actions on all points")
    print()
    return fig


def example_4_sparse_interactions():
    """Example 4: Sparse interactions with selective actions.

    Demonstrates using None to disable specific dimensions,
    and having different actions for click vs select.
    """
    print("Example 4: Sparse interactions")

    # Generate sample data
    n_steps = 5
    n_particles = 10
    force_mags = np.random.rand(n_steps * n_particles)
    particle_ids = np.tile(np.arange(n_particles), n_steps)
    steps = np.repeat(np.arange(n_steps), n_particles)

    fig = px.scatter(x=particle_ids, y=force_mags, title="Sparse Interactions Example")

    # 3D customdata: [step, particle_id, force_magnitude]
    customdata = np.column_stack([steps, particle_ids, force_mags])

    fig.update_traces(
        customdata=customdata,
        meta={
            "interactions": [
                {"click": "step"},  # dimension 0: click only
                {"select": "particles"},  # dimension 1: select only
                None,  # dimension 2: no interaction
            ]
        },
    )

    print("  - Clicking a point only sets the frame (no particle selection)")
    print("  - Selecting a region only selects particles (no frame selection)")
    print("  - Force magnitude dimension has no interaction")
    print()
    return fig


def create_example_figures():
    """Create all example figures."""
    print("\n=== Plotly Interactions Examples ===\n")

    fig1 = example_1_step_mapping()
    fig2 = example_2_geometry_mapping()
    fig3 = example_3_multi_dimensional()
    fig4 = example_4_sparse_interactions()

    print("=== Usage in ZnDraw ===\n")
    print("To use these figures in ZnDraw:")
    print()
    print("  from zndraw import ZnDraw")
    print("  vis = ZnDraw()")
    print()
    print("  # Add an interactive figure")
    print("  vis.figures['energy_plot'] = fig1")
    print()
    print("  # Now the figure responds to interactions!")
    print("  # Clicking points in the plot updates the 3D view")
    print("  # And 3D selections highlight points in the plot")
    print()

    return fig1, fig2, fig3, fig4


if __name__ == "__main__":
    fig1, fig2, fig3, fig4 = create_example_figures()

    # Optionally show the figures
    # fig1.show()
    # fig2.show()
    # fig3.show()
    # fig4.show()
