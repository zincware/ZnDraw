"""
Tests for Plotly figure interactions and validation.

Tests the validate_interaction_schema function from figures_manager module.
"""

import numpy as np
import pytest

from zndraw.figures_manager import validate_interaction_schema

# ============================================================================
# Test: Basic Validation - Valid Cases
# ============================================================================


def test_validate_1d_step_schema_valid():
    """Test valid 1D schema with step interaction."""
    customdata = np.arange(10)
    interactions = [{"click": "step", "select": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
    assert errors == []


def test_validate_1d_geometry_schema_valid():
    """Test valid 1D schema mapping to geometry."""
    customdata = np.arange(10)
    interactions = [{"click": "particles"}]
    geometries = {"particles": {}, "forces": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_2d_schema_valid():
    """Test valid 2D schema with step and geometry."""
    customdata = np.column_stack([np.arange(10), np.arange(10, 20)])
    interactions = [
        {"click": "step", "select": "step"},
        {"click": "particles", "select": "particles"},
    ]
    geometries = {"particles": {}, "forces": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_sparse_schema_with_none():
    """Test schema with None entries (no interaction)."""
    customdata = np.column_stack([np.arange(10), np.arange(10, 20), np.arange(20, 30)])
    interactions = [{"click": "step"}, None, {"click": "forces"}]
    geometries = {"particles": {}, "forces": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_selective_actions_click_only():
    """Test schema with only click action (no select)."""
    customdata = np.arange(10)
    interactions = [{"click": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
    assert errors == []


def test_validate_selective_actions_select_only():
    """Test schema with only select action (no click)."""
    customdata = np.arange(10)
    interactions = [{"select": "particles"}]
    geometries = {"particles": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_hover_action():
    """Test schema with hover action."""
    customdata = np.arange(10)
    interactions = [{"hover": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
    assert errors == []


def test_validate_all_three_actions():
    """Test schema with click, select, and hover."""
    customdata = np.arange(10)
    interactions = [{"click": "step", "select": "step", "hover": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
    assert errors == []


def test_validate_list_customdata():
    """Test with Python list instead of numpy array."""
    customdata = [[i, i + 10] for i in range(10)]  # List of lists
    interactions = [{"click": "step"}, {"click": "particles"}]
    geometries = {"particles": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


# ============================================================================
# Test: Validation Errors - Invalid Cases
# ============================================================================


def test_validate_dimension_mismatch():
    """Test error when schema dimensions don't match customdata."""
    customdata = np.column_stack([np.arange(10), np.arange(10, 20)])  # 2D
    interactions = [{"click": "step"}]  # Only 1 schema entry
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert len(errors) == 1
    assert "Dimension mismatch" in errors[0]
    assert "2 dimension(s)" in errors[0]
    assert "schema has 1" in errors[0]


def test_validate_unknown_geometry():
    """Test error when schema references unknown geometry."""
    customdata = np.arange(10)
    interactions = [{"click": "unknown_geometry"}]
    geometries = {"particles": {}, "forces": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is False
    assert len(errors) == 1
    assert "unknown_geometry" in errors[0]
    assert "not found in vis.geometries" in errors[0]


def test_validate_multiple_unknown_geometries():
    """Test error with multiple unknown geometries in schema."""
    customdata = np.column_stack([np.arange(10), np.arange(10, 20)])
    interactions = [
        {"click": "unknown1"},
        {"click": "unknown2"},
    ]
    geometries = {"particles": {}, "forces": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is False
    assert len(errors) == 2
    assert all("not found" in e for e in errors)


def test_validate_none_customdata():
    """Test error when customdata is None."""
    customdata = None
    interactions = [{"click": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert len(errors) == 1
    assert "None" in errors[0]


def test_validate_empty_customdata():
    """Test error when customdata is empty."""
    customdata = np.array([])
    interactions = [{"click": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert len(errors) == 1
    assert "empty" in errors[0]


def test_validate_empty_list_customdata():
    """Test error when customdata is empty list."""
    customdata = []
    interactions = [{"click": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert len(errors) == 1


def test_validate_completely_empty_schema():
    """Test error when schema is completely empty (all None)."""
    customdata = np.arange(10)
    interactions = [None]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert any("completely empty" in e for e in errors)


def test_validate_invalid_action_value_type():
    """Test error when action value is not a string."""
    customdata = np.arange(10)
    interactions = [{"click": 123}]  # Should be string
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert len(errors) == 1
    assert "must be a string" in errors[0]


def test_validate_invalid_interaction_type():
    """Test error when interaction is not a dict or None."""
    customdata = np.arange(10)
    interactions = ["invalid"]  # Should be dict or None
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert len(errors) == 1
    assert "must be a dict or None" in errors[0]


# ============================================================================
# Test: Edge Cases & Special Scenarios
# ============================================================================


def test_validate_large_customdata():
    """Test validation with large customdata array."""
    customdata = np.arange(100000)
    interactions = [{"click": "step"}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
    assert errors == []


def test_validate_3d_customdata():
    """Test validation with 3D customdata."""
    customdata = np.column_stack(
        [
            np.arange(10),
            np.arange(10, 20),
            np.arange(20, 30),
        ]
    )
    interactions = [
        {"click": "step"},
        {"click": "particles"},
        None,
    ]
    geometries = {"particles": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_without_vis_geometries():
    """Test validation without providing vis.geometries (should fail for geometry actions)."""
    customdata = np.arange(10)
    interactions = [{"click": "particles"}]
    # No vis_geometries provided
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is False
    assert any("not found" in e for e in errors)


def test_validate_multiple_errors():
    """Test that multiple validation errors are collected."""
    customdata = np.column_stack([np.arange(10), np.arange(10, 20)])  # 2D
    interactions = [
        {"click": "unknown1"},
        {"click": "unknown2"},
    ]
    geometries = {}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is False
    assert len(errors) >= 2  # At least geometry errors


def test_validate_step_is_always_valid():
    """Test that 'step' action is always valid (no geometries check needed)."""
    customdata = np.arange(10)
    interactions = [{"click": "step"}]
    # No geometries provided, but step should still be valid
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
    assert errors == []


def test_validate_mixed_step_and_geometry():
    """Test schema mixing step and geometry actions."""
    customdata = np.column_stack(
        [
            np.arange(10),
            np.arange(10, 20),
            np.arange(20, 30),
        ]
    )
    interactions = [
        {"click": "step", "select": "step"},  # step - always valid
        {"click": "particles", "select": "forces"},  # geometry - needs to exist
        {"hover": "step"},  # step - always valid
    ]
    geometries = {"particles": {}, "forces": {}}
    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_null_action_values():
    """Test schema with explicit null action values."""
    customdata = np.arange(10)
    interactions = [{"click": None, "select": None}]
    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
    assert errors == []


# ============================================================================
# Test: Real-World Examples
# ============================================================================


def test_validate_energy_plot_example():
    """Test real-world example: energy prediction plot."""
    n_frames = 100
    frame_indices = np.arange(n_frames)
    interactions = [{"click": "step", "select": "step"}]

    is_valid, errors = validate_interaction_schema(frame_indices, interactions)
    assert is_valid is True
    assert errors == []


def test_validate_force_distribution_example():
    """Test real-world example: force distribution histogram."""
    n_forces = 100
    force_ids = np.arange(n_forces)
    interactions = [{"click": "forces", "select": "forces"}]
    geometries = {"forces": {}}

    is_valid, errors = validate_interaction_schema(force_ids, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_particle_trajectory_example():
    """Test real-world example: particle trajectories over time."""
    n_steps = 10
    n_particles = 50
    steps = np.repeat(np.arange(n_steps), n_particles)
    particle_ids = np.tile(np.arange(n_particles), n_steps)

    customdata = np.column_stack([steps, particle_ids])
    interactions = [
        {"click": "step", "select": "step"},
        {"click": "particles", "select": "particles"},
    ]
    geometries = {"particles": {}}

    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


def test_validate_sparse_example():
    """Test real-world example: sparse interactions."""
    customdata = np.column_stack(
        [
            np.arange(10),
            np.arange(10, 20),
            np.random.rand(10),  # force magnitudes (unused)
        ]
    )
    interactions = [
        {"click": "step"},  # click only
        {"select": "particles"},  # select only
        None,  # no interaction
    ]
    geometries = {"particles": {}}

    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True
    assert errors == []


# ============================================================================
# Test: Parametrized Tests (Testing Multiple Combinations)
# ============================================================================


@pytest.mark.parametrize(
    "action_name",
    ["step", "particles", "forces", "bonds", "atoms"],
)
def test_validate_various_geometry_names(action_name):
    """Test that custom geometry names are accepted when they exist."""
    customdata = np.arange(10)
    interactions = [{"click": action_name}]
    geometries = {action_name: {}}

    is_valid, errors = validate_interaction_schema(customdata, interactions, geometries)
    assert is_valid is True


@pytest.mark.parametrize(
    "action_type",
    ["click", "select", "hover"],
)
def test_validate_each_action_type(action_type):
    """Test that each action type works independently."""
    customdata = np.arange(10)
    interactions = [{action_type: "step"}]

    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True


@pytest.mark.parametrize(
    "dimension",
    [1, 2, 3, 4, 5],
)
def test_validate_various_dimensions(dimension):
    """Test validation with different dimensionalities."""
    if dimension == 1:
        customdata = np.arange(10)
    else:
        customdata = np.column_stack([np.arange(10) for _ in range(dimension)])

    interactions = [{"click": "step"} if i == 0 else None for i in range(dimension)]

    is_valid, errors = validate_interaction_schema(customdata, interactions)
    assert is_valid is True
