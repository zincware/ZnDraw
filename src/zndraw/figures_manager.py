import typing as t
from collections.abc import MutableMapping

import plotly.graph_objs as go
import plotly.io as pio

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Figures(MutableMapping):
    """Manage Plotly figures with interactive customdata and meta.interactions schema.

    This class enables bidirectional synchronization between Plotly plots and 3D geometry:
    - Clicking/selecting points in plots updates frame, particle selection, or other state
    - State changes (frame, selection) highlight corresponding points in plots

    ## Interactive Figures

    To enable interactions, set `customdata` on traces and define an `interactions` schema
    in `meta`. The frontend uses this schema to dispatch state updates.

    ### Schema Format: `meta.interactions`

    Each element in the `interactions` array corresponds to a dimension in `customdata`:

    ```python
    import numpy as np
    import plotly.express as px

    # Example 1: Step mapping only
    fig = px.scatter(x=pred_energy, y=true_energy)
    fig.update_traces(
        customdata=frame_indices,
        meta={
            "interactions": [
                {"click": "step", "select": "step"}  # dimension 0: customdata[0]
            ]
        }
    )
    # Clicking a point sets current frame
    # Selecting region sets frame_selection

    # Example 2: Multi-dimensional (step + particles geometry)
    fig.update_traces(
        customdata=np.column_stack([steps, particle_ids]),
        meta={
            "interactions": [
                {"click": "step", "select": "step"},              # dimension 0
                {"click": "particles", "select": "particles"}     # dimension 1
            ]
        }
    )
    # Clicking point selects both the frame and the particle
    ```

    ### Standard Action Names

    - `"step"`: Controls the current frame or frame selection
      - `click: "step"` - Jump to frame
      - `select: "step"` - Add frames to selection
      - `hover: "step"` - Hover feedback

    - `"<geometry_name>"`: Maps to a geometry in vis.geometries (e.g., `"particles"`, `"forces"`)
      - `click: "<geometry_name>"` - Select this item exclusively
      - `select: "<geometry_name>"` - Select items matching lasso/box
      - `hover: "<geometry_name>"` - Hover feedback

    ### Validation

    Use `validate_interaction_schema()` to verify customdata and schema alignment:

    ```python
    from zndraw import ZnDraw
    import numpy as np

    vis = ZnDraw()
    fig = px.scatter(...)
    customdata = np.column_stack([steps, particle_ids])
    interactions = [{"click": "step"}, {"click": "particles"}]

    is_valid, errors = validate_interaction_schema(
        customdata, interactions, vis.geometries
    )
    if not is_valid:
        print(f"Validation errors: {errors}")
    ```
    """

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def _dict_to_figure(self, figure: dict) -> go.Figure:
        if not isinstance(figure, dict):
            raise ValueError("Figure data must be a dictionary")

        if figure.get("type") == "plotly":
            fig = pio.from_json(figure["data"])
            return fig

        raise ValueError("Unsupported figure type or invalid figure data")

    def _figure_to_dict(self, figure: go.Figure) -> dict:
        if not isinstance(figure, go.Figure):
            raise ValueError("Only plotly.graph_objs.Figure instances are supported")
        response = figure.to_json()
        if response is None:
            raise ValueError("Failed to convert figure to JSON")
        return {"data": response, "type": "plotly"}

    def __getitem__(self, key: str) -> go.Figure:
        if key not in self.vis._figures:
            response = self.vis.api.get_figure(key=key)
            if response is None:
                raise KeyError(f"Figure with key '{key}' does not exist")
            self.vis._figures[key] = response
        if key not in self.vis._figures:
            raise KeyError(f"Figure with key '{key}' does not exist")
        return self._dict_to_figure(figure=self.vis._figures[key])

    def __setitem__(self, key: str, value: t.Any) -> None:
        self.vis._figures[key] = self._figure_to_dict(figure=value)
        self.vis.api.add_figure(
            key=key,
            figure=self.vis._figures[key],
        )

    def __delitem__(self, key: str) -> None:
        self.vis._figures.pop(key, None)
        self.vis.api.delete_figure(key=key)

    def __iter__(self):
        return iter(self.vis.api.list_figures())

    def __len__(self) -> int:
        return len(self.vis.api.list_figures())


def validate_interaction_schema(
    customdata: t.Any,
    interactions_meta: t.Any,
    vis_geometries: dict[str, t.Any] | None = None,
) -> tuple[bool, list[str]]:
    """Validate that meta.interactions aligns with customdata and available geometries.

    Parameters
    ----------
    customdata
        The array passed to fig.update_traces(customdata=...). Can be 1D or 2D.
    interactions_meta
        The schema from meta={"interactions": [...]}.
    vis_geometries
        Available geometries from vis.geometries. If None, skips geometry validation.

    Returns
    -------
    is_valid, errors_list
        Tuple of (success: bool, error_messages: list[str])

    Raises
    ------
    ValueError
        If the schema is invalid.

    Examples
    --------
    >>> import numpy as np
    >>> customdata = np.column_stack([[0, 1, 2], [10, 11, 12]])
    >>> interactions = [{"click": "step"}, {"click": "particles"}]
    >>> is_valid, errors = validate_interaction_schema(customdata, interactions)
    >>> print(is_valid, errors)
    (True, [])

    >>> # Dimension mismatch
    >>> interactions = [{"click": "step"}]  # Missing second dimension
    >>> is_valid, errors = validate_interaction_schema(customdata, interactions)
    >>> print(is_valid)
    False
    """
    if vis_geometries is None:
        vis_geometries = {}

    errors = []

    # Check 1: customdata is not empty
    try:
        if customdata is None:
            errors.append("customdata is None")
            return False, errors
        if len(customdata) == 0:
            errors.append("customdata is empty")
            return False, errors
    except TypeError:
        errors.append("customdata is not a sequence")
        return False, errors

    # Check 2: Determine expected dimensions
    try:
        # Handle 2D customdata (numpy array or list of lists/tuples)
        first_element = customdata[0]
        # Check if it's a sequence (list, tuple, or numpy array)
        if isinstance(first_element, (list, tuple)):
            expected_dims = len(first_element)
        elif hasattr(first_element, "__len__") and not isinstance(first_element, str):
            # numpy array or similar
            try:
                expected_dims = len(first_element)
            except TypeError:
                # Scalar value, 1D customdata
                expected_dims = 1
        else:
            # Scalar value, 1D customdata
            expected_dims = 1
    except (TypeError, IndexError):
        errors.append("customdata has invalid structure")
        return False, errors

    # Check 3: Dimension alignment
    if len(interactions_meta) != expected_dims:
        errors.append(
            f"Dimension mismatch: customdata has {expected_dims} dimension(s), "
            f"schema has {len(interactions_meta)}"
        )

    # Check 4: Valid action names
    valid_actions = {"click", "select", "hover"}
    for idx, interaction in enumerate(interactions_meta):
        if interaction is None:
            continue

        if not isinstance(interaction, dict):
            errors.append(
                f"Schema dimension {idx}: interaction must be a dict or None, got {type(interaction)}"
            )
            continue

        for action_type in valid_actions:
            action_name = interaction.get(action_type)
            if action_name is None:
                continue

            if not isinstance(action_name, str):
                errors.append(
                    f"Schema dimension {idx}, action '{action_type}': "
                    f"action value must be a string, got {type(action_name)}"
                )
                continue

            # Check if action is "step" (always valid) or a known geometry
            if action_name == "step":
                continue  # Valid

            if action_name not in vis_geometries:
                available = list(vis_geometries.keys())
                errors.append(
                    f"Schema dimension {idx}, action '{action_type}': "
                    f"'{action_name}' not found in vis.geometries. "
                    f"Available: {available}"
                )

    # Check 5: Schema not completely empty
    if not interactions_meta or all(x is None for x in interactions_meta):
        errors.append("Interaction schema is completely empty (all None)")

    return len(errors) == 0, errors
