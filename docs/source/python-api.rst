Python Interface
================

The ``zndraw`` package provides a Python interface to interact with the visualisation tool.
To use this API, you need to have a running instance of the ZnDraw web server.

Getting Started
---------------

Start a local webserver using the command line interface:

.. code:: console

    $ zndraw file.xyz --port 1234

Then connect from Python:

.. code:: python

    from zndraw import ZnDraw

    vis = ZnDraw(url="http://localhost:1234", room="my-room")

.. note::

    Each visualisation is associated with a room name (visible in the URL).
    Use this room name to interact via Python API or share with others.

Click the connection info button in the UI to see how to connect from Python:

.. image:: /_static/screenshots/lightmode/python_connection.png
   :class: only-light
   :alt: Connection info dialog

.. image:: /_static/screenshots/darkmode/python_connection.png
   :class: only-dark
   :alt: Connection info dialog


Authentication
--------------

ZnDraw supports optional user authentication:

.. code:: python

    vis = ZnDraw(
        url="http://localhost:1234",
        room="my-room",
        user="my-username",
        password="my-password"
    )

If no user is provided, the server will assign a guest username.


Working with Frames
-------------------

.. image:: /_static/screenshots/lightmode/timeline.png
   :class: only-light
   :alt: Timeline with multiple frames

.. image:: /_static/screenshots/darkmode/timeline.png
   :class: only-dark
   :alt: Timeline with multiple frames

The ``vis`` object behaves like a Python list of `ase.Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_ objects.
Modifying the list updates the visualisation in real-time.

.. code:: python

    import ase.io as aio

    # Load and display frames
    frames = aio.read("file.xyz", index=":")
    vis.extend(frames)

    # Access frames
    atoms = vis[vis.step]  # Current frame
    subset = vis[10:20]    # Slice of frames

    # Iterate
    for atoms in vis:
        print(atoms)

    # Navigate to a specific frame
    vis.step = 25


Selections
----------

.. image:: /_static/screenshots/lightmode/selection.png
   :class: only-light
   :alt: Selection tools

.. image:: /_static/screenshots/darkmode/selection.png
   :class: only-dark
   :alt: Selection tools

Select atoms by index using ``vis.selection`` (shortcut for particles geometry):

.. code:: python

    # Set selection for particles
    vis.selection = [0, 1, 2, 3]

    # Get selection
    selected = vis.selection

    # Clear selection
    vis.selection = []

Select by geometry type using ``vis.selections`` (dict-like interface):

.. code:: python

    # Set selection for specific geometry
    vis.selections["particles"] = [0, 1, 2, 3]
    vis.selections["forces"] = [5, 6, 7]

    # Get selection for a geometry
    particle_selection = vis.selections["particles"]

    # Clear selection for a geometry
    del vis.selections["particles"]

    # Iterate over geometries with selections
    for geometry in vis.selections:
        print(f"{geometry}: {vis.selections[geometry]}")

Use selection groups to save and restore named selections across multiple geometries:

.. code:: python

    # Create named group (maps geometry names to indices)
    vis.selection_groups["backbone"] = {
        "particles": [0, 1, 2],
        "forces": [0, 1, 2]
    }

    # Create group with only particles
    vis.selection_groups["active_site"] = {"particles": [10, 11, 12]}

    # Get a group
    group = vis.selection_groups["backbone"]
    print(group)  # {"particles": [0, 1, 2], "forces": [0, 1, 2]}

    # Load a group (apply it to current selections)
    vis.load_selection_group("backbone")

    # List all groups
    for group_name in vis.selection_groups:
        print(f"{group_name}: {vis.selection_groups[group_name]}")

    # Delete a group
    del vis.selection_groups["backbone"]


Bookmarks
---------

.. image:: /_static/screenshots/lightmode/bookmarks.png
   :class: only-light
   :alt: Frame bookmarks

.. image:: /_static/screenshots/darkmode/bookmarks.png
   :class: only-dark
   :alt: Frame bookmarks

Label important frames with bookmarks:

.. code:: python

    # Add bookmark to frame 0
    vis.bookmarks[0] = "Initial State"

    # Add bookmark to frame 50
    vis.bookmarks[50] = "Transition"

    # List all bookmarks
    for frame, label in vis.bookmarks.items():
        print(f"Frame {frame}: {label}")

    # Delete bookmark
    del vis.bookmarks[0]


Geometries
----------

.. image:: /_static/screenshots/lightmode/geometries.png
   :class: only-light
   :alt: Geometries in scene

.. image:: /_static/screenshots/darkmode/geometries.png
   :class: only-dark
   :alt: Geometries in scene

Add 3D geometries to the scene using ``vis.geometries``:

.. code:: python

    from zndraw.geometries import Box, Sphere, Curve, Arrow, Floor

    # Add a floor
    vis.geometries["floor"] = Floor(active=True, height=-2.0, color="#808080")

    # Add a red box with cartoon material
    vis.geometries["box"] = Box(
        position=[(0, 2, 0)],
        size=[(4, 4, 4)],
        color=["#e74c3c"],
        material="MeshToonMaterial"
    )

    # Add a blue sphere with glass material
    vis.geometries["sphere"] = Sphere(
        position=[(8, 2, 0)],
        radius=[2.0],
        color=["#3498db"],
        material="MeshPhysicalMaterial_glass"
    )

    # Add a green curve
    vis.geometries["curve"] = Curve(
        position=[(-6, 0, -6), (-3, 4, -3), (0, 0, 0), (3, 4, 3), (6, 0, 6)],
        color="#2ecc71"
    )

    # Add an orange arrow with shiny material
    vis.geometries["arrow"] = Arrow(
        position=[(12, 0, 0)],
        direction=[(0, 5, 0)],
        color=["#f39c12"],
        material="MeshPhysicalMaterial_shiny"
    )

    # List geometries
    print(list(vis.geometries.keys()))

    # Delete geometry
    del vis.geometries["box"]

Available materials: ``MeshPhysicalMaterial_matt`` (default), ``MeshPhysicalMaterial_glass``,
``MeshPhysicalMaterial_shiny``, ``MeshToonMaterial``, ``MeshStandardMaterial_metallic``, and more.


Curve Customization
^^^^^^^^^^^^^^^^^^^

Curves use CatmullRom spline interpolation between control points. Customize the curve
appearance with additional parameters:

.. code:: python

    from zndraw.geometries import Curve, CurveMarker

    vis.geometries["curve"] = Curve(
        position=[(-6, 0, -6), (-3, 4, -3), (0, 0, 0), (3, 4, 3), (6, 0, 6)],
        color="#2ecc71",
        divisions=100,  # Interpolation smoothness (1-200, default 50)
        thickness=3.0,  # Line thickness (0.5-10, default 2.0)
        marker=CurveMarker(
            enabled=True,    # Show control point markers
            size=0.15,       # Marker size (0.01-1.0)
            color="default", # Use curve color, or specify hex
            opacity=1.0,     # Marker opacity (0-1)
        ),
        virtual_marker=CurveMarker(
            enabled=True,    # Show markers between control points
            size=0.1,        # Smaller than main markers
            opacity=0.5,     # Semi-transparent
        ),
    )

**Marker Settings:**

- ``marker``: Settings for control point markers (the main editable points)
- ``virtual_marker``: Settings for markers between control points (shown in editing mode,
  click to insert new control point)

Both marker types support:

- ``enabled``: Show or hide markers
- ``size``: Marker size (0.01 to 1.0)
- ``color``: Hex color or ``"default"`` to use curve color
- ``opacity``: Transparency (0 = invisible, 1 = opaque)
- ``selecting``: Appearance when selected (color, opacity)
- ``hovering``: Appearance when hovered (color, opacity)

Manage geometries through the UI panel:

.. image:: /_static/screenshots/lightmode/geometry_viewer.png
   :class: only-light
   :alt: Geometry management panel

.. image:: /_static/screenshots/darkmode/geometry_viewer.png
   :class: only-dark
   :alt: Geometry management panel


Camera Control
^^^^^^^^^^^^^^

.. image:: /_static/screenshots/lightmode/camera.png
   :class: only-light
   :alt: Camera with curve attachment

.. image:: /_static/screenshots/darkmode/camera.png
   :class: only-dark
   :alt: Camera with curve attachment

Control the camera programmatically. Cameras can use direct coordinates for static positions
or ``CurveAttachment`` to follow curve paths for animations:

.. code:: python

    from zndraw.geometries import Camera, Curve
    from zndraw.transformations import CurveAttachment

    # Static camera with direct coordinates
    vis.geometries["camera"] = Camera(
        position=(0, 5, 10),
        target=(0, 0, 0),
        fov=60
    )

    # Animated camera following curves
    vis.geometries["cam_path"] = Curve(
        position=[(25, 15, 25), (30, 20, 15), (25, 15, 5)],
        color="#3498db"
    )
    vis.geometries["target_path"] = Curve(
        position=[(10, 10, 10), (12, 10, 10), (10, 10, 10)],
        color="#e74c3c"
    )

    vis.geometries["camera"] = Camera(
        position=CurveAttachment(geometry_key="cam_path", progress=0.5),
        target=CurveAttachment(geometry_key="target_path", progress=0.5),
        fov=60,
        helper_visible=True,
        helper_color="#00ff00"
    )

You can mix direct coordinates with ``CurveAttachment``:

.. code:: python

    # Camera follows curve but always looks at origin
    vis.geometries["camera"] = Camera(
        position=CurveAttachment(geometry_key="flight_path", progress=0.0),
        target=(0, 0, 0),  # Fixed target
        fov=60
    )

Camera parameters:

- ``position``: Direct ``(x, y, z)`` coordinates or ``CurveAttachment``
- ``target``: Direct ``(x, y, z)`` coordinates or ``CurveAttachment``
- ``fov``: Field of view in degrees (1-179, default 50)
- ``camera_type``: ``CameraType.PERSPECTIVE`` or ``CameraType.ORTHOGRAPHIC``
- ``helper_visible``: Show camera cone visualization
- ``helper_color``: Color of the helper (hex or named)
- ``near``, ``far``: Clipping planes
- ``zoom``: Camera zoom factor


Drawing Mode
^^^^^^^^^^^^

.. image:: /_static/screenshots/lightmode/drawing_mode.png
   :class: only-light
   :alt: Drawing mode

.. image:: /_static/screenshots/darkmode/drawing_mode.png
   :class: only-dark
   :alt: Drawing mode

Draw curve control points interactively in the 3D view.

**Entering Drawing Mode:**

1. Press ``X`` to enter drawing mode (from view mode)
2. A drawing marker appears at your cursor position
3. Click to add control points to the active curve
4. Press ``X`` again to exit and return to view mode

The drawing marker shows where a new point will be added. It turns red when
the cursor is over an invalid position.

.. note::

   Drawing mode only works with Curve geometries that have static positions.
   See :doc:`keyboard-shortcuts` for the complete list of keyboard controls.


Editing Mode
^^^^^^^^^^^^

.. image:: /_static/screenshots/lightmode/editing_mode.png
   :class: only-light
   :alt: Editing mode

.. image:: /_static/screenshots/darkmode/editing_mode.png
   :class: only-dark
   :alt: Editing mode

Transform geometries interactively using translate, rotate, and scale controls.

**Entering Editing Mode:**

1. Press ``E`` to enter editing mode (from view mode)
2. Select geometry instances by clicking on them
3. Use the transform gizmo to manipulate selected objects
4. Press ``T`` to cycle between translate, rotate, and scale modes
5. Hold ``X``, ``Y``, or ``Z`` to constrain movement to a single axis
6. Press ``S`` to save changes
7. Press ``E`` again to exit and return to view mode

.. image:: /_static/screenshots/lightmode/editing_axis_constraint.png
   :class: only-light
   :alt: Axis constraint indicator

.. image:: /_static/screenshots/darkmode/editing_axis_constraint.png
   :class: only-dark
   :alt: Axis constraint indicator

When holding an axis key, a colored chip indicates the active constraint.

**Editing Curves:**

.. image:: /_static/screenshots/lightmode/curve_editing.png
   :class: only-light
   :alt: Curve editing with virtual markers

.. image:: /_static/screenshots/darkmode/curve_editing.png
   :class: only-dark
   :alt: Curve editing with virtual markers

In editing mode, curves display virtual markers between control points.
Click a virtual marker to insert a new control point at that position.
Use ``Delete`` or ``Backspace`` to remove selected markers.

.. note::

   See :doc:`keyboard-shortcuts` for all controls.


Dynamic Properties
^^^^^^^^^^^^^^^^^^

.. image:: /_static/screenshots/lightmode/dynamic_properties_dropdown.png
   :class: only-light
   :alt: Dynamic properties dropdown

.. image:: /_static/screenshots/darkmode/dynamic_properties_dropdown.png
   :class: only-dark
   :alt: Dynamic properties dropdown

Geometry properties like ``position`` and ``direction`` can reference atom data dynamically.
Instead of specifying fixed coordinates, use string references to atom arrays:

.. code:: python

    import numpy as np
    from zndraw.geometries import Arrow

    # Create atoms with calculated forces
    atoms = ase.Atoms("H4", positions=[(0, 0, 0), (2, 0, 0), (0, 2, 0), (2, 2, 0)])
    atoms.arrays["forces"] = np.array([
        [0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0]
    ], dtype=float)
    vis.append(atoms)

    # Create arrows showing forces at each atom position
    vis.geometries["force_arrows"] = Arrow(
        position="arrays.positions",  # Reference atom positions
        direction="arrays.forces",     # Reference force vectors
        color=["#ff6600"],
        radius=0.1,
    )

Available dynamic property references are computed from the ``atoms.info``, ``atoms.arrays``, and if available, ``atoms.calc.results`` dictionaries:

- ``arrays.positions`` - Atom positions
- ``arrays.numbers`` - Atomic numbers
- ``arrays.colors`` - Per-atom colors
- ``arrays.radii`` - Per-atom radii
- ``arrays.forces`` - Calculated forces (if available)
- ``calc.energy`` - Calculated energy
- ``info.connectivity`` - Bond connectivity


Analysis & Figures
------------------

.. image:: /_static/screenshots/lightmode/analysis_1d.png
   :class: only-light
   :alt: 1D analysis plot

.. image:: /_static/screenshots/darkmode/analysis_1d.png
   :class: only-dark
   :alt: 1D analysis plot

Display interactive Plotly figures with ``vis.figures``:

.. code:: python

    import plotly.express as px
    import pandas as pd

    # Create figure from data
    df = pd.DataFrame({
        "frame": range(len(vis)),
        "energy": [atoms.get_potential_energy() for atoms in vis]
    })

    fig = px.line(df, x="frame", y="energy", title="Energy vs Frame")

    # Display in ZnDraw
    vis.figures["energy_plot"] = fig

    # Remove figure
    del vis.figures["energy_plot"]

.. image:: /_static/screenshots/lightmode/analysis_2d.png
   :class: only-light
   :alt: 2D analysis scatter plot

.. image:: /_static/screenshots/darkmode/analysis_2d.png
   :class: only-dark
   :alt: 2D analysis scatter plot

2D analysis with scatter plots:

.. code:: python

    # 2D scatter plot
    df = pd.DataFrame({
        "ml_energy": [...],
        "dft_energy": [...]
    })

    fig = px.scatter(df, x="ml_energy", y="dft_energy", title="ML vs DFT Energy")

    vis.figures["comparison"] = fig


Molecule Building
-----------------

.. image:: /_static/screenshots/lightmode/molecule_builder.png
   :class: only-light
   :alt: Molecule builder

.. image:: /_static/screenshots/darkmode/molecule_builder.png
   :class: only-dark
   :alt: Molecule builder

Build molecules from SMILES strings using the molecule builder:

- Add molecules from SMILES notation
- Use the Ketcher molecular editor
- Pack molecules into simulation boxes

**Ketcher Editor:**

.. image:: /_static/screenshots/lightmode/molecule_builder_editor.png
   :class: only-light
   :alt: Ketcher molecular editor

.. image:: /_static/screenshots/darkmode/molecule_builder_editor.png
   :class: only-dark
   :alt: Ketcher molecular editor

.. note::

   The Ketcher editor currently does not support dark mode. See `Ketcher issue #5353 <https://github.com/epam/ketcher/issues/5353>`_ for more information.


Chat & Logging
--------------

.. image:: /_static/screenshots/lightmode/chat.png
   :class: only-light
   :alt: Chat panel

.. image:: /_static/screenshots/darkmode/chat.png
   :class: only-dark
   :alt: Chat panel

Send messages to the chat panel:

.. code:: python

    # Send a message
    vis.log("Analysis complete!")

    # Messages support Markdown and LaTeX
    vis.log("Energy: $E = mc^2$")

    # Get chat history
    messages = vis.get_messages(limit=10)


Frame References
^^^^^^^^^^^^^^^^

.. image:: /_static/screenshots/lightmode/chat_frame_reference.png
   :class: only-light
   :alt: Chat with frame references

.. image:: /_static/screenshots/darkmode/chat_frame_reference.png
   :class: only-dark
   :alt: Chat with frame references

Reference frames in chat messages using ``@{frame}`` syntax. Frame references
become clickable chips that navigate to the referenced frame:

.. code:: python

    # Reference specific frames in messages
    vis.log("Initial structure at @0")
    vis.log("Compare @10 with @15 to see the transition")

Clicking a frame reference chip navigates directly to that frame.


Markdown & Code Blocks
^^^^^^^^^^^^^^^^^^^^^^

Chat messages support full Markdown rendering including:

- **Text formatting**: ``**bold**``, ``*italic*``, ``~~strikethrough~~``
- **Lists**: Ordered and unordered lists
- **Links**: ``[text](url)``
- **LaTeX math**: Inline ``$E = mc^2$`` or block ``$$\sum_{i=1}^n x_i$$``
- **Code blocks**: Syntax-highlighted code with language specification

.. code:: python

    vis.log("""
    ## Results Summary

    The simulation converged after **1000 steps**.

    Energy: $E = -42.5$ eV

    ```python
    for atom in atoms:
        print(atom.symbol)
    ```
    """)


Progress Bars
^^^^^^^^^^^^^

.. image:: /_static/screenshots/lightmode/chat_progress.png
   :class: only-light
   :alt: Chat with progress bars

.. image:: /_static/screenshots/darkmode/chat_progress.png
   :class: only-dark
   :alt: Chat with progress bars

Display progress bars in chat using the ``progress`` code block syntax:

.. code:: python

    vis.log("""
    ```progress
    description: Processing frames
    value: 75
    max: 100
    color: success
    ```
    """)

**Parameters:**

- ``value``: Current progress value. If omitted, shows an indeterminate spinner.
- ``min``: Minimum value (default: ``0``)
- ``max``: Maximum value (default: ``100``)
- ``description``: Label text displayed above the progress bar
- ``color``: MUI color - ``primary``, ``secondary``, ``success``, ``error``, ``warning``, ``info`` (default: ``primary``)

The progress bar displays the percentage and the current value relative to max.
For long-running operations, consider using the :meth:`~zndraw.ZnDraw.progress_tracker`
context manager instead, which provides real-time updates.


Molecule Structures
^^^^^^^^^^^^^^^^^^^

.. image:: /_static/screenshots/lightmode/chat_smiles.png
   :class: only-light
   :alt: Chat with SMILES molecule rendering

.. image:: /_static/screenshots/darkmode/chat_smiles.png
   :class: only-dark
   :alt: Chat with SMILES molecule rendering

Display molecule structures in chat using SMILES notation with the ``smiles`` code block syntax:

.. code:: python

    vis.log("""
    ```smiles
    CCO
    ```
    """)

The SMILES string is rendered as a 2D molecule structure image using RDKit.


Property Inspector
------------------

.. image:: /_static/screenshots/lightmode/property_inspector.png
   :class: only-light
   :alt: Property inspector info boxes

.. image:: /_static/screenshots/darkmode/property_inspector.png
   :class: only-dark
   :alt: Property inspector info boxes

The property inspector displays frame properties in floating info boxes.
Press ``i`` to toggle visibility. Configure which properties to display
via settings:

.. code:: python

    # Enable properties in the inspector
    vis.sessions["<sessionId>"].settings.property_inspector.enabled_properties = [
        "calc.energy",
        "calc.forces",
    ]

Two info boxes are available:

- **Scene Info** (top-right): Displays global properties like ``calc.energy``
- **Hover Info** (follows cursor): Shows per-particle properties when hovering over atoms

Properties are automatically categorized based on their shape:

- **Global**: Scalar values or arrays not matching particle count
- **Per-particle**: Arrays with first dimension equal to particle count


Browser Sessions
----------------

Access connected browser windows via ``vis.sessions``. Each frontend session
has its own camera and rendering settings:

.. code:: python

    # List all connected browser sessions
    for session_id in vis.sessions:
        print(session_id)

    # Access a specific session
    session = vis.sessions["abc-123"]

    # Get/set camera for that browser window
    cam = session.camera
    print(cam.position, cam.target)

    # Update camera position
    from zndraw.geometries import Camera
    session.camera = Camera(position=(10, 5, 10), target=(0, 0, 0), fov=60)

    # Access session settings
    settings = session.settings
    settings.studio_lighting.key_light = 1.5  # adjust settings

.. note::

    Only frontend browser windows appear in ``vis.sessions``.
    Python API clients do not create entries here.


Progress Tracking
-----------------

.. image:: /_static/screenshots/lightmode/progress_tracker.png
   :class: only-light
   :alt: Progress tracker

.. image:: /_static/screenshots/darkmode/progress_tracker.png
   :class: only-dark
   :alt: Progress tracker

Track long-running operations with ``vis.progress_tracker()``:

.. code:: python

    with vis.progress_tracker("Processing data") as tracker:
        for i, item in enumerate(items):
            process(item)
            tracker.update(
                f"Step {i + 1}/{len(items)}",
                progress=(i + 1) / len(items) * 100  # 0-100 percentage
            )

The progress bar appears in the UI with the current message and completion percentage.


Lock Mechanism
--------------

.. image:: /_static/screenshots/lightmode/locked_room.png
   :class: only-light
   :alt: Locked room indicator

.. image:: /_static/screenshots/darkmode/locked_room.png
   :class: only-dark
   :alt: Locked room indicator

Use ``vis.get_lock()`` for safe batch operations that prevent concurrent modifications:

.. code:: python

    # Lock the room during bulk operations
    with vis.get_lock(msg="Uploading trajectory..."):
        for atoms in trajectory:
            vis.append(atoms)

    # Lock specific targets for fine-grained control
    with vis.get_lock(target="step"):
        vis.step = 42

While a lock is held, other clients see a locked indicator (shown above) and cannot
modify the locked resources. The lock message is displayed in the UI so users know
what operation is in progress. This is useful when uploading large trajectories or
performing multi-step operations that should not be interrupted.


Custom Extensions
-----------------

.. image:: /_static/screenshots/lightmode/custom_modifier.png
   :class: only-light
   :alt: Custom modifier interface

.. image:: /_static/screenshots/darkmode/custom_modifier.png
   :class: only-dark
   :alt: Custom modifier interface

ZnDraw supports custom extensions for modifiers, selections, and analysis.
Register your own Python classes to extend the UI:

.. code:: python

    from pydantic import Field
    from zndraw.extensions import Extension, Category

    class ScaleAtoms(Extension):
        """Scale atom positions by a factor."""

        category = Category.MODIFIER  # or SELECTION or ANALYSIS
        factor: float = Field(
            1.5, ge=0.1, le=5.0,
            description="Scale factor",
            json_schema_extra={"format": "range"},
        )
        center_first: bool = Field(
            True,
            description="Center atoms before scaling",
        )

        def run(self, vis, **kwargs):
            atoms = vis.atoms.copy()
            if self.center_first:
                atoms.positions -= atoms.get_center_of_mass()
            atoms.positions *= self.factor
            vis.append(atoms)
            vis.step = len(vis) - 1

    # Register the extension
    vis.register_extension(ScaleAtoms)


Extension Categories
^^^^^^^^^^^^^^^^^^^^

Extensions are categorized by their purpose:

- ``Category.MODIFIER``: Modify atomic structures (e.g., delete, rotate, translate)
- ``Category.SELECTION``: Select atoms (e.g., by type, neighbors, random)
- ``Category.ANALYSIS``: Analyze data and create plots (e.g., properties, correlations)


Schema Customization
^^^^^^^^^^^^^^^^^^^^

Use ``json_schema_extra`` and ``Field`` options to customize how fields appear in the UI:

**Slider Input**

.. code:: python

    factor: float = Field(
        1.0, ge=0.0, le=10.0,
        json_schema_extra={"format": "range"},
    )

**Dynamic Dropdowns**

Populate dropdowns at runtime from available data:

.. code:: python

    # Dropdown from geometry names (filtered to Curves)
    curve: str = Field(
        "curve",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-geometries"],
            "x-geometry-filter": "Curve",
        },
    )

    # Dropdown from atom/frame properties
    property: str = Field(
        ...,
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props"],
        },
    )

**SMILES Input with Ketcher Editor**

.. code:: python

    smiles: str = Field(
        ...,
        json_schema_extra={"x-custom-type": "smiles"},
    )

**Available Options**

+----------------------------------+-------------------------------------------+
| Option                           | Description                               |
+==================================+===========================================+
| ``"format": "range"``            | Render as slider (requires ge/le bounds)  |
+----------------------------------+-------------------------------------------+
| ``"x-custom-type": "smiles"``    | SMILES input with Ketcher editor button   |
+----------------------------------+-------------------------------------------+
| ``"x-custom-type":               | Runtime-populated dropdown                |
| "dynamic-enum"``                 |                                           |
+----------------------------------+-------------------------------------------+
| ``"x-features":                  | Populate from geometry names              |
| ["dynamic-geometries"]``         |                                           |
+----------------------------------+-------------------------------------------+
| ``"x-features":                  | Populate from atom/frame properties       |
| ["dynamic-atom-props"]``         |                                           |
+----------------------------------+-------------------------------------------+
| ``"x-geometry-filter":           | Filter geometries by type                 |
| "Curve"``                        |                                           |
+----------------------------------+-------------------------------------------+


Custom Filesystems
------------------

Register any `fsspec <https://filesystem-spec.readthedocs.io/>`_-compatible
filesystem to browse and load files from the UI:

.. code:: python

    from fsspec.implementations.dirfs import DirFileSystem

    vis.register_filesystem(DirFileSystem(path="."), "local")

All fsspec-compatible filesystems are supported, including S3, GCS, Azure, HDFS, and more.
