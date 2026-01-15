"""Screenshot tests for ZnDraw documentation."""

import ase

from zndraw import ZnDraw
from zndraw.geometries import Arrow, Box, Camera, Curve, Floor, Sphere
from zndraw.transformations import CurveAttachment


def test_overview(server, page, capture, bmim_bf4, request):
    """Capture the main overview screenshot."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_selection(server, page, capture, bmim_bf4, request):
    """Capture screenshot with selection sidebar open."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)
    vis.selection = list(range(100))
    vis.selection_groups["bf4"] = {"particles": list(range(32 * 5))}
    vis.selection_groups["bmim"] = {"particles": list(range(32 * 25))}

    page.goto(f"{server}/room/{request.node.name}")

    page.get_by_role("button", name="Selection tools and groups").click()
    # blur tooltip, by clicking at (0, 0)
    page.mouse.click(0, 0)
    page.wait_for_timeout(2500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_drawing_mode(server, page, capture, request):
    """Capture drawing mode with geometries."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(ase.Atoms())
    vis.geometries["floor"] = Floor(active=True, height=-5.0, color="#A0A0A0")
    vis.geometries["box"] = Box(position=[(0, 0, 0)], size=[(10, 10, 10)])
    vis.geometries["curve"] = Curve(
        position=[
            (-4.789177573713647, 5.474554893697157, 5.646802136489912),
            (-4.272226349849859, -4.421872431387516, 5.1448314636945724),
            (5.55918672834876, -4.537280536828369, 5.1448314636945724),
            (5.58483963937306, 5.279062292631962, 5.1448314636945724),
            (5.651912331581116, 5.290737777511481, -4.622419531622475),
            (5.651912331581116, -4.54817732764105, -4.742356696752667),
            (12.155146865890176, 1.703348269955838, 1.9575992621569651),
        ]
    )
    vis.selections["curve"] = [0]

    page.goto(f"{server}/room/{request.node.name}")

    page.wait_for_timeout(2500)

    page.keyboard.press("e")

    capture.light()
    capture.toggle()
    page.keyboard.press("e")
    capture.dark()


def test_geometries(server, page, capture, request):
    """Capture geometries (Box, Sphere, Curve, Arrow) with different materials."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(ase.Atoms())

    # Remove default geometries
    for key in list(vis.geometries.keys()):
        del vis.geometries[key]

    # Add floor for context
    vis.geometries["floor"] = Floor(active=True, height=-2.0, color="#808080")

    # Add a red box with toon/cartoon material
    vis.geometries["box"] = Box(
        position=[(0, 2, 0)],
        size=[(4, 4, 4)],
        color=["#e74c3c"],
        material="MeshToonMaterial",
    )

    # Add a blue sphere with glass material
    vis.geometries["sphere"] = Sphere(
        position=[(8, 2, 0)],
        radius=[2.0],
        color=["#3498db"],
        material="MeshPhysicalMaterial_glass",
    )

    # Add a green curve
    vis.geometries["curve"] = Curve(
        position=[
            (-6, 0, -6),
            (-3, 4, -3),
            (0, 0, 0),
            (3, 4, 3),
            (6, 0, 6),
        ],
        color="#2ecc71",
    )

    # Add an orange arrow with shiny material
    vis.geometries["arrow"] = Arrow(
        position=[(12, 0, 0)],
        direction=[(0, 5, 0)],
        color=["#f39c12"],
        material="MeshPhysicalMaterial_shiny",
    )

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(2000)

    capture.light()
    capture.toggle()
    capture.dark()


def test_geometry_viewer(server, page, capture, bmim_bf4, request):
    """Capture the geometry management panel."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # Add a floor and curve to show different geometry types in panel
    vis.geometries["floor"] = Floor(active=True, height=-5.0, color="#808080")
    vis.geometries["my_curve"] = Curve(
        position=[
            (-10, 0, -10),
            (0, 10, 0),
            (10, 0, 10),
        ],
        color="#2ecc71",
    )

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(2500)

    page.get_by_role("button", name="Manage geometries").click()

    capture.light()
    capture.toggle()
    capture.dark()


def test_bookmarks(server, page, capture, bmim_bf4, request):
    """Capture bookmarks on the timeline."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.extend([bmim_bf4 for _ in range(20)])
    vis.bookmarks[5] = "Frame 5"
    vis.bookmarks[17] = "Frame 17"
    vis.step = 5
    page.goto(f"{server}/room/{request.node.name}")

    # Hover and capture light mode
    page.get_by_role("button", name="Bookmark: Frame 5 (Frame 5)").hover()
    page.wait_for_timeout(400)
    capture.light()

    # Toggle, re-hover, capture dark mode
    capture.toggle()
    page.get_by_role("button", name="Bookmark: Frame 5 (Frame 5)").hover()
    page.wait_for_timeout(400)
    capture.dark()

    # Reset to light
    capture.toggle()


def test_timeline(server, page, capture, bmim_bf4, request):
    """Capture the timeline with multiple frames and playback controls."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.extend([bmim_bf4 for _ in range(50)])
    vis.step = 25  # Set to middle frame

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1000)

    capture.light()
    capture.toggle()
    capture.dark()


def test_analysis_1d(server, page, capture, bmim_bf4_e_f, request):
    """Capture Properties1D analysis screenshot."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.extend(bmim_bf4_e_f)
    vis.step = 16

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1000)

    # Open analysis panel
    page.get_by_role("button", name="Analysis tools").click()
    page.wait_for_timeout(500)

    # Select Properties1D from the dropdown
    page.get_by_label("analysis Method").click()
    page.get_by_role("option", name="Properties1D").click()
    page.wait_for_timeout(500)

    # Select calc.energy from the value dropdown
    page.get_by_label("Value").click()
    page.get_by_role("option", name="calc.energy").click()
    page.wait_for_timeout(500)

    # Run the extension
    page.get_by_role("button", name="Run Extension").click()
    page.wait_for_timeout(2000)

    capture.light()
    capture.toggle()
    capture.dark()


def test_analysis_2d(server, page, capture, bmim_bf4_e_f, request):
    """Capture Properties2D analysis screenshot."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.extend(bmim_bf4_e_f)
    vis.step = 16

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1000)

    # Open analysis panel
    page.get_by_role("button", name="Analysis tools").click()
    page.wait_for_timeout(500)

    # Select Properties1D from the dropdown
    page.get_by_label("analysis Method").click()
    page.get_by_role("option", name="Properties1D").click()
    page.wait_for_timeout(500)

    # Select calc.energy from the value dropdown
    page.get_by_label("Value").click()
    page.get_by_role("option", name="calc.energy").click()
    page.wait_for_timeout(500)

    # Run the extension
    page.get_by_role("button", name="Run Extension").click()
    page.wait_for_timeout(2000)

    # Select Properties2D from the dropdown
    page.get_by_label("analysis Method").click()
    page.get_by_role("option", name="Properties2D").click()
    page.wait_for_timeout(500)

    # Set x_data to calc.energy
    page.get_by_label("X data").click()
    page.get_by_role("option", name="calc.energy").click()
    page.wait_for_timeout(300)

    # Set y_data to calc.dft_energy
    page.get_by_label("Y data").click()
    page.get_by_role("option", name="calc.dft_energy").click()
    page.wait_for_timeout(300)

    # Set color to calc.energy
    page.get_by_label("Color").click()
    page.get_by_role("option", name="calc.energy").click()
    page.wait_for_timeout(300)

    # Run the extension
    page.get_by_role("button", name="Run Extension").click()
    page.wait_for_timeout(2000)

    capture.light()
    capture.toggle()
    capture.dark()


def test_molecule_builder(server, page, capture, request):
    """Capture AddFromSMILES interface with amino acid SMILES."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(ase.Atoms())

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(500)

    # Open modifiers panel
    page.get_by_role("button", name="Modifier tools").click()
    page.wait_for_timeout(500)

    # Select AddFromSMILES
    page.get_by_label("modifiers Method").click()
    page.get_by_role("option", name="AddFromSMILES").click()
    page.wait_for_timeout(500)

    # Enter Phenylalanine SMILES
    page.get_by_placeholder("Enter SMILES notation").fill("NC(Cc1ccccc1)C(=O)O")
    page.wait_for_timeout(500)

    # Run extension to add the molecule
    page.get_by_role("button", name="Run Extension").click()
    page.wait_for_timeout(2000)

    capture.light()
    capture.toggle()
    page.wait_for_timeout(500)
    capture.dark()


def test_molecule_builder_editor(server, page, capture, request):
    """Capture Ketcher molecular editor with amino acid."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(ase.Atoms())

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(500)

    # Open modifiers panel
    page.get_by_role("button", name="Modifier tools").click()
    page.wait_for_timeout(500)

    # Select AddFromSMILES
    page.get_by_label("modifiers Method").click()
    page.get_by_role("option", name="AddFromSMILES").click()
    page.wait_for_timeout(500)

    # Enter Phenylalanine SMILES and run extension first
    page.get_by_placeholder("Enter SMILES notation").fill("NC(Cc1ccccc1)C(=O)O")
    page.wait_for_timeout(500)

    # Run extension to add the molecule
    page.get_by_role("button", name="Run Extension").click()
    page.wait_for_timeout(2000)

    # Enter Alanine SMILES for the editor view
    page.get_by_placeholder("Enter SMILES notation").fill("CC(N)C(=O)O")
    page.wait_for_timeout(500)

    # Open Ketcher editor (lazy loaded)
    page.get_by_role("button", name="Draw", exact=True).click()
    # Wait for the Ketcher dialog to appear
    page.get_by_role("dialog", name="Molecular Structure Editor").wait_for(
        state="visible"
    )
    page.wait_for_timeout(1000)  # Extra delay for molecule canvas to render

    capture.light()
    capture.toggle()

    # Reopen Ketcher editor for dark mode (dialog closes on theme toggle)
    page.get_by_role("button", name="Draw", exact=True).click()
    page.get_by_role("dialog", name="Molecular Structure Editor").wait_for(
        state="visible"
    )
    page.wait_for_timeout(1000)

    capture.dark()


def test_rooms(server, page, capture, bmim_bf4):
    """Capture the rooms management page with multiple rooms."""
    # Create multiple rooms with different frame counts
    vis1 = ZnDraw(url=server, room="demo-simulation")
    vis1.extend([bmim_bf4 for _ in range(50)])

    vis2 = ZnDraw(url=server, room="analysis-results")
    vis2.extend([bmim_bf4 for _ in range(25)])

    vis3 = ZnDraw(url=server, room="molecule-library")
    vis3.extend([bmim_bf4 for _ in range(10)])

    # First visit a room to access theme toggle
    page.goto(f"{server}/room/demo-simulation")
    page.wait_for_timeout(500)

    # Navigate to rooms page for light mode capture
    page.goto(f"{server}/rooms")
    page.wait_for_timeout(1500)
    capture.light()

    # Go back to room to toggle theme
    page.goto(f"{server}/room/demo-simulation")
    page.wait_for_timeout(500)
    capture.toggle()

    # Navigate back to rooms for dark mode capture
    page.goto(f"{server}/rooms")
    page.wait_for_timeout(1500)
    capture.dark()


def test_file_browser(server, page, capture, bmim_bf4, request):
    """Capture the file browser page."""
    # Create a room first to access theme toggle
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # First visit a room to access theme toggle
    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(500)

    # Navigate to file browser for light mode capture
    page.goto(f"{server}/file-browser")
    page.wait_for_timeout(1500)
    capture.light()

    # Go back to room to toggle theme
    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(500)
    capture.toggle()

    # Navigate back to file browser for dark mode capture
    page.goto(f"{server}/file-browser")
    page.wait_for_timeout(1500)
    capture.dark()


def test_camera(server, page, capture, bmim_bf4, request):
    """Capture camera with curve attachment from Python."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # Create curve for camera position path (flight path around molecule)
    vis.geometries["cam_pos"] = Curve(
        position=[
            (25, 15, 25),
            (30, 20, 15),
            (25, 15, 5),
        ],
        color="#3498db",
    )

    # Create curve for camera target path (looking at center of molecule)
    vis.geometries["cam_target"] = Curve(
        position=[
            (10, 10, 10),
            (12, 10, 10),
            (10, 10, 10),
        ],
        color="#e74c3c",
    )

    # Create camera attached to curves
    vis.geometries["camera"] = Camera(
        position=CurveAttachment(geometry_key="cam_pos", progress=0.5),
        target=CurveAttachment(geometry_key="cam_target", progress=0.5),
        fov=60,
        helper_visible=True,
        helper_color="#00ff00",
    )

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1500)

    # Open geometry viewer to show camera settings
    page.get_by_role("button", name="Manage geometries").click()
    page.wait_for_timeout(1000)

    capture.light()
    capture.toggle()
    capture.dark()


def test_python_connection(server, page, capture, bmim_bf4, request):
    """Capture Python code connection info dialog."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(500)

    # Open connection info dialog and capture light mode
    page.get_by_role("button", name="show connection info").click()
    page.wait_for_timeout(500)
    capture.light()

    # Close dialog, toggle to dark mode, reopen and capture
    page.keyboard.press("Escape")
    page.wait_for_timeout(200)
    capture.toggle()
    page.get_by_role("button", name="show connection info").click()
    page.wait_for_timeout(500)
    capture.dark()


def test_chat(server, page, capture, bmim_bf4, request):
    """Capture chat panel with code examples and formulas."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # Add chat messages with different content types
    vis.log("Welcome to ZnDraw! This chat supports **markdown** formatting.")

    vis.log('Here\'s a code example:\n\n```py\nprint("Hello ZnDraw!")\n```')

    vis.log(
        "The potential energy is given by:\n\n"
        "$$E = \\sum_{i<j} 4\\varepsilon \\left[ "
        "\\left(\\frac{\\sigma}{r_{ij}}\\right)^{12} - "
        "\\left(\\frac{\\sigma}{r_{ij}}\\right)^{6} \\right]$$"
    )

    vis.log("You can also use inline math like $E = mc^2$ in your messages.")

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1000)

    # Open chat panel
    page.get_by_role("button", name="toggle chat").click()
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_custom_modifier(server, page, capture, bmim_bf4, request):
    """Capture custom modifier registered with vis.register_extension()."""
    from pydantic import Field

    from zndraw.extensions import Category, Extension

    class ScaleAtoms(Extension):
        """Scale atom positions by a factor."""

        category = Category.MODIFIER
        factor: float = Field(
            1.5,
            ge=0.1,
            le=5.0,
            description="Scale factor",
            json_schema_extra={"format": "range"},
        )
        center_first: bool = Field(
            True,
            description="Center atoms before scaling",
            json_schema_extra={"format": "checkbox"},
        )

        def run(self, vis, **kwargs):
            atoms = vis.atoms.copy()
            if self.center_first:
                atoms.positions -= atoms.get_center_of_mass()
            atoms.positions *= self.factor
            vis.append(atoms)
            vis.step = len(vis) - 1
            vis.log(f"Scaled atoms by factor {self.factor}")

    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)
    vis.register_extension(ScaleAtoms)

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1000)

    # Open modifiers panel
    page.get_by_role("button", name="Modifier tools").click()
    page.wait_for_timeout(500)

    # Select ScaleAtoms from the dropdown
    page.get_by_label("modifiers Method").click()
    page.get_by_role("option", name="ScaleAtoms").click()
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_progress_tracker(server, page, capture, bmim_bf4, request):
    """Capture progress tracker showing active progress state."""
    import threading
    import time

    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # Hold progress tracker active in background thread
    progress_active = threading.Event()
    progress_done = threading.Event()

    def show_progress():
        with vis.progress_tracker("Geometry Optimization") as tracker:
            tracker.update("Step 42/100 - Energy: -1234.56 eV", progress=42)
            progress_active.set()
            # Hold progress visible until screenshots taken
            for _ in range(100):
                if progress_done.is_set():
                    break
                time.sleep(0.1)

    progress_thread = threading.Thread(target=show_progress, daemon=True)
    progress_thread.start()

    # Wait for progress to be set
    progress_active.wait(timeout=5)
    time.sleep(0.3)

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(2000)

    capture.light()
    capture.toggle()
    capture.dark()

    # Release progress tracker
    progress_done.set()
    progress_thread.join(timeout=2)


def test_locked_room(server, page, capture, bmim_bf4, request):
    """Capture a room with an active lock showing the lock indicator."""
    import threading
    import time

    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # Start a background thread that holds a lock
    lock_active = threading.Event()
    lock_done = threading.Event()

    def hold_lock():
        with vis.get_lock(msg="Uploading trajectory data..."):
            lock_active.set()
            # Hold lock for 10 seconds or until signaled
            for _ in range(100):
                if lock_done.is_set():
                    break
                time.sleep(0.1)

    lock_thread = threading.Thread(target=hold_lock, daemon=True)
    lock_thread.start()

    # Wait for lock to be acquired
    lock_active.wait(timeout=5)
    time.sleep(0.5)  # Give time for UI to update

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(2000)

    capture.light()
    capture.toggle()
    capture.dark()

    # Release the lock
    lock_done.set()
    lock_thread.join(timeout=2)


def test_dynamic_properties_dropdown(server, page, capture, request):
    """Capture the dynamic properties dropdown for Arrow geometry."""
    import numpy as np

    vis = ZnDraw(url=server, room=request.node.name)

    # Create atoms with forces for demonstration
    atoms = ase.Atoms(
        "H4",
        positions=[(0, 0, 0), (2, 0, 0), (0, 2, 0), (2, 2, 0)],
    )
    # Add calculated forces to demonstrate dynamic properties
    atoms.arrays["forces"] = np.array(
        [
            [0, 0, 1],
            [0, 0, -1],
            [1, 0, 0],
            [-1, 0, 0],
        ],
        dtype=float,
    )
    vis.append(atoms)

    # Add arrow geometry using dynamic properties
    vis.geometries["force_arrows"] = Arrow(
        position="arrays.positions",
        direction="arrays.forces",
        color=["#ff6600"],
        radius=0.1,
    )

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(2000)

    # Open geometry manager
    page.get_by_role("button", name="Manage geometries").click()
    page.wait_for_timeout(500)

    # Click on the row with force_arrows to open edit mode (MUI DataGrid row click)
    page.get_by_role("row", name="force_arrows").click()
    page.wait_for_timeout(500)

    # Click on the position dropdown to show dynamic options
    page.get_by_label("Position").click()
    page.wait_for_timeout(500)

    capture.light()

    # Close dropdown, toggle theme, then re-open dropdown
    page.keyboard.press("Escape")
    page.wait_for_timeout(200)
    capture.toggle()
    page.get_by_label("Position").click()
    page.wait_for_timeout(500)

    capture.dark()


def test_property_inspector(server, page, capture, bmim_bf4_e_f, request):
    """Capture property inspector settings panel with info boxes visible."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.extend(bmim_bf4_e_f[:10])

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(2000)

    # Press 'i' to show info boxes
    page.keyboard.press("i")
    page.wait_for_timeout(500)

    # Open settings panel
    page.get_by_role("button", name="Application settings").click()
    page.wait_for_timeout(500)

    # Select "Property Inspector" from dropdown
    page.get_by_label("Settings Category").click()
    page.wait_for_timeout(300)
    page.get_by_role("option", name="Property Inspector").click()
    page.wait_for_timeout(500)

    # Select arrays.positions (per-particle property) by clicking its list item
    page.get_by_text("arrays.positions").click()
    page.wait_for_timeout(300)

    # Select calc.energy (global property) by clicking its list item
    page.get_by_text("calc.energy").click()
    page.wait_for_timeout(500)

    # Hover over a particle - particles are on right side around x=950, y=280
    page.mouse.move(950, 280)
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()

    # Re-hover after theme toggle to ensure HoverInfoBox is visible
    page.mouse.move(950, 280)
    page.wait_for_timeout(500)
    capture.dark()


def test_chat_frame_reference(server, page, capture, bmim_bf4, request):
    """Capture chat with clickable @frame references."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.extend([bmim_bf4 for _ in range(20)])

    # Send messages with frame references that become clickable chips
    vis.log("Initial structure loaded at @0")
    vis.log("Check the transition at @10 and compare with @15")

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1000)

    # Open chat panel
    page.get_by_role("button", name="toggle chat").click()
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_chat_progress(server, page, capture, bmim_bf4, request):
    """Capture chat with progress bar examples."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # Show determinate progress bar
    vis.log(
        "Processing simulation data:\n\n"
        "```progress\n"
        "description: Analyzing frames\n"
        "value: 75\n"
        "max: 100\n"
        "color: success\n"
        "```"
    )

    # Show indeterminate progress spinner
    vis.log(
        "Waiting for calculation:\n\n"
        "```progress\n"
        "description: Optimizing geometry...\n"
        "color: primary\n"
        "```"
    )

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1000)

    page.get_by_role("button", name="toggle chat").click()
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_editing_mode(server, page, capture, request):
    """Capture editing mode with transform controls and editing indicator."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(ase.Atoms())

    # Add geometries to edit
    vis.geometries["box"] = Box(
        position=[(0, 2, 0)],
        size=[(4, 4, 4)],
        color=["#3498db"],
    )
    vis.geometries["sphere"] = Sphere(
        position=[(8, 2, 0)],
        radius=[2.0],
        color=["#e74c3c"],
    )

    # Select both geometries to show transform controls
    vis.selections["box"] = [0]
    vis.selections["sphere"] = [0]

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1500)

    # Enter editing mode
    page.keyboard.press("e")
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_editing_axis_constraint(server, page, capture, request):
    """Capture editing mode with axis constraint indicator."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(ase.Atoms())

    vis.geometries["box"] = Box(
        position=[(0, 2, 0)],
        size=[(4, 4, 4)],
        color=["#3498db"],
    )

    # Select the box to show transform controls
    vis.selections["box"] = [0]

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1500)

    # Enter editing mode
    page.keyboard.press("e")
    page.wait_for_timeout(500)

    # Hold X key to show axis constraint (using keydown without keyup)
    page.keyboard.down("x")
    page.wait_for_timeout(300)

    capture.light()

    page.keyboard.up("x")
    capture.toggle()

    # Re-hold X for dark mode capture
    page.keyboard.down("x")
    page.wait_for_timeout(300)
    capture.dark()
    page.keyboard.up("x")


def test_curve_editing(server, page, capture, request):
    """Capture curve in editing mode showing virtual markers."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(ase.Atoms())

    # Add a curve with visible markers
    vis.geometries["curve"] = Curve(
        position=[
            (-6, 0, -6),
            (-3, 4, -3),
            (0, 0, 0),
            (3, 4, 3),
            (6, 0, 6),
        ],
        color="#2ecc71",
        marker={"enabled": True, "size": 0.15, "opacity": 1.0},
        virtual_marker={"enabled": True, "size": 0.1, "opacity": 0.6},
    )

    # Select the two peak control points to show transform controls
    vis.selections["curve"] = [1, 3]

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1500)

    # Enter editing mode to show virtual markers
    page.keyboard.press("e")
    page.wait_for_timeout(500)

    capture.light()
    capture.toggle()
    capture.dark()


def test_chat_smiles(server, page, capture, bmim_bf4, request):
    """Capture chat with SMILES molecule rendering."""
    vis = ZnDraw(url=server, room=request.node.name)
    vis.append(bmim_bf4)

    # Send messages with SMILES code blocks
    vis.log("Here's ethanol:\n\n```smiles\nCCO\n```")
    vis.log("And here's caffeine:\n\n```smiles\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C\n```")

    page.goto(f"{server}/room/{request.node.name}")
    page.wait_for_timeout(1500)

    # Open chat panel
    page.get_by_role("button", name="toggle chat").click()
    page.wait_for_timeout(1000)

    capture.light()
    capture.toggle()
    capture.dark()
