"""Pytest fixtures for screenshot generation."""

import subprocess
import time
from pathlib import Path

import ase
import molify
import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator
from playwright.sync_api import Page, sync_playwright

MISC_DIR = Path(__file__).parent
SCREENSHOTS_DIR = MISC_DIR.parent / "docs" / "source" / "_static" / "screenshots"
PORT = 5000
BASE_URL = f"http://localhost:{PORT}"


def _wait_for_server(url: str, timeout: float = 30.0) -> bool:
    """Wait for the server to be ready."""
    import urllib.request

    start = time.time()
    while time.time() - start < timeout:
        try:
            urllib.request.urlopen(url, timeout=1)
            return True
        except Exception:
            time.sleep(0.5)
    return False


def _toggle_color_mode(page: Page) -> None:
    """Toggle color mode by clicking the button."""
    page.get_by_label("toggle theme").click(force=True)
    page.get_by_label("toggle theme").click(force=True)


@pytest.fixture(scope="session")
def server():
    """Start ZnDraw server for all screenshot tests."""
    print("\nStarting ZnDraw server...")
    # Use DEVNULL to avoid blocking when pipe buffer fills up
    process = subprocess.Popen(
        [
            "uv",
            "run",
            "zndraw",
            "--no-browser",
            "--file-browser",
            "--file-browser-root",
            str(MISC_DIR.parent),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    if not _wait_for_server(BASE_URL):
        process.kill()
        raise TimeoutError("Server did not start in time")

    print("Server is ready")

    yield BASE_URL

    print("\nShutting down server...")
    subprocess.run(["uv", "run", "zndraw", "--shutdown"], check=True)
    process.wait(timeout=10)


@pytest.fixture
def page():
    """Create a fresh browser and page for each test."""
    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=["--use-angle=gl", "--no-sandbox"],
        )
        # browser = p.chromium.launch(
        #     headless=False,
        #     args=["--no-sandbox"],
        # )
        page = browser.new_page(viewport={"width": 1280, "height": 720})
        yield page
        page.close()
        browser.close()


@pytest.fixture
def capture(page: Page, request):
    """Return functions to capture screenshots in light and dark mode separately."""
    # Ensure output directories exist
    (SCREENSHOTS_DIR / "lightmode").mkdir(parents=True, exist_ok=True)
    (SCREENSHOTS_DIR / "darkmode").mkdir(parents=True, exist_ok=True)

    # Derive screenshot name from test name (test_overview -> overview)
    test_name = request.node.name
    screenshot_name = test_name.removeprefix("test_")

    class Capture:
        def light(self) -> None:
            """Capture light mode screenshot."""
            page.screenshot(
                path=SCREENSHOTS_DIR / "lightmode" / f"{screenshot_name}.png"
            )
            print(f"  Saved: lightmode/{screenshot_name}.png")

        def dark(self) -> None:
            """Capture dark mode screenshot."""
            page.screenshot(
                path=SCREENSHOTS_DIR / "darkmode" / f"{screenshot_name}.png"
            )
            print(f"  Saved: darkmode/{screenshot_name}.png")

        def toggle(self) -> None:
            """Toggle between light and dark mode."""
            _toggle_color_mode(page)

    return Capture()


@pytest.fixture(scope="session")
def bmim_bf4() -> ase.Atoms:
    bf4 = molify.smiles2conformers("[B-](F)(F)(F)F", numConfs=100)
    bmim = molify.smiles2conformers("CCCCN1C=C[N+](=C1)C", numConfs=100)
    mols = [
        molify.pack([[bf4[i]], [bmim[i]]], counts=[1, 1], density=1200)
        for i in range(32)
    ]
    box = molify.pack([mols], counts=[32], density=1200)
    return box


@pytest.fixture
def bmim_bf4_e_f(bmim_bf4) -> list[ase.Atoms]:
    frames = []
    rng = np.random.default_rng(42)
    n_frames = 100

    # Generate correlated energies using sine curve + noise
    t = np.linspace(0, 4 * np.pi, n_frames)
    base_energy = 10 * np.sin(t)
    noise_energy = rng.normal(0, 1.5, n_frames)
    noise_dft = rng.normal(0, 1.5, n_frames)

    for i in range(n_frames):
        atoms = bmim_bf4.copy()
        atoms.calc = SinglePointCalculator(
            atoms,
            energy=base_energy[i] + noise_energy[i],
            forces=rng.uniform(-0.1, 0.1, size=(len(atoms), 3)),
        )
        atoms.calc.results["dft_energy"] = base_energy[i] + noise_dft[i]
        atoms.calc.results["dft_forces"] = rng.uniform(-0.1, 0.1, size=(len(atoms), 3))
        frames.append(atoms)
    return frames
