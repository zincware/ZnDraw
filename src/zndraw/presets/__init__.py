"""Bundled visual presets shipped with ZnDraw."""

from pathlib import Path

PRESETS_DIR = Path(__file__).parent


def list_bundled_presets() -> list[Path]:
    """Return paths to all bundled preset JSON files."""
    return sorted(PRESETS_DIR.glob("*.json"))
