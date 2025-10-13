"""Screenshot management for ZnDraw rooms."""

import dataclasses
from pathlib import Path
from typing import Literal


@dataclasses.dataclass
class Screenshot:
    """Screenshot metadata.

    Attributes
    ----------
    id : int
        Unique identifier for the screenshot (incrementing).
    format : str
        Image format (png, jpeg, or webp).
    size : int
        File size in bytes.
    width : int | None
        Image width in pixels, if known.
    height : int | None
        Image height in pixels, if known.
    """

    id: int
    format: Literal["png", "jpeg", "webp"]
    size: int
    width: int | None = None
    height: int | None = None


class ScreenshotManager:
    """Manages screenshot storage for a room.

    Uses incrementing integer IDs for simplicity and performance.

    Parameters
    ----------
    room_id : str
        The room identifier.
    storage_path : str
        Base storage path for screenshots.
    """

    VALID_FORMATS = {"png", "jpeg", "webp"}

    def __init__(self, room_id: str, storage_path: str):
        self.room_id = room_id
        self.screenshot_dir = Path(storage_path) / room_id / "screenshots"
        self.screenshot_dir.mkdir(parents=True, exist_ok=True)

    def _get_next_id(self) -> int:
        """Get next available screenshot ID.

        Returns
        -------
        int
            Next sequential ID.
        """
        existing = [
            int(f.stem) for f in self.screenshot_dir.glob("*.*") if f.stem.isdigit()
        ]
        return max(existing, default=0) + 1

    def _get_filepath(self, screenshot_id: int, format: str) -> Path:
        """Get file path for screenshot ID and format.

        Parameters
        ----------
        screenshot_id : int
            Screenshot identifier.
        format : str
            Image format extension.

        Returns
        -------
        Path
            Full file path.
        """
        return self.screenshot_dir / f"{screenshot_id}.{format}"

    def save(
        self,
        image_data: bytes,
        format: str,
        width: int | None = None,
        height: int | None = None,
    ) -> Screenshot:
        """Save screenshot and return metadata.

        Parameters
        ----------
        image_data : bytes
            The image file data.
        format : str
            Image format (png, jpeg, or webp).
        width : int | None
            Image width in pixels.
        height : int | None
            Image height in pixels.

        Returns
        -------
        Screenshot
            Metadata for the saved screenshot.

        Raises
        ------
        ValueError
            If format is not supported.
        """
        if format not in self.VALID_FORMATS:
            raise ValueError(
                f"Invalid format '{format}'. Must be one of {self.VALID_FORMATS}"
            )

        screenshot_id = self._get_next_id()
        filepath = self._get_filepath(screenshot_id, format)
        filepath.write_bytes(image_data)

        return Screenshot(
            id=screenshot_id,
            format=format,
            size=len(image_data),
            width=width,
            height=height,
        )

    def get(self, screenshot_id: int) -> tuple[Path, Screenshot] | None:
        """Get screenshot file path and metadata.

        Parameters
        ----------
        screenshot_id : int
            Unique screenshot identifier.

        Returns
        -------
        tuple[Path, Screenshot] | None
            File path and metadata, or None if not found.
        """
        for filepath in self.screenshot_dir.glob(f"{screenshot_id}.*"):
            return filepath, Screenshot(
                id=screenshot_id,
                format=filepath.suffix[1:],
                size=filepath.stat().st_size,
                width=None,
                height=None,
            )
        return None

    def list(self, limit: int = 20, offset: int = 0) -> list[Screenshot]:
        """List screenshots sorted by ID (newest first).

        Parameters
        ----------
        limit : int
            Maximum number of screenshots to return.
        offset : int
            Number of screenshots to skip.

        Returns
        -------
        list[Screenshot]
            List of screenshot metadata.
        """
        # Get all files with numeric stems
        files = [f for f in self.screenshot_dir.glob("*.*") if f.stem.isdigit()]
        # Sort by ID descending
        files.sort(key=lambda f: int(f.stem), reverse=True)

        # Apply pagination
        result = []
        for filepath in files[offset : offset + limit]:
            screenshot_id = int(filepath.stem)
            result.append(
                Screenshot(
                    id=screenshot_id,
                    format=filepath.suffix[1:],
                    size=filepath.stat().st_size,
                    width=None,
                    height=None,
                )
            )
        return result

    def count(self) -> int:
        """Count total screenshots.

        Returns
        -------
        int
            Total number of screenshots.
        """
        return len([f for f in self.screenshot_dir.glob("*.*") if f.stem.isdigit()])

    def delete(self, screenshot_id: int) -> bool:
        """Delete screenshot by ID.

        Parameters
        ----------
        screenshot_id : int
            Unique screenshot identifier.

        Returns
        -------
        bool
            True if deleted, False if not found.
        """
        for filepath in self.screenshot_dir.glob(f"{screenshot_id}.*"):
            filepath.unlink()
            return True
        return False
