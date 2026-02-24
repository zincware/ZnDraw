"""tqdm subclass that reports progress to a ZnDraw room."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any
from uuid import uuid4

from tqdm.auto import tqdm as std_tqdm

if TYPE_CHECKING:
    from zndraw.client import ZnDraw


class ZnDrawTqdm(std_tqdm):  # type: ignore[type-arg]
    """tqdm subclass that sends progress updates to a ZnDraw room.

    Uses tqdm's built-in ``mininterval`` throttling to limit REST calls.
    Terminal output is suppressed (headless).

    Parameters
    ----------
    *args : Any
        Positional arguments forwarded to tqdm.
    vis : ZnDraw
        The ZnDraw client instance (must have ``vis.api``).
    description : str
        Human-readable label shown in the UI.
    mininterval : float
        Minimum seconds between REST updates (default 0.5).
    **kwargs : Any
        Keyword arguments forwarded to tqdm.

    Examples
    --------
    >>> for atoms in ZnDrawTqdm(frames, vis=vis, description="Loading"):
    ...     vis.append(atoms)
    """

    def __init__(
        self,
        *args: Any,
        vis: ZnDraw,
        description: str = "Processing...",
        mininterval: float = 0.5,
        **kwargs: Any,
    ) -> None:
        kwargs["mininterval"] = mininterval
        super().__init__(*args, **kwargs)
        self._vis_api = vis.api
        self._progress_id = str(uuid4())
        self._vis_api.progress_start(
            self._progress_id, description, unit=kwargs.get("unit", "it")
        )

    def display(self, *args: Any, **kwargs: Any) -> None:
        """Suppress terminal output (headless mode)."""

    def update(self, n: float | None = 1) -> bool | None:
        """Advance the progress bar and send a throttled REST update."""
        displayed = super().update(n)
        if displayed:
            d = self.format_dict
            self._vis_api.progress_update(
                self._progress_id,
                n=d["n"],
                total=d.get("total"),
                elapsed=d["elapsed"],
                unit=d.get("unit", "it"),
            )
        return displayed

    def close(self) -> None:
        """Complete the progress tracker and clean up."""
        if not self.disable:
            self._vis_api.progress_complete(self._progress_id)
        super().close()
