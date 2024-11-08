import typing as t

import ase

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw

from . import UpdateScene


class NewScene(UpdateScene):
    """Clear the scene, deleting all atoms and points."""

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        del vis[vis.step + 1 :]
        vis.points = []
        vis.append(ase.Atoms())
        vis.selection = []
        step = len(vis) - 1
        vis.step = step
        vis.bookmarks = vis.bookmarks | {step: "New Scene"}
        vis.camera = {"position": [0, 0, 20], "target": [0, 0, 0]}


class ClearTools(UpdateScene):
    """Clear the tools, removing all guiding points and undoing any selection."""

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        vis.points = []
        vis.selection = []
