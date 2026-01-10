<div align="center">

![Logo](https://raw.githubusercontent.com/zincware/ZnDraw/refs/heads/main/docs/source/_static/zndraw-light.svg#gh-light-mode-only)
![Logo](https://raw.githubusercontent.com/zincware/ZnDraw/refs/heads/main/docs/source/_static/zndraw-dark.svg#gh-dark-mode-only)

[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)
[![PyPI version](https://badge.fury.io/py/zndraw.svg)](https://badge.fury.io/py/zndraw)
[![DOI](https://img.shields.io/badge/DOI-10.1038/s41467--025--60963--3-blue)](https://doi.org/10.1038/s41467-025-60963-3)
[![codecov](https://codecov.io/gh/zincware/ZnDraw/graph/badge.svg?token=3GPCKH1BBX)](https://codecov.io/gh/zincware/ZnDraw)
[![Discord](https://img.shields.io/discord/1034511611802689557)](https://discord.gg/7ncfwhsnm4)
[![Documentation Status](https://readthedocs.org/projects/zndraw/badge/?version=latest)](https://zndraw.readthedocs.io/en/latest/?badge=latest)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/zincware/ZnDraw)

</div>

# ZnDraw

A Python-first visualization and editing tool for atomic structures with real-time collaboration.

## Installation

You can install `zndraw` into your Python environment via pip:
```bash
pip install zndraw
```
or set it up as uv tool to run anywhere:
```bash
uvx zndraw <file>
```

## Quick Start

Visualize your trajectories with a single command:

```bash
zndraw <file>
```

> [!NOTE]
> ZnDraw's webapp-based approach allows you to use port forwarding to work with trajectories on remote systems.

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/darkmode/overview.png#gh-dark-mode-only "ZnDraw UI")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/lightmode/overview.png#gh-light-mode-only "ZnDraw UI")

## Python API

ZnDraw supports multiple clients connecting to the same visualization. Each visualization session is identified by a **room name** visible in the URL.

```python
from zndraw import ZnDraw

vis = ZnDraw(url="http://localhost:1234", room="my-room")
```

### Authentication

For protected deployments, provide credentials:

```python
vis = ZnDraw(
    url="http://localhost:1234",
    room="my-room",
    user="username",
    password="password"
)
```

If no credentials are provided, the server assigns a guest user.

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/darkmode/python_connection.png#gh-dark-mode-only "ZnDraw Python Client")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/lightmode/python_connection.png#gh-light-mode-only "ZnDraw Python Client")

## Working with Frames

The `vis` object behaves like a Python list of [ase.Atoms](https://wiki.fysik.dtu.dk/ase/ase/atoms.html) objects. Changes are synchronized in real-time across all connected clients.

```python
from ase.collections import s22

# Add structures
vis.extend(list(s22))

# Access current frame
atoms = vis[vis.step]

# Iterate over frames
for atoms in vis:
    print(atoms)

# Slice operations
subset = vis[10:20]
```

## Scene Properties

Control various aspects of the visualization:

```python
vis.selection   # Currently selected atoms
vis.step        # Current frame index
vis.figures     # Plotly figures
vis.bookmarks   # Saved frame annotations
vis.geometries  # 3D geometry overlays (dict-like)
vis.sessions    # Session configuration
```

## Geometries

Add 3D geometry overlays to your visualization:

```python
from zndraw.geometries import Box, Sphere, Arrow, Camera, Curve

vis.geometries["box"] = Box(position=(0, 1, 2))
```

Available geometry types: `Sphere`, `Arrow`, `Bond`, `Curve`, `Cell`, `Floor`, `Box`, `Plane`, `Shape`, `Camera`.

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/darkmode/geometries.png#gh-dark-mode-only "ZnDraw Geometries")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/lightmode/geometries.png#gh-light-mode-only "ZnDraw Geometries")

## Analysis

ZnDraw integrates with [Plotly](https://plotly.com/) for interactive data visualization. It automatically detects available properties and provides selection menus.

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/darkmode/analysis_1d.png#gh-dark-mode-only "ZnDraw Analysis")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/lightmode/analysis_1d.png#gh-light-mode-only "ZnDraw Analysis")

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/darkmode/analysis_2d.png#gh-dark-mode-only "ZnDraw Analysis 2D")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/docs/source/_static/screenshots/lightmode/analysis_2d.png#gh-light-mode-only "ZnDraw Analysis 2D")

## Extensions

Create custom tools accessible via the ZnDraw UI:

```python
from zndraw import Extension
from molify import smiles2atoms

class AddMolecule(Extension):
    smiles: str

    def run(self, vis, **kwargs) -> None:
        vis.append(smiles2atoms(self.smiles))
        vis.step = len(vis) - 1

vis.register(AddMolecule, public=True)
vis.wait()
```

Extensions can be registered as `public` (available to all clients) or private (only to the registering client).

## Hosted Version

A hosted version is available at https://zndraw.icp.uni-stuttgart.de

```bash
zndraw <file> --url https://zndraw.icp.uni-stuttgart.de
```

## Self-Hosting

ZnDraw can be deployed using Docker:

See the [docker/](docker/) directory for complete deployment configurations:
- [Standalone](docker/standalone/README.md) - Simple single-instance deployment for personal use or small teams
- [Production](docker/production/README.md) - Horizontal scaling with nginx load balancer for high load

## References

If you use ZnDraw in your research, please cite:

```bibtex
@article{elijosius2025zero,
  title = {Zero-shot molecular generation via similarity kernels},
  author = {Elijo{\v s}ius, Rokas and Zills, Fabian and Batatia, Ilyes and Norwood, Sam Walton and Kov{\'a}cs, D{\'a}vid P{\'e}ter and Holm, Christian and Cs{\'a}nyi, G{\'a}bor},
  journal = {Nature Communications},
  volume = {16},
  pages = {5479},
  year = {2025},
  doi = {10.1038/s41467-025-60963-3},
  url = {https://doi.org/10.1038/s41467-025-60963-3},
}
```

## Acknowledgements

The creation of ZnDraw was supported by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) in the framework of the priority program SPP 2363, "Utilization and Development of Machine Learning for Molecular Applications - Molecular Machine Learning" Project No. 497249646. Further funding through the DFG under Germany's Excellence Strategy - EXC 2075 - 390740016 and the Stuttgart Center for Simulation Science (SimTech) was provided.
