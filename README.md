[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)
[![PyPI version](https://badge.fury.io/py/zndraw.svg)](https://badge.fury.io/py/zndraw)
[![DOI](https://img.shields.io/badge/arXiv-2402.08708-red)](https://arxiv.org/abs/2402.08708)
[![codecov](https://codecov.io/gh/zincware/ZnDraw/graph/badge.svg?token=3GPCKH1BBX)](https://codecov.io/gh/zincware/ZnDraw)
!['Threejs](https://img.shields.io/badge/threejs-black?style=for-the-badge&logo=three.js&logoColor=white)

# ZnDraw

Welcome to ZnDraw, a powerful tool for visualizing and interacting with your trajectories.

## Installation

It is recommended to install ZnDraw from PyPi via:

```bash
pip install zndraw
```

## Quick Start

Visualize your trajectories with a single command:

```bash
zndraw <file>
```

> [!NOTE]
> ZnDraw's webapp-based approach allows you to use port forwarding to work with trajectories on remote systems.

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/darkmode/overview.png#gh-dark-mode-only "ZnDraw UI")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/lightmode/overview.png#gh-light-mode-only "ZnDraw UI")

## Multi-User and Multi-Client Support

ZnDraw supports multiple users and clients. Connect one or more Python clients to your ZnDraw instance:

1. Click on `Python access` in the ZnDraw UI.
2. Connect using the following code:

```python
from zndraw import ZnDraw

vis = ZnDraw(url="http://localhost:1234", token="<your-token>")
```

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/darkmode/python.png#gh-dark-mode-only "ZnDraw Python Client")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/lightmode/python.png#gh-light-mode-only "ZnDraw Python Client")

The `vis` object provides direct access to your visualized scene. It inherits from `abc.MutableSequence`, so any changes you make are reflected for all connected clients.

```python
from ase.collections import s22
vis.extend(list(s22))
```

## Additional Features

You can modify various aspects of the visualization:

- `vis.camera`
- `vis.points`
- `vis.selection`
- `vis.step`
- `vis.figures`
- `vis.bookmarks`
- `vis.geometries`

For example, to add a geometry:

```python
from zndraw import Box

vis.geometries = [Box(position=[0, 1, 2])]
```

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/darkmode/box.png#gh-dark-mode-only "ZnDraw Geometries")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/lightmode/box.png#gh-light-mode-only "ZnDraw Geometries")

## Analyzing Data

ZnDraw enables you to analyze your data and generate plots using [Plotly](https://plotly.com/). It automatically detects available properties and offers a convenient drop-down menu for selection.

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/darkmode/analysis.png#gh-dark-mode-only "ZnDraw Analysis")
![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/lightmode/analysis.png#gh-light-mode-only "ZnDraw Analysis")

ZnDraw will look for the `step` and `atom` index in the [customdata](https://plotly.com/python/reference/scatter/#scatter-customdata)`[0]` and `[1]` respectively to highlight the steps and atoms.

## Writing Extensions

Make your tools accessible via the ZnDraw UI by writing an extension:

```python
from zndraw import Extension

class AddMolecule(Extension):
    name: str

    def run(self, vis, **kwargs) -> None:
        structures = kwargs["structures"]
        vis.append(structures[self.name])
        vis.step = len(vis) - 1

vis.register(AddMolecule, run_kwargs={"structures": s22}, public=True)
vis.socket.wait()  # This can be ignored when using Jupyter
```

The `AddMolecule` extension will appear for all `tokens` and can be used by any client.

# Hosted Version

A hosted version of ZnDraw is available at https://zndraw.icp.uni-stuttgart.de . To upload data, use:

```bash
zndraw <file> --url https://zndraw.icp.uni-stuttgart.de
```

## Self-Hosting

To host your own version of ZnDraw, use the following `docker-compose.yaml` setup:

```yaml
version: "3.9"

services:
  zndraw:
    image: pythonf/zndraw:latest
    command: --no-standalone /src/file.xyz
    volumes:
      - /path/to/files:/src
    restart: unless-stopped
    ports:
      - 5003:5003
    depends_on:
      - redis
      - worker
    environment:
      - FLASK_STORAGE=redis://redis:6379/0
      - FLASK_AUTH_TOKEN=super-secret-token

  worker:
    image: pythonf/zndraw:latest
    entrypoint: celery -A zndraw_app.make_celery worker --loglevel=info -P eventlet
    volumes:
      - /path/to/files:/src
    restart: unless-stopped
    depends_on:
      - redis
    environment:
      - FLASK_STORAGE=redis://redis:6379/0
      - FLASK_SERVER_URL="http://zndraw:5003"
      - FLASK_AUTH_TOKEN=super-secret-token

  redis:
    image: redis:latest
    restart: always
    environment:
      - REDIS_PORT=6379
```

If you want to host zndraw as subdirectory `domain.com/zndraw` you need to adjust the environmental variables as well as update `base: "/",` in the `app/vite.config.ts` before building the ap..

# References

If you use ZnDraw in your research and find it helpful please cite us.

```bibtex
@misc{elijosiusZeroShotMolecular2024,
  title = {Zero {{Shot Molecular Generation}} via {{Similarity Kernels}}},
  author = {Elijo{\v s}ius, Rokas and Zills, Fabian and Batatia, Ilyes and Norwood, Sam Walton and Kov{\'a}cs, D{\'a}vid P{\'e}ter and Holm, Christian and Cs{\'a}nyi, G{\'a}bor},
  year = {2024},
  eprint = {2402.08708},
  archiveprefix = {arxiv},
}
```

# Acknowledgements

The creation of ZnDraw was supported by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) in the framework of the priority program SPP 2363, “Utilization and Development of Machine Learning for Molecular Applications - Molecular Machine Learning” Project No. 497249646. Further funding though the DFG under Germany's Excellence Strategy - EXC 2075 - 390740016 and the Stuttgart Center for Simulation Science (SimTech) was provided.
