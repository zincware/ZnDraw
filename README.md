[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)
[![PyPI version](https://badge.fury.io/py/zndraw.svg)](https://badge.fury.io/py/zndraw)
[![DOI](https://img.shields.io/badge/arXiv-2402.08708-red)](https://arxiv.org/abs/2402.08708)
[![codecov](https://codecov.io/gh/zincware/ZnDraw/graph/badge.svg?token=3GPCKH1BBX)](https://codecov.io/gh/zincware/ZnDraw)
!['Threejs](https://img.shields.io/badge/threejs-black?style=for-the-badge&logo=three.js&logoColor=white)

# ZnDraw

Install via `pip install zndraw`. If you have `pip install pywebview` installed,
ZnDraw will open in a dedicated window instead of your default browser.

## CLI

You can use ZnDraw to view a file using the CLI `zndraw traj.xyz`. Supported
file formats include everything that `ase.io` can read and additionally `h5`
files in the H5MD standard.

If you want to view the frames while they are added to the scene you can use
`zndraw -mp traj.xyz`. See `zndraw --help` for more CLI options.

## Python

ZnDraw provides a Python interface. The `zndraw.ZnDraw` object offers `append`,
`extend` as well as assignment operations. More information is available in the
example notebook.

```python
from zndraw import ZnDraw
import ase

vis = ZnDraw()

vis.socket.sleep(2) # give it some time to fully connect
vis[0] = ase.Atoms(
  "H2O", positions=[[0.75, -0.75, 0], [0.75, 0.75, 0], [0, 0, 0]]
  )
```

ZnDraw also provides an interface to the Python
[logging](https://docs.python.org/3/library/logging.html) library, including
support for formatters and different logging levels.

```python
import logging

log = logging.getLogger(__name__)
log.addHandler(vis.get_logging_handler())
log.critical("Critical Message")
```

## Modifier

You can register `modifier` to change the scene via the interactions menu.

```python
import typing as t

from zndraw import ZnDraw, Extension
import ase

vis = ZnDraw()

class MyModifier(UpdateScene):
  name: str = "H2O"

  def run(self, vis: ZnDraw, **kwargs) -> None:
    vis.append(molecule(self.name))

vis.register_modifier(
  MyModifier, public=True, run_kwargs={}
)
```

## User Interface

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_ui.png "ZnDraw UI")

![ZnDraw UI3](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_draw.png "ZnDraw UI3")

# Development

ZnDraw is developed using https://python-poetry.org/. Furthermore, the
javascript packages have to be installed using https://www.npmjs.com/.

```bash
cd zndraw/static/
npm install
```

# Docker

This package is available as a docker container and can be used via

```bash
docker run -pythonf/zndraw
```

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
