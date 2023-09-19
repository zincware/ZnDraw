[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)
[![PyPI version](https://badge.fury.io/py/zndraw.svg)](https://badge.fury.io/py/zndraw)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8304530.svg)](https://doi.org/10.5281/zenodo.8304530)
[![codecov](https://codecov.io/gh/zincware/ZnDraw/graph/badge.svg?token=3GPCKH1BBX)](https://codecov.io/gh/zincware/ZnDraw)
!['Threejs](https://img.shields.io/badge/threejs-black?style=for-the-badge&logo=three.js&logoColor=white)

# ZnDraw

Install via `pip install zndraw`. If you have `pip install pywebview` installed,
ZnDraw will open in a dedicated window.

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
