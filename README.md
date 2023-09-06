[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)
[![PyPI version](https://badge.fury.io/py/zndraw.svg)](https://badge.fury.io/py/zndraw)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8304530.svg)](https://doi.org/10.5281/zenodo.8304530)
!['Threejs](https://img.shields.io/badge/threejs-black?style=for-the-badge&logo=three.js&logoColor=white)

# ZnDraw

Install via `pip install zndraw`. If you have `pip install pywebview` installed,
ZnDraw will open in a dedicated window.

## CLI

You can use ZnDraw to view a file `zndraw traj.xyz`.
Supported file formats include everything that `ase.io` can read and additionally `h5` files in the H5MD standard.

If you want to view the frames while they are added to the scene you can use `zndraw -mp traj.xyz`.
See `zndraw --help` for more CLI options.

## Python
ZnDraw provides a Python interface.
The `zndraw.ZnDraw` object offers `append`, `extend` as well as assignment operations. More information is available in the example notebook.

```python
from zndraw import ZnDraw
import ase

zndraw = ZnDraw()

zndraw.socket.sleep(2) # give it some time to fully connect
zndraw[0] = ase.Atoms(
  "H2O", positions=[[0.75, -0.75, 0], [0.75, 0.75, 0], [0, 0, 0]]
  )
```

## User Interface

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_ui.png "ZnDraw UI")

![ZnDraw UI2](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_protein.gif "ZnDraw Protein GIF")

![ZnDraw UI3](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_draw.png "ZnDraw UI3")

# Development

ZnDraw is developed using https://python-poetry.org/. Furthermore, the
javascript packages have to be installed using https://www.npmjs.com/.

```bash
cd zndraw/static/
npm install
```
