[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)
[![PyPI version](https://badge.fury.io/py/zndraw.svg)](https://badge.fury.io/py/zndraw)
!['Threejs](https://img.shields.io/badge/threejs-black?style=for-the-badge&logo=three.js&logoColor=white)


# ZnDraw

Install via ``pip install zndraw`` or ``pip install zndraw[webview]`` to open zndraw in a dedicated window.

## CLI

You can use ZnDraw with the CLI ``zndraw atoms.xyz``.
For a full list of arguments use `zndraw --help`.

ZnDraw is designed to work with your Python scripts.
To interface you can inherit from `zndraw.examples.UpdateScene` or follow this base class:

```python
import abc
from pydantic import BaseModel

class UpdateScene(BaseModel, abc.ABC):
    @abc.abstractmethod
    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        pass
```

The ``run`` method expects as inputs
- atom_ids: list[int], the ids of the currently selected atoms
- atoms: ase.Atoms, the configuration as `ase.Atoms` file where atom_ids where selected.
- kwargs: dict could be additional information from the scene

and as an output:
- list[ase.Atoms], a list of ase Atoms objects to display.


You can define the parameters using `pydantic.Field` which will be displayed in the UI.

```python
class MyUpdateCls(UpdateScene):
    steps: int = Field(100, le=1000, ge=1)
    x: float = Field(0.5, le=5, ge=0)
    symbol: str = Field("same")
```

To add your method click on the `+` on the right side of the window.
Your should be able to add your method from the working directory via `module.MyUpdateCls` as long
as it can be imported via `from module import MyUpdateCls`.

![ZnDraw UI](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_ui.png "ZnDraw UI")

![ZnDraw UI2](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_protein.png "ZnDraw UI2")

![ZnDraw UI3](https://raw.githubusercontent.com/zincware/ZnDraw/main/misc/zndraw_draw.png "ZnDraw UI3")
