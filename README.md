[![zincware](https://img.shields.io/badge/Powered%20by-zincware-darkcyan)](https://github.com/zincware)

# ZnDraw

Install via ``pip install zndraw``.

## CLI

You can use ZnDraw with the CLI ``zndraw atoms.xyz``.
For a full list of arguments use `zndraw --help`.

To interface with ``zndraw --update-function module.function`` you need to be able to import via ``from module import function``.

The ZnDraw function expects as inputs
- atom_ids: list[int], the ids of the currently selected atoms
- atoms: ase.Atoms, the configuration as `ase.Atoms` file where atom_ids where selected.
and as an output:
- list[ase.Atoms], a list of ase Atoms objects to display.

```python
def function(atom_ids: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
    ...
```