We have

```python
vis.geometries = {
    "particles": {
        "type": "Sphere",
        "entity": "atoms",  # ← Links to atom indices
        "position": "arrays.positions",
        "color": "arrays.colors"
    },
    "forces": {
        "type": "Arrow",
        "entity": "atoms",  # ← SAME entity type, shares atom indices
        "position": "arrays.positions",
        "direction": "calc.forces",
        "color": "#FF0000"
    },
    "bonds": {
        "type": "Cylinder",
        "entity": "bonds",  # ← Different entity type, separate indices
        "indices": "arrays.bonds",  # Bond connectivity pairs
        "color": "#888888"
    },
    "boxes": {
        "type": "Box",
        "entity": "boxes",  # ← Different entity type
        "corners": [[0,0,0], [10,10,10]],
        "color": "#0000FF",
        "selectable": False  # Decorations usually not selectable
    }
}
```

and then `vis.selection = [1, 2, 3]` operates on `entity: "atoms"` and assumes that all geometries sharing that entity type use the same atom indices.
There is `vis.selection_groups = {"my-selection": {"atoms" : [0, 1, 2]}}` to define named selection groups and we can do `vis.selection = "my-selection"` to select that group.

