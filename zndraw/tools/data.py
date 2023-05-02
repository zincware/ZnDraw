from zndraw import shared


def serialize_atoms(start: int, stop: int):
    result = {"position": [], "force": [], "box": []}
    try:
        for step in range(start, stop):
            atoms = shared.config.get_atoms(step=int(step))
            result["position"].append(atoms.get_positions().tolist())
            result["box"].append(atoms.get_cell().diagonal().tolist())
            # TODO MAKE THIS OPTIONAL!!, also energy, etc.
            # try:
            #     result["force"].append(atoms.get_forces().tolist())
            # except:
            #     result["force"].append(np.zeros_like(atoms.get_positions()).tolist())
        return result
    except KeyError:
        return result
