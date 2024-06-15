import json

import znjson

from zndraw.utils import ASEConverter


def test_ase_converter(s22):
    atoms_json = json.dumps(s22[0], cls=znjson.ZnEncoder.from_converters([ASEConverter]))
    atoms = json.loads(atoms_json, cls=znjson.ZnDecoder.from_converters([ASEConverter]))
    assert atoms == s22[0]
