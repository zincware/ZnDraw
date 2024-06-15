import json

import pytest
import znjson
from ase.calculators.singlepoint import SinglePointCalculator

from zndraw.utils import ASEConverter


def test_ase_converter(s22):
    s22[0].connectivity = [[0, 1, 1], [1, 2, 1], [2, 3, 1]]
    s22[3].calc = SinglePointCalculator(s22[3])
    s22[3].calc.results = {"energy": 0.0, "predicted_energy": 1.0}
    s22[4].info = {"key": "value"}

    structures_json = json.dumps(
        s22, cls=znjson.ZnEncoder.from_converters([ASEConverter])
    )
    structures = json.loads(
        structures_json, cls=znjson.ZnDecoder.from_converters([ASEConverter])
    )
    for s1, s2 in zip(s22, structures):
        assert s1 == s2

    assert structures[0].connectivity == [[0, 1, 1], [1, 2, 1], [2, 3, 1]]
    with pytest.raises(AttributeError):
        _ = structures[1].connectivity

    assert structures[3].calc.results == {"energy": 0.0, "predicted_energy": 1.0}

    assert "colors" in structures[0].arrays
    assert "radii" in structures[0].arrays

    assert structures[4].info == {"key": "value"}
