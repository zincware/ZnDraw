import ase.build
import pytest


@pytest.fixture
def water():
    return ase.build.molecule("H2O")
