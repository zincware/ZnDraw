from zndraw.tasks import get_generator_from_filename, FileIO
import ase


def test_get_generator_from_filename():
    file = FileIO(name="https://raw.githubusercontent.com/LarsSchaaf/Guaranteed-Non-Local-Molecular-Dataset/main/gnl-dataset/GNL-v0.2/gnl-v0.2-test.xyz")
    generator = get_generator_from_filename(file)
    
    assert isinstance(next(iter(generator)), ase.Atoms)

