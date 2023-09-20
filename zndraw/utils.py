import socket
import znjson
import ase


def get_port(default: int = 1234) -> int:
    """Get an open port."""
    try:
        sock = socket.socket()
        sock.bind(("", default))
        port = 1234
    except OSError:
        sock = socket.socket()
        sock.bind(("", 0))
        port = sock.getsockname()[1]
    finally:
        sock.close()
    return port


class AtomsConverter(znjson.ConverterBase):
    """Encode and decode 'ase.Atoms'."""

    representation = "ase.Atoms"
    instance = ase.Atoms

    def encode(self, obj: ase.Atoms) -> dict:
        """Encode 'ase.Atoms'."""
        atoms_dict = obj.todict()
        from zndraw.bonds import ASEComputeBonds
        ase_bond_calculator = ASEComputeBonds()

        try:
            calc_data = {}
            for key in obj.calc.results:
                value = obj.calc.results[key]
                # if isinstance(value, np.ndarray):
                #     value = value.tolist()
                calc_data[key] = value

            atoms_dict["calc"] = calc_data
        except (RuntimeError, AttributeError):
            pass
        try:
            atoms_dict["connectivity"] = ase_bond_calculator.get_bonds(obj)
        except AttributeError:
            atoms_dict["connectivity"] = []
        return atoms_dict
    
    
    def decode(self, obj: dict) -> ase.Atoms:
        """Decode 'ase.Atoms'."""
        return ase.Atoms(
            symbols=obj["symbols"],
            positions=obj["positions"],
            cell=obj["cell"],
            pbc=obj["pbc"],
        )
