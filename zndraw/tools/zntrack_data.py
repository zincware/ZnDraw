"""Handling ZnTrack data requests."""
import ase
import zntrack


def get_atoms_via_dvc(name_and_attr, remote, rev) -> list[ase.Atoms]:
    """Load Atoms via ZnTrack.

    Attributes
    ----------
    name_and_attr : str
        A string of the form "node_name:attr" where "node_name.attr : list[Atoms]"
    remote : str
        The remote to use
    rev: str
        The revision to use
    """
    node_name, attr = name_and_attr.split(":")
    node = zntrack.from_rev(node_name, remote=remote, rev=rev)
    return getattr(node, attr)
