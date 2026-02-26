"""Filesystem extension for loading files from registered providers."""

import typing as t

from zndraw_joblib.client import Category, Extension


class LoadFile(Extension):
    """Load a file from a registered filesystem provider into the room."""

    category: t.ClassVar[Category] = Category.MODIFIER

    provider_name: str
    path: str
    start: int | None = None
    stop: int | None = None
    step: int | None = None

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        """Read file via the filesystem handler and extend the room."""
        import ase.io
        from ase.io.formats import filetype

        providers: dict[str, t.Any] = kwargs.get("providers") or {}
        fs = providers.get(self.provider_name)
        if fs is None:
            available = list(providers)
            raise ValueError(
                f"Provider '{self.provider_name}' not found. "
                f"Available providers: {available}"
            )
        # Infer format from path extension — fsspec file objects lack .name,
        # so ASE's auto-detection falls back to magic-byte sniffing which
        # fails on text-mode streams.
        fmt = filetype(self.path, read=False)
        with fs.open(self.path, "r") as f:
            atoms_list = ase.io.read(
                f, index=slice(self.start, self.stop, self.step), format=fmt
            )

        if not isinstance(atoms_list, list):
            atoms_list = [atoms_list]

        vis.extend(atoms_list)
