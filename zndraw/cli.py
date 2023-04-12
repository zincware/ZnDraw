import typer
from zndraw import io, viewer

app = typer.Typer()


@app.command()
def main(
    file: str,
    select_canvas: bool = typer.Option(True),
    select_atoms: bool = typer.Option(False),
    bonds: bool = typer.Option(False),
) -> None:
    """Dask4DVC CLI callback.
    Run the DVC graph or DVC experiments in parallel using dask.
    """
    graph = io.read_file(file)

    app = viewer.DashApp(
        graph=graph, select_canvas=select_canvas, select_atoms=select_atoms
    )
    app.update_layout()
    app.add_canvas_slider()
    app.plot_atoms()
    if bonds:
        app.plot_bonds()
    app.run_dash_server_click()
