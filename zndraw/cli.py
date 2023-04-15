import typer
from zndraw import app, globals

cli = typer.Typer()


# @app.command()
# def main(
#     file: str,
#     select_canvas: bool = typer.Option(True),
#     select_atoms: bool = typer.Option(False),
#     bonds: bool = typer.Option(False),
# ) -> None:
#     """ZnDraw: Visualize Molecules

#     With 'zndraw molecule.xyz' you can visualize the molecule in your browser.

#     Attributes
#     ----------
#     file : str
#         Path to the xyz file.
#     select_canvas : bool
#         Add a canvas that helps to select points in space.
#     select_atoms : bool
#         Allow selection of atoms.
#     """
#     graph = io.read_file(file)

#     app = viewer.DashApp(
#         graph=graph, select_canvas=select_canvas, select_atoms=select_atoms
#     )
#     app.update_layout()
#     app.add_canvas_slider()
#     app.plot_atoms()
#     if bonds:
#         app.plot_bonds()
#     app.run_dash_server_click()


@cli.command()
def main(
    file: str = typer.Argument(..., help="Trajectory File"),
    port: int = typer.Option(5123, help="Port to run the server on"),
    animate: bool = typer.Option(False, help="Animate the trajectory"),
    sphere_size: float = typer.Option(1.0, help="size of the hydrogen sphere"),
    bond_size: float = typer.Option(1.0, help="size of a bond"),
    max_fps: int = typer.Option(1000, help="Maximum frames per second"),
    update_function: str = typer.Option(
        None, help="Path to a python file with an update function 'module.function'."
    ),
    repeat: int = typer.Option(1, help="Repeat the trajectory n times"),
):
    globals.config.file = file
    globals.config.animate = animate
    globals.config.sphere_size = sphere_size
    globals.config.bond_size = bond_size
    globals.config.max_fps = max_fps
    globals.config.update_function = update_function
    globals.config.repeat = (repeat, repeat, repeat)

    app.run(port=port)
