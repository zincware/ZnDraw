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
def main(file: str, port: int = typer.Option(5123)):
    globals.file = file
    app.run(port=port)
