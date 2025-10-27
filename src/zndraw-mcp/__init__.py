from fastmcp import FastMCP

from zndraw import ZnDraw

mcp = FastMCP("ZnDraw Connection")


@mcp.tool
def get_atoms(url, room, user, step: int) -> str:
    """Retrieve the atoms representation from the ZnDraw server.

    Parameter
    ----------
    url : str
        The URL of the ZnDraw server.
    room : str
        The room name to connect to.
    user : str
        The user name for the connection.
    step : int
        The frame number to retrieve.

    """
    vis = ZnDraw(url=url, room=room, user=user)
    atoms = vis[step]
    return str(atoms)


@mcp.tool
def get_step(url, room, user) -> int:
    """Retrieve the current step/frame number from the ZnDraw server.

    Parameter
    ----------
    url : str
        The URL of the ZnDraw server.
    room : str
        The room name to connect to.
    user : str
        The user name for the connection.

    """
    vis = ZnDraw(url=url, room=room, user=user)
    return vis.step


@mcp.tool
def set_step(url, room, user, step: int) -> None:
    """Set the current step/frame number on the ZnDraw server.

    Parameter
    ----------
    url : str
        The URL of the ZnDraw server.
    room : str
        The room name to connect to.
    user : str
        The user name for the connection.
    step : int
        The frame number to set.

    """
    vis = ZnDraw(url=url, room=room, user=user)
    vis.step = step


@mcp.tool
def get_length(url, room, user) -> int:
    """Retrieve the total number of steps/frames from the ZnDraw server.

    Parameter
    ----------
    url : str
        The URL of the ZnDraw server.
    room : str
        The room name to connect to.
    user : str
        The user name for the connection.

    """
    vis = ZnDraw(url=url, room=room, user=user)
    return len(vis)


if __name__ == "__main__":
    mcp.run()
