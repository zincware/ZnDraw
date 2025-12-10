Python interface
================

The ``zndraw`` package provides a Python interface to interact with the visualisation tool.
To use this API, you need to have a running instance of the ZnDraw web server.

.. note::

    You can run a local webserver by using the command line interface.

    .. code:: console

        $ zndraw file.xyz --port 1234


.. code:: python

    from zndraw import ZnDraw

    vis = ZnDraw(url="http://localhost:1234", room="my-room")


.. note::

    In ZnDraw each visualisation is associated with a room name.
    You find this room name in the URL of the visualisation.
    This room name can be used to interact with the visualisation using the Python API.
    Additionally, you can use the room name to share the visualisation with others or view
    the visualisation from different angles in different browser tabs.


Authentication
--------------

ZnDraw supports optional user authentication. You can provide a username and password:

.. code:: python

    vis = ZnDraw(
        url="http://localhost:1234",
        room="my-room",
        user="my-username",
        password="my-password"
    )

If no user is provided, the server will assign a guest username.


Working with Frames
-------------------

The ``vis`` object provides a Python interface to interact with the visualisation.
Most basically, it behaves like a Python list of `ase.Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_ objects.
Modifying the list will update the visualisation in real-time.

.. code:: python

    import ase.io as aio

    frames = aio.read("file.xyz", index=":")
    vis.extend(frames)

You can also read frames from the visualisation:

.. code:: python

    # Get the current frame
    atoms = vis[vis.step]

    # Iterate over all frames
    for atoms in vis:
        print(atoms)

    # Get a slice of frames
    subset = vis[10:20]
