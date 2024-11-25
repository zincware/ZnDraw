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

    vis = ZnDraw(url="http://localhost:1234", token="c91bb84f")


.. note::

    In ZnDraw each visualisation is associated with a unique token.
    You find this token in the URL of the visualisation.
    This token can be used to interact with the visualisation using the Python API.
    Additionally, you can use the token to share the visualisation with others or view
    the visualisation from different angles in different browser tabs.


The ``vis`` object provides a Python interface to interact with the visualisation.
Most basically, it behaves like a Python list of `ase.Atoms <https://wiki.fysik.dtu.dk/ase/ase/atoms.html>`_ objects.
Modifying the list will update the visualisation in real-time.

.. code:: python

    import ase.io as aio

    frames = aio.read("file.xyz", index=":")
    vis.extend(frames)

