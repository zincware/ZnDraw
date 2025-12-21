User Interface Guide
====================

This guide covers the ZnDraw web interface features that complement the Python API.


Rooms Management
----------------

ZnDraw supports multiple visualization rooms. Each room is an independent workspace.

.. image:: /_static/screenshots/lightmode/rooms.png
   :class: only-light
   :alt: Rooms management page

.. image:: /_static/screenshots/darkmode/rooms.png
   :class: only-dark
   :alt: Rooms management page

Access the rooms page at ``/rooms`` to:

- View all active rooms
- See frame counts per room
- Create new rooms
- Navigate between workspaces


File Browser
------------

When started with ``--file-browser``, ZnDraw provides a file browser interface.

.. code:: console

    $ zndraw --file-browser --file-browser-root /path/to/data

.. image:: /_static/screenshots/lightmode/file_browser.png
   :class: only-light
   :alt: File browser

.. image:: /_static/screenshots/darkmode/file_browser.png
   :class: only-dark
   :alt: File browser

The file browser allows you to:

- Navigate directories
- Load trajectory files directly
- Access remote data


Theme Support
-------------

ZnDraw supports both light and dark themes. Toggle between themes using the theme button in the interface.
The theme preference is persisted across sessions.

All UI components including Plotly figures automatically adapt to the selected theme.
