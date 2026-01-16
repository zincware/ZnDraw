Keyboard Shortcuts
==================

ZnDraw provides keyboard shortcuts for the navigation and interaction.
Shortcuts are organized by interaction mode.

.. note::

   Keyboard shortcuts are disabled when typing in input fields (text boxes, search, etc.).


Global Shortcuts
----------------

These shortcuts work in any mode:

+------------------------+-----------------------------------------------+
| Shortcut               | Action                                        |
+========================+===============================================+
| ``←`` / ``→``          | Navigate to previous / next frame             |
+------------------------+-----------------------------------------------+
| ``Shift + ←`` /        | Jump to previous / next bookmark              |
| ``Shift + →``          |                                               |
+------------------------+-----------------------------------------------+
| ``Space``              | Toggle playback                               |
+------------------------+-----------------------------------------------+
| ``B``                  | Add bookmark to current frame                 |
+------------------------+-----------------------------------------------+
| ``C``                  | Center camera on selection (or all particles) |
+------------------------+-----------------------------------------------+
| ``R``                  | Rotate camera 15° around view axis            |
+------------------------+-----------------------------------------------+
| ``Shift + R``          | Rotate camera 45° around view axis            |
+------------------------+-----------------------------------------------+
| ``Ctrl + R``           | Rotate camera in opposite direction           |
+------------------------+-----------------------------------------------+
| ``I``                  | Toggle info boxes (property inspector)        |
+------------------------+-----------------------------------------------+


View Mode Shortcuts
-------------------

These shortcuts are available in the default view mode:

+------------------------+-----------------------------------------------+
| Shortcut               | Action                                        |
+========================+===============================================+
| ``Ctrl/Cmd + A``       | Select all particles                          |
+------------------------+-----------------------------------------------+
| ``X``                  | Enter drawing mode                            |
+------------------------+-----------------------------------------------+
| ``E``                  | Enter editing mode                            |
+------------------------+-----------------------------------------------+


Drawing Mode
------------

Drawing mode allows interactive curve creation. Press ``X`` to toggle.

.. image:: /_static/screenshots/lightmode/drawing_mode.png
   :class: only-light
   :alt: Drawing mode

.. image:: /_static/screenshots/darkmode/drawing_mode.png
   :class: only-dark
   :alt: Drawing mode

+------------------------+-----------------------------------------------+
| Shortcut               | Action                                        |
+========================+===============================================+
| ``X``                  | Exit drawing mode (return to view)            |
+------------------------+-----------------------------------------------+
| ``Click``              | Add point to active curve                     |
+------------------------+-----------------------------------------------+

**How Drawing Mode Works:**

1. Press ``X`` to enter drawing mode
2. A drawing marker appears at your cursor position
3. Click to add control points to the active curve
4. The marker turns red when over an invalid position
5. Press ``X`` again to exit and return to view mode

.. note::

   Drawing mode only works with Curve geometries that have static positions.
   The room is temporarily locked while drawing to prevent conflicts.


Editing Mode
------------

Editing mode allows interactive transformation of geometries. Press ``E`` to toggle.

.. image:: /_static/screenshots/lightmode/editing_mode.png
   :class: only-light
   :alt: Editing mode

.. image:: /_static/screenshots/darkmode/editing_mode.png
   :class: only-dark
   :alt: Editing mode

+------------------------+-----------------------------------------------+
| Shortcut               | Action                                        |
+========================+===============================================+
| ``E``                  | Exit editing mode (return to view)            |
+------------------------+-----------------------------------------------+
| ``T``                  | Cycle transform mode                          |
|                        | (translate → rotate → scale)                  |
+------------------------+-----------------------------------------------+
| ``X`` (hold)           | Constrain to X axis                           |
+------------------------+-----------------------------------------------+
| ``Y`` (hold)           | Constrain to Y axis                           |
+------------------------+-----------------------------------------------+
| ``Z`` (hold)           | Constrain to Z axis                           |
+------------------------+-----------------------------------------------+
| ``↑`` / ``↓``          | Adjust transform sensitivity                  |
+------------------------+-----------------------------------------------+
| ``Shift + ↑/↓``        | Larger sensitivity adjustment                 |
+------------------------+-----------------------------------------------+
| ``S``                  | Save frame edits                              |
+------------------------+-----------------------------------------------+
| ``Delete`` /           | Delete selected curve markers                 |
| ``Backspace``          |                                               |
+------------------------+-----------------------------------------------+

**How Editing Mode Works:**

1. Press ``E`` to enter editing mode
2. The editing indicator appears showing the current transform mode
3. Select geometry instances by clicking on them
4. Use transform gizmo to translate, rotate, or scale
5. Press ``T`` to cycle between transform modes
6. Hold ``X``, ``Y``, or ``Z`` to constrain movement to a single axis
7. Press ``S`` to save your changes
8. Press ``E`` again to exit and return to view mode

**Axis Constraints:**

.. image:: /_static/screenshots/lightmode/editing_axis_constraint.png
   :class: only-light
   :alt: Axis constraint indicator

.. image:: /_static/screenshots/darkmode/editing_axis_constraint.png
   :class: only-dark
   :alt: Axis constraint indicator

When holding an axis key, a colored chip appears indicating the active constraint:

- **Red chip**: X axis constraint
- **Green chip**: Y axis constraint
- **Blue chip**: Z axis constraint

**Virtual Markers (Curves):**

.. image:: /_static/screenshots/lightmode/curve_editing.png
   :class: only-light
   :alt: Curve editing with virtual markers

.. image:: /_static/screenshots/darkmode/curve_editing.png
   :class: only-dark
   :alt: Curve editing with virtual markers

In editing mode, curves display virtual markers between control points.
Click a virtual marker to insert a new control point at that position.
