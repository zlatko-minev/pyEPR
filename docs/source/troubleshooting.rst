.. _troubleshooting:

Troubleshooting
===============

.. contents:: On this page
   :local:
   :depth: 2

This page covers the most common errors encountered when installing and
running pyEPR.  For bugs not listed here, please search or open an issue at
`github.com/zlatko-minev/pyEPR/issues
<https://github.com/zlatko-minev/pyEPR/issues>`_.


Installation issues
-------------------

ImportError: cannot import 'attrdict'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**::

    ImportError: Please install python package `AttrDict`.

**Cause:** Old pyEPR versions used ``attrdict``, which is unmaintained and
incompatible with Python 3.10+.

**Fix:** Upgrade to pyEPR ≥ 0.8.5 (uses ``addict`` instead):

.. code-block:: bash

    pip install --upgrade pyEPR-quantum

If you are working from source, pull the latest master and reinstall:

.. code-block:: bash

    git pull origin master
    pip install -e .


pint error: ``system='mks' unknown``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**::

    ValueError: system='mks' unknown

**Fix:** Update pint to a version newer than 0.7.2:

.. code-block:: bash

    pip install --upgrade pint


COM / HFSS connection issues
-----------------------------

COM error on opening HFSS
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**::

    pywintypes.com_error: (-2147352567, 'Exception occurred.', ...)

**Common causes and fixes:**

* The project path contains special characters (apostrophes, colons, spaces
  with unusual encoding). Move the project to a simple path like
  ``C:\sims\my_project``.
* HFSS has popped an error dialog (e.g. "File locked" or "License error").
  Open HFSS manually, dismiss any dialogs, then re-run the script.
* A previous pyEPR session did not release HFSS cleanly. Kill HFSS in
  Task Manager and restart it.

.. note::

   On macOS and Linux, the COM interface is not available natively.
   See :ref:`without-hfss` for how to use pyEPR's quantum analysis layer
   without a live HFSS connection.

HFSS refuses to close after the script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Call ``pinfo.disconnect()`` at the end of your script, or use pyEPR as a
context manager.  If HFSS is already frozen, kill it from Task Manager
(Windows) or Activity Monitor (macOS).

AEDT 2024.1 creates a Hybrid design unexpectedly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom:** The design solution type is ``"HFSS Hybrid Modal Network"``
instead of ``"DrivenModal"``.

**Fix:** Use ``project.new_dm_design(name)`` instead of ``new_design()``.
pyEPR ≥ 0.9.2 calls ``SetSolutionType`` automatically after creation to
enforce the non-hybrid type. See :ref:`hfss-setup` for details.


EPR analysis issues
--------------------

``validate_junction_info()`` raises "object not found"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The HFSS object or variable name does not match exactly (including
capitalisation). Double-check:

* The rectangle name (e.g. ``rect_jj1``) exists in the HFSS 3D Modeler.
* The inductance variable (e.g. ``Lj1``) exists in the design — not a global
  variable prefixed with ``$``.
* The polyline name (e.g. ``line_jj1``) exists and spans the junction rectangle.

Use ``pinfo.validate_junction_info()`` early; it raises a descriptive error
if anything is missing.

Participation ratios are near zero for all modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* The junction rectangle does not have a **Lumped RLC** boundary condition,
  or the inductance variable evaluates to zero or infinity.
* The rectangle is too small relative to the mesh, so very little energy is
  captured in the integration.
* Check in HFSS that the boundary condition is assigned correctly and that
  the inductance variable has a sensible value (e.g. 10 nH).

Negative or greater-than-one participation ratios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This usually indicates a meshing or convergence issue.  Tighten the adaptive
delta-F criterion (try 0.01 %) and re-solve.  Also check that the energy
normalisation is correct — the total energy should equal the sum across all
junctions plus the non-junction energy.

Parametric sweep gives results for only the nominal variation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

"Save Fields and Mesh" was not enabled for the parametric sweep.  To fix:

1. Right-click your Parametric setup in HFSS → **Properties → Options**.
2. Enable **Save Fields and Mesh** for each variation.
3. Re-run the sweep before calling ``do_EPR_analysis()``.

Setup name not found / wrong setup used
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom** (pyEPR ≥ 0.9.4)::

    ValueError: Setup 'MySetup' not found in design. Available setups: [...]

Pass the correct ``setup_name`` to ``ProjectInfo``, or omit it to use the
first setup automatically.  In older versions the specified name was silently
ignored — upgrade to ≥ 0.9.4.


Quantum analysis issues
------------------------

``AttributeError: 'complex' object has no attribute 'norm'``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom**::

    AttributeError: 'complex' object has no attribute 'norm'

**Cause:** QuTiP 5 changed the return type of ``ket.dag() * ket`` from a
1×1 ``Qobj`` to a plain complex scalar.

**Fix:** Upgrade to pyEPR ≥ 0.9.4, which includes the QuTiP 5 compatibility
fix throughout ``back_box_numeric.py``.  pyEPR now requires QuTiP ≥ 5.0.

``KeyError: '0' is not in variations``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This occurs when the HDF5 data file has no variation labelled ``'0'``, which
can happen with some older data files or when mode selection was applied
inconsistently.  Reload the file and inspect ``epra.variations`` to see the
available variation keys.

``[1] not in index`` when selecting a subset of modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When passing a subset of modes (e.g. ``modes=[0, 4]``) to
``analyze_variation``, the internal indices must match the mode labels in the
data. If you simulated 5 modes and selected only modes 0 and 4, the resulting
DataFrame has indices 0 and 4 — not 0 and 1. Pass the same mode list
consistently to ``analyze_all_variations`` and ``report_results``.

Convergence plot raises ``ValueError: Data has no positive values``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Fix:** Upgrade to pyEPR ≥ 0.9.4.  The convergence plotting functions now
fall back to linear scale when delta-F data contains no positive values
(e.g. when the simulation converges on the first pass, reporting Δf = 0).

``FutureWarning`` from pandas in ``analyze_all_variations``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These warnings about multi-dimensional indexing are harmless in current
pandas versions and do not affect results.  They will be addressed in a
future release.


Plotting and reporting issues
------------------------------

Junction sign always shows ``(+)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Fix:** Upgrade to pyEPR ≥ 0.9.4.  The participation print statement
previously used ``if _Smj`` instead of ``if _Smj > 0``; since both ``1``
and ``-1`` are truthy, the sign was always shown as ``(+)``.  The stored
sign value was always correct — only the display was wrong.

``AttributeError: module 'pandas.compat' has no attribute 'StringIO'``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a pandas version incompatibility from pandas < 0.25.  Update pandas:

.. code-block:: bash

    pip install --upgrade pandas

``ValueError: cannot set WRITEABLE flag to True of this array``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This occurs when reading an HDF5 file with numpy 1.16.  Upgrade numpy:

.. code-block:: bash

    pip install --upgrade numpy h5py

``AttributeError: module 'numpy' has no attribute '__config__'``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Update numpy and qutip together:

.. code-block:: bash

    pip install --upgrade numpy qutip


Getting help
------------

If your issue is not listed here:

1. Search `existing issues <https://github.com/zlatko-minev/pyEPR/issues>`_
   — it may already be reported and answered.
2. Open a new issue with the full traceback, your pyEPR version
   (``import pyEPR; print(pyEPR.__version__)``), Python version, operating
   system, and Ansys AEDT version.
