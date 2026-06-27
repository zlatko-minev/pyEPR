.. _pyaedt-backend:

PyAEDT (gRPC) HFSS backend
==========================

.. contents:: On this page
   :local:
   :depth: 2

pyEPR has always driven Ansys HFSS through a hand-rolled **COM** layer
(:mod:`pyEPR.ansys`), which is Windows-only and tied to a private scripting
interface.  The :mod:`pyEPR.ansys_pyaedt` module performs the **same**
Energy-Participation-Ratio field extraction through **PyAEDT**
(``pyaedt``) — Ansys's official, maintained Python API — **entirely
over gRPC, with no COM**.

The physics is unchanged: the extracted participations feed pyEPR's own
:meth:`~pyEPR.calcs.basic.CalcsBasic.epr_to_zpf` and
:func:`~pyEPR.calcs.back_box_numeric.epr_numerical_diagonalization` exactly as
the COM path does, and the results match the COM path digit-for-digit
(:math:`p_{mj} = 0.9755` on the demo transmon).

For COM users nothing changes — :mod:`pyEPR.ansys` and
:class:`~pyEPR.DistributedAnalysis` are untouched, and PyAEDT remains an
optional dependency.


Why a gRPC backend
------------------

* **No COM / pywin32.**  The transport is gRPC, the default for AEDT since
  2022 R2, instead of COM (.NET, Windows-only).
* **Attaches to a running session.**  PyAEDT reads the project's
  ``.aedt.lock`` and attaches to the AEDT instance that owns it, avoiding the
  stale-session and *project-locked* errors common with COM's
  Running-Object-Table lookup.
* **Maintenance moves to Ansys.**  ``pyaedt`` is officially versioned
  and supported, rather than a private COM interface pyEPR must track by hand
  across AEDT releases.


Installation
------------

PyAEDT is declared as an optional extra, so ``import pyEPR`` never requires it:

.. code-block:: bash

   pip install "pyEPR-quantum[pyaedt]"


Usage
-----

Describe the project and its junctions with a :class:`~pyEPR.ProjectInfo`, just
as in the COM workflow, but pass ``do_connect=False`` — PyAEDT does the
connecting.  Each junction uses the same keys pyEPR already understands:
``Lj_variable`` (or ``Lj_henries``) and ``line``.

.. code-block:: python

   import pyEPR as epr
   from pyEPR.ansys_pyaedt import PyaedtDistributedAnalysis

   pinfo = epr.ProjectInfo(
       project_path=r"C:\HFSS_Projects",
       project_name="EPR_Demo_Project",
       design_name="EPR_Sample_Demo",
       setup_name="EPR_Scan",
       do_connect=False,            # PyAEDT does the connecting, not COM
   )
   pinfo.junctions["j1"] = {
       "Lj_variable": "Lj_Transmon",   # HFSS variable, e.g. '12.9201nH'
       "line": "Junction_line",        # polyline across the junction
   }

   eprd = PyaedtDistributedAnalysis(pinfo, aedt_version="2026.1")
   eprd.do_EPR_analysis()              # pure gRPC field extraction

   f_ND, chi_ND = eprd.analyze()       # pyEPR's own diagonalizer

After :meth:`~pyEPR.ansys_pyaedt.PyaedtDistributedAnalysis.do_EPR_analysis` the
participations, signs, eigenfrequencies, and junction inductances are available
as plain NumPy arrays on the object (``PJ``, ``SJ``, ``freqs_GHz``, ``Ljs``) —
the same shapes pyEPR's diagonalizer expects.

.. note::

   This backend reads ``Cj_farads`` directly (a float, in farads) rather than
   ``Cj_variable`` (an HFSS variable name used by the COM backend). Omit both to
   treat the junction capacitance as zero.


How it works over gRPC
----------------------

The field-calculator token stack (``EnterQty`` / ``ClcMaterial`` / ``EnterVol``
/ ``EnterLine`` / ``Integrate``) and the eigenmode normalization
(``Solutions.EditSources``) are identical to the COM path — they all run fine
over gRPC.  Only two things differ:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * -
     - COM (:mod:`pyEPR.ansys`)
     - PyAEDT (:mod:`pyEPR.ansys_pyaedt`)
   * - field-calculator handle
     - ``design.GetModule('FieldsReporter')``
     - ``hfss.ofieldsreporter``
   * - read a result back
     - ``ClcEval`` + ``GetTopEntryValue``
     - ``CalculatorWrite`` to a ``.fld`` file

The single operation that does **not** survive gRPC is the stateful
``ClcEval`` / ``GetTopEntryValue`` read-back round-trip; writing the calculator
result to a file with ``CalculatorWrite`` and reading it back is the equivalent
that works over gRPC.


See also
--------

* `Tutorial 7 — EPR through PyAEDT (gRPC) <_tutorial_notebooks/Tutorial 7. EPR through PyAEDT (gRPC) — no COM.html>`__ — a full worked example against a solved transmon.
* :mod:`pyEPR.ansys` — the original COM backend, unchanged.
