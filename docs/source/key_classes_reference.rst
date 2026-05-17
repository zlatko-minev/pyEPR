Key Classes
===========

.. contents:: On this page
   :local:
   :depth: 2

pyEPR's public API is organized around three main classes that map directly
onto the three stages of the EPR analysis pipeline:

.. code-block:: text

   ProjectInfo            — configure: paths, junctions, dissipative elements
        │
        ▼
   DistributedAnalysis    — simulate: connect to HFSS, extract EPR fields, save HDF5
        │
        ▼
   QuantumAnalysis        — quantize: load HDF5, diagonalize Hamiltonian, report results

All three are importable from the top-level ``pyEPR`` namespace:

.. code-block:: python

   import pyEPR as epr

   pinfo = epr.ProjectInfo(...)
   eprd  = epr.DistributedAnalysis(pinfo)
   epra  = epr.QuantumAnalysis(eprd.data_filename)

``DistributedAnalysis`` and ``QuantumAnalysis`` can also be used with their
legacy names ``pyEPR_HFSSAnalysis`` and ``pyEPR_Analysis`` (kept for
backwards compatibility).


.. _project-info:

ProjectInfo
-----------

``ProjectInfo`` is the configuration object.  It stores:

* Paths to the Ansys project, design, and setup.
* Junction descriptions — HFSS variable names, rectangle and polyline object
  names, and optional capacitance variables.
* Dissipative element names (substrate bodies, metal surfaces, seams) for
  loss-tangent analysis.
* A reference to the live ``HfssDesign`` COM wrapper once connected.

Key attributes:

.. code-block:: python

   pinfo.junctions          # OrderedDict of junction parameter dicts
   pinfo.dissipative        # Dict with keys 'dielectrics_bulk', 'dielectric_surfaces', etc.
   pinfo.setup              # HfssSetup wrapper (eigenmode or driven-modal)
   pinfo.design             # HfssDesign wrapper

Full API: :class:`pyEPR.project_info.ProjectInfo`


.. _distributed-analysis:

DistributedAnalysis
--------------------

``DistributedAnalysis`` performs the microwave / field-extraction half of the
EPR pipeline.  It:

1. Connects to the HFSS design via ``ProjectInfo``.
2. Iterates over all solved variations (nominal + optimetric sweep points).
3. For each variation, integrates the electric and magnetic energy density
   over each junction rectangle to compute the inductive participation
   ratio :math:`p_{mj}`.
4. Saves all results (frequencies, participation ratios, convergence data)
   to an HDF5 file.

Key methods:

.. code-block:: python

   eprd.do_EPR_analysis()             # main entry point — runs the full extraction
   eprd.quick_plot_frequencies(swp)   # plot solved eigenfrequencies vs. sweep variable
   eprd.hfss_report_full_convergence()# print/plot adaptive-pass convergence
   eprd.data_filename                 # path to the saved HDF5 results file

Full API: :class:`pyEPR.core_distributed_analysis.DistributedAnalysis`


.. _quantum-analysis:

QuantumAnalysis
----------------

``QuantumAnalysis`` performs the quantum / Hamiltonian half of the pipeline.
It loads saved EPR data and:

1. Constructs the Josephson Hamiltonian from the participation ratios and
   junction parameters.
2. Diagonalizes it numerically using QuTiP (``analyze_all_variations``).
3. Returns and plots the dressed frequencies, anharmonicities
   :math:`\alpha_m`, and dispersive shifts :math:`\chi_{mm'}`.

Key methods:

.. code-block:: python

   epra.analyze_all_variations(cos_trunc=8, fock_trunc=7)
   epra.report_results(swp_variable='Lj1', numeric=True)
   epra.plot_hamiltonian_results(swp_variable='Lj1')
   epra.quick_plot_mode(mode_idx, variation, ...)

**Truncation parameters**

* ``cos_trunc`` — number of terms kept in the cosine expansion of the
  Josephson potential (:math:`\cos\hat\phi \approx \sum_n`).
  Higher values improve accuracy for highly anharmonic modes.
  Typical range: 6–10.
* ``fock_trunc`` — Fock-space dimension per mode.  Total Hilbert space
  size grows as ``fock_trunc ** n_modes``, so keep this small for
  many-mode systems.  Typical range: 6–10.

Full API: :class:`pyEPR.core_quantum_analysis.QuantumAnalysis`


solution_types module
---------------------

``pyEPR.solution_types`` exposes canonical solution-type constants and
helpers for handling the AEDT 2021.2+ renamed solution-type strings:

.. code-block:: python

   from pyEPR.solution_types import normalize, DRIVEN_MODAL_NAMES, is_drivenmodal

   normalize("HFSS Modal Network")           # → "DrivenModal"
   normalize("HFSS Hybrid Terminal Network") # → "DrivenTerminal"
   is_drivenmodal("HFSS Modal Network")      # → True

Full API: :mod:`pyEPR.solution_types`


calcs subpackage
-----------------

``pyEPR.calcs`` contains the low-level analytic and numeric routines:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Module
     - Contents
   * - ``calcs.basic``
     - General EM / qubit parameter formulas
   * - ``calcs.constants``
     - Physical constants in SI units
   * - ``calcs.convert``
     - Unit conversion helpers (frequency, energy, etc.)
   * - ``calcs.hamiltonian``
     - Hamiltonian matrix construction and manipulation
   * - ``calcs.transmon``
     - Analytic transmon formulas (charge-basis diagonalization)
   * - ``calcs.back_box_numeric``
     - Numerical black-box diagonalization (EPR numerical method)

These are used internally by ``QuantumAnalysis`` but are also importable
directly for custom calculations:

.. code-block:: python

   from pyEPR.calcs.transmon import transmon_get_spectrum_charge_basis
   from pyEPR.calcs.convert import Convert
