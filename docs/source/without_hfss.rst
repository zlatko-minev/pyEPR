.. _without-hfss:

Using pyEPR without Ansys HFSS
===============================

.. contents:: On this page
   :local:
   :depth: 2

pyEPR has two independent functional layers.  The **quantum analysis layer**
(``DistributedAnalysis``, ``QuantumAnalysis``) requires only eigenfrequencies
and energy-participation ratios as inputs — it does **not** need a live Ansys
HFSS connection or even an Ansys installation.

This page explains the three main ways to use pyEPR's quantum analysis
without HFSS:

1. :ref:`load-saved-hdf5` — re-run quantum analysis on a previously saved dataset.
2. :ref:`bring-your-own-data` — supply participation ratios from any solver.
3. :ref:`future-solvers` — upcoming and third-party solver integrations.


.. _load-saved-hdf5:

Re-loading saved HDF5 data
---------------------------

After ``DistributedAnalysis.do_EPR_analysis()`` completes, pyEPR saves all
field-integral results to an HDF5 file.  This file is fully self-contained:
you can close HFSS, move to another machine, and run the quantum analysis at
any time.

.. code-block:: python

   import pyEPR as epr

   # Path to the HDF5 file written by do_EPR_analysis()
   data_file = r'C:\sim_results\my_device.hdf5'   # Windows path
   # data_file = '/home/user/sim_results/my_device.hdf5'  # Linux/macOS

   epra = epr.QuantumAnalysis(data_file)
   epra.analyze_all_variations(cos_trunc=8, fock_trunc=7)
   epra.report_results(numeric=True)
   epra.plot_hamiltonian_results()

This is the recommended workflow for:

* Sweeping ``cos_trunc`` / ``fock_trunc`` without re-running simulations.
* Sharing results with collaborators who do not have HFSS.
* Running the quantum analysis on macOS or Linux.
* Archiving simulation results for reproducibility.


.. _bring-your-own-data:

Bring your own participation ratios
-------------------------------------

If you have eigenfrequencies and participation ratios from another source
(a different FEM solver, an analytic calculation, or a circuit model), you
can feed them directly into ``QuantumAnalysis`` by constructing the
``HamiltonianResultsContainer`` that it expects.

The key quantities pyEPR needs are:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Quantity
     - Description
   * - ``freqs_hfss_GHz``
     - Eigenfrequencies in GHz (one per mode)
   * - ``Qs``
     - Quality factors (use ``np.inf`` if not computing dissipation)
   * - ``Ljs``
     - Junction inductances in Henries (one per junction)
   * - ``Cjs``
     - Junction capacitances in Farads (one per junction; use ``0`` if unknown)
   * - ``p_mj``
     - DataFrame of inductive participation ratios: rows = modes, columns = junctions
   * - ``sign_mj``
     - DataFrame of participation signs ±1 (same shape as ``p_mj``)

Below is a minimal worked example for a single transmon qubit coupled to a
readout cavity:

.. code-block:: python

   import numpy as np
   import pandas as pd
   import pyEPR as epr
   from pyEPR.core_quantum_analysis import QuantumAnalysis, HamiltonianResultsContainer

   # ── Device parameters ────────────────────────────────────────────────────
   modes     = ['qubit', 'cavity']
   junctions = ['jj']

   freqs_GHz = np.array([5.0,  7.2])   # GHz
   Qs        = np.array([1e6,  1e4])   # quality factors

   Lj_H = 10e-9   # 10 nH junction inductance
   Cj_F = 2e-15   # 2 fF junction capacitance

   # Participation ratios (rows=modes, cols=junctions)
   # Qubit is nearly all junction; cavity has tiny participation
   p_mj = pd.DataFrame(
       [[0.95],
        [0.02]],
       index=modes, columns=junctions,
   )
   sign_mj = pd.DataFrame(
       [[+1],
        [+1]],
       index=modes, columns=junctions,
   )

   # ── Build the results container ───────────────────────────────────────────
   # HamiltonianResultsContainer is an OrderedDict keyed by "variation" label.
   # A single variation (no parametric sweep) uses key '0'.
   results = HamiltonianResultsContainer()
   results['0'] = {
       'freqs_hfss_GHz': pd.Series(freqs_GHz, index=modes),
       'Qs':             pd.Series(Qs,         index=modes),
       'Ljs':            pd.Series([Lj_H],     index=junctions),
       'Cjs':            pd.Series([Cj_F],     index=junctions),
       'p_mj':           p_mj,
       'sign_mj':        sign_mj,
   }

   # ── Save and reload via QuantumAnalysis ───────────────────────────────────
   import tempfile, os
   tmp = tempfile.mktemp(suffix='.hdf5')
   results.save(tmp)

   epra = QuantumAnalysis(tmp)
   epra.analyze_all_variations(cos_trunc=8, fock_trunc=7)
   epra.report_results(numeric=True)

   os.unlink(tmp)   # clean up temp file

.. note::

   The sign convention for ``sign_mj`` follows the current orientation of the
   polyline defined in HFSS.  For custom data, a consistent choice of ``+1``
   is fine unless you need to track relative phase between junctions.


.. _future-solvers:

Third-party and future solver integrations
------------------------------------------

pyEPR's quantum analysis layer is deliberately solver-agnostic.  The only
hard requirement is a set of eigenfrequencies and field-energy integrals.
Below is the current status of integrations beyond Ansys HFSS.

**PyAEDT** (`github.com/ansys/pyaedt <https://github.com/ansys/pyaedt>`_)
   Ansys's official cross-platform Python scripting library for AEDT.
   PyAEDT can drive HFSS simulations on Windows, and in principle the field
   results it extracts can be passed into pyEPR's quantum analysis layer.
   Direct integration is not yet packaged; contributions welcome.

**Palace** (`github.com/awslabs/palace <https://github.com/awslabs/palace>`_)
   AWS open-source, GPU-accelerated FEM solver with an eigenmode solver.
   Palace outputs field data in a format that can be post-processed to yield
   participation ratios.  A pyEPR–Palace bridge is a natural future
   contribution; the :ref:`bring-your-own-data` approach above is the
   recommended path in the meantime.

**Qiskit Metal**
   Metal uses pyEPR internally.  If you are using Metal's EPR analysis
   renderer, the data flows through pyEPR automatically.

**Custom / analytic**
   For simple geometries (e.g. lumped transmon estimates), you can compute
   participation ratios analytically and plug them directly into
   ``HamiltonianResultsContainer`` as shown in :ref:`bring-your-own-data`.

.. tip::

   If you have added a solver integration or have working Palace / OpenEMS /
   COMSOL scripts that produce pyEPR-compatible output, please open a pull
   request or issue at
   `github.com/zlatko-minev/pyEPR <https://github.com/zlatko-minev/pyEPR/issues>`_.
   The core maintainers are happy to help integrate and document new solver
   backends.
