Quick-start examples
=====================

.. contents:: On this page
   :local:
   :depth: 2

.. note::

   For interactive, step-by-step walkthroughs see the
   `Jupyter notebook tutorials
   <https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks>`_
   and the `example HFSS files
   <https://github.com/zlatko-minev/pyEPR/tree/master/_example_files>`_.


Example 1 — Full HFSS-to-Hamiltonian pipeline
-----------------------------------------------

This is the canonical pyEPR workflow: connect to a live Ansys HFSS session,
extract EPR data, and compute the quantum Hamiltonian.

The device here is a two-qubit (Alice & Bob) plus one cavity chip, already
simulated in HFSS in eigenmode.  See :ref:`hfss-setup` for how to prepare
the HFSS project.

.. code-block:: python

   import pyEPR as epr

   # ── Step 1: Connect to HFSS ───────────────────────────────────────────────
   pinfo = epr.ProjectInfo(
       project_path = r'C:\sim_folder',
       project_name = r'cavity_with_two_qubits',
       design_name  = r'Alice_Bob',
   )

   # ── Step 2: Describe the Josephson junctions ──────────────────────────────
   # Each junction needs four HFSS object/variable names:
   #   Lj_variable : HFSS variable holding the junction inductance (e.g. "Lj_alice")
   #   rect        : name of the lumped-RLC rectangle in HFSS
   #   line        : polyline spanning the junction (defines current orientation / ZPF sign)
   #   Cj_variable : HFSS variable holding the junction capacitance (optional but recommended)
   pinfo.junctions['jAlice'] = {
       'Lj_variable': 'Lj_alice', 'rect': 'rect_alice',
       'line': 'line_alice',      'Cj_variable': 'Cj_alice',
   }
   pinfo.junctions['jBob'] = {
       'Lj_variable': 'Lj_bob', 'rect': 'rect_bob',
       'line': 'line_bob',      'Cj_variable': 'Cj_bob',
   }
   pinfo.validate_junction_info()  # raises if object/variable names are not found in HFSS

   # ── Step 2b (optional): Dissipative elements ──────────────────────────────
   # Supply HFSS object names to compute participation in dielectric / surface losses
   pinfo.dissipative['dielectrics_bulk']    = ['si_substrate']
   pinfo.dissipative['dielectric_surfaces'] = ['substrate_top_interface']

   # ── Step 3: EPR field extraction ─────────────────────────────────────────
   eprd = epr.DistributedAnalysis(pinfo)
   eprd.quick_plot_frequencies('Lj_alice')   # sanity-check: plot solved eigenfrequencies
   eprd.hfss_report_full_convergence()        # print/plot convergence metrics
   eprd.do_EPR_analysis()                     # compute p_mj for every mode × junction

   # ── Step 4: Quantum Hamiltonian ───────────────────────────────────────────
   epra = epr.QuantumAnalysis(eprd.data_filename)
   epra.analyze_all_variations(cos_trunc=8, fock_trunc=7)

   # ── Step 5: Report results ────────────────────────────────────────────────
   swp = 'Lj_alice'   # optimetric sweep variable
   epra.plot_hamiltonian_results(swp_variable=swp)
   epra.report_results(swp_variable=swp, numeric=True)
   epra.quick_plot_mode(0, 0, 1, numeric=True, swp_variable=swp)


Example 2 — Post-processing saved data (no live HFSS needed)
--------------------------------------------------------------

Once ``do_EPR_analysis()`` has run and saved an HDF5 file, you can close
HFSS entirely and re-run the quantum analysis at any time — on any machine,
any platform.

.. code-block:: python

   import pyEPR as epr

   # Point directly at the saved HDF5 data file
   data_file = r'C:\sim_folder\cavity_with_two_qubits.hdf5'

   epra = epr.QuantumAnalysis(data_file)
   epra.analyze_all_variations(cos_trunc=8, fock_trunc=7)
   epra.report_results(numeric=True)

This is the recommended workflow for:

* Sweeping truncation parameters (``cos_trunc``, ``fock_trunc``) without
  re-running the simulation.
* Running on macOS or Linux where the HFSS COM interface is not available.
* Sharing results with collaborators who do not have HFSS installed.


Example 3 — Using pyEPR without Ansys HFSS
-------------------------------------------

See :ref:`without-hfss` for the full guide.  The short version: pyEPR's
``QuantumAnalysis`` only needs a dictionary of participation ratios and
frequencies — not an HFSS connection.  You can populate that dictionary from
any source (Palace, custom FEM code, analytic estimates).

.. code-block:: python

   import pyEPR as epr
   import pandas as pd
   import numpy as np

   # Suppose your solver gives you these participation ratios:
   #   modes: ['qubit', 'cavity']
   #   junctions: ['junction_1']
   freqs_GHz = pd.Series({'qubit': 5.0, 'cavity': 7.2})   # GHz
   Lj_nH     = pd.Series({'junction_1': 10.0})             # nH

   # Build the p_mj DataFrame (rows = modes, columns = junctions)
   p_mj = pd.DataFrame(
       {'junction_1': [0.95, 0.01]},
       index=['qubit', 'cavity'],
   )

   # See the without_hfss guide for how to pass these into QuantumAnalysis
   # via the HamiltonianResultsContainer interface.

See :ref:`without-hfss` for a complete working example.


Video tutorials
---------------

.. raw:: html

   <div style="overflow-x:auto;">
   <table>
   <tr>
     <td style="padding:8px; text-align:center;">
       <a href="https://www.youtube.com/watch?v=fSRYvD-ITnQ&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=1">
         <img src="https://img.youtube.com/vi/fSRYvD-ITnQ/mqdefault.jpg"
              alt="Tutorial 1 – Overview" width="240"><br>
         <b>Tutorial 1 – Overview</b>
       </a>
     </td>
     <td style="padding:8px; text-align:center;">
       <a href="https://www.youtube.com/watch?v=ZTi1pb6wSbE&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=2">
         <img src="https://img.youtube.com/vi/ZTi1pb6wSbE/mqdefault.jpg"
              alt="Tutorial 2 – Setup" width="240"><br>
         <b>Tutorial 2 – Environment setup</b>
       </a>
     </td>
     <td style="padding:8px; text-align:center;">
       <a href="https://www.youtube.com/watch?v=L79nlXY2w4s&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=3">
         <img src="https://img.youtube.com/vi/L79nlXY2w4s/mqdefault.jpg"
              alt="Tutorial 3 – Packages & Config" width="240"><br>
         <b>Tutorial 3 – Packages &amp; config</b>
       </a>
     </td>
   </tr>
   </table>
   </div>


Jupyter notebooks
-----------------

The most thorough way to learn pyEPR is through the interactive notebooks:

* `View notebooks online (nbviewer)
  <https://nbviewer.jupyter.org/github/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks/>`_
* `Download notebooks + example HFSS files
  <https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks>`_
