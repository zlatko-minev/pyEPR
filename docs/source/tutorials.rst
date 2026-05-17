.. _tutorials:

Tutorial Notebooks
==================

These tutorials are Jupyter notebooks that can be run interactively.
Tutorials 3, 5, and 6 require only ``pip install pyEPR-quantum`` — no Ansys licence needed.
Tutorials 1, 2, and 4 require a live Ansys HFSS session.

No-HFSS tutorials
-----------------

.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Tutorial 3 — Circuit QED Parameters

      E_J, E_C, L_J, I_c conversions; transmon model. **No HFSS required.**

   .. grid-item-card:: Tutorial 5 — Generic Junction Potential & Fluxonium

      Exact cosine diagonalization for fluxonium; custom V(φ); asymmetric SQUIDs. **No HFSS required.**

   .. grid-item-card:: Tutorial 6 — EPR without HFSS

      Supply freqs, Ljs, φ_zpf directly. Full χ matrix without any EM solver. **No HFSS required.**

HFSS tutorials
--------------

.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Tutorial 1 — Startup Example

      End-to-end workflow: HFSS eigenmode simulation → EPR extraction → χ matrix. *Requires Ansys HFSS.*

   .. grid-item-card:: Tutorial 2 — Dielectric Loss EPR

      Dielectric energy participation, loss rates, HFSS fields calculator. *Requires Ansys HFSS.*

   .. grid-item-card:: Tutorial 4 — Parametric Sweeps

      HFSS Optimetrics: linear, log, and file-based parametric sweeps. *Requires Ansys HFSS.*

.. toctree::
   :hidden:
   :maxdepth: 1

   _tutorial_notebooks/Tutorial 3.  toolbox_circuits
   _tutorial_notebooks/Tutorial 5. Generic junction potential and fluxonium EPR
   _tutorial_notebooks/Tutorial 6. EPR without HFSS — purely numerical workflow
   _tutorial_notebooks/Tutorial 1.  Startup example
   _tutorial_notebooks/Tutorial 2.  Field calculations - dielectric energy participation ratios (EPRs)
   _tutorial_notebooks/Tutorial 4. Parametric sweep options
