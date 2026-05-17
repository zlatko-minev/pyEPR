.. _tutorials:

Tutorial Notebooks
==================

These tutorials are Jupyter notebooks that can be run interactively.
Tutorials 3, 5, and 6 require only ``pip install pyEPR-quantum``
and run entirely without Ansys HFSS.

.. grid:: 1
   :gutter: 2

   .. grid-item-card:: Tutorial 3 — Circuit QED Parameters

      E_J, E_C, L_J, I_c conversions; transmon model. **No HFSS required.**

   .. grid-item-card:: Tutorial 5 — Generic Junction Potential & Fluxonium

      Exact cosine diagonalization for fluxonium; custom V(φ); asymmetric SQUIDs. **No HFSS required.**

   .. grid-item-card:: Tutorial 6 — EPR without HFSS

      Supply freqs, Ljs, φ_zpf directly. Full χ matrix without any EM solver. **No HFSS required.**

.. toctree::
   :hidden:
   :maxdepth: 1

   _tutorial_notebooks/Tutorial 3.  toolbox_circuits
   _tutorial_notebooks/Tutorial 5. Generic junction potential and fluxonium EPR
   _tutorial_notebooks/Tutorial 6. EPR without HFSS — purely numerical workflow
