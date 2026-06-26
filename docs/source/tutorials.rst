.. _tutorials:

Tutorial Notebooks
==================

Seven Jupyter notebook tutorials covering the full pyEPR workflow.
Tutorials 1, 2, 4, and 7 require a live Ansys HFSS session.
Tutorials 3, 5, and 6 run entirely with ``pip install pyEPR-quantum`` — no Ansys licence needed.

.. grid:: 1
   :gutter: 3

   .. grid-item-card:: Tutorial 1 — End-to-End HFSS Workflow
      :link: _tutorial_notebooks/Tutorial 1.  Startup example.html
      :link-type: url

      :bdg-warning:`Ansys HFSS required`

      Connect to HFSS, extract EPR participation ratios, and diagonalize the full Josephson Hamiltonian to get qubit frequencies, anharmonicities, and the χ matrix.

   .. grid-item-card:: Tutorial 2 — Dielectric Loss & EPR Fields
      :link: _tutorial_notebooks/Tutorial 2.  Field calculations - dielectric energy participation ratios (EPRs).html
      :link-type: url

      :bdg-warning:`Ansys HFSS required`

      Compute dielectric energy participation ratios, loss rates, and use the HFSS fields calculator for surface and volume integrals.

   .. grid-item-card:: Tutorial 3 — Circuit QED Parameters
      :link: _tutorial_notebooks/Tutorial 3.  toolbox_circuits.html
      :link-type: url

      :bdg-success:`No HFSS required`

      Convert between E_J, E_C, L_J, and I_c; explore the transmon charge-basis model and energy spectrum.

   .. grid-item-card:: Tutorial 4 — Parametric Sweeps
      :link: _tutorial_notebooks/Tutorial 4. Parametric sweep options.html
      :link-type: url

      :bdg-warning:`Ansys HFSS required`

      Set up and run HFSS Optimetrics sweeps (linear, log, file-based), save fields, and batch-process results across sweep points.

   .. grid-item-card:: Tutorial 5 — Fluxonium & Generic Junction Potentials
      :link: _tutorial_notebooks/Tutorial 5. Generic junction potential and fluxonium EPR.html
      :link-type: url

      :bdg-success:`No HFSS required`

      Diagonalize the exact cosine potential for fluxonium; define custom V(φ); handle asymmetric SQUIDs with large zero-point fluctuations.

   .. grid-item-card:: Tutorial 6 — Numerical EPR without HFSS
      :link: _tutorial_notebooks/Tutorial 6. EPR without HFSS — purely numerical workflow.html
      :link-type: url

      :bdg-success:`No HFSS required`

      Supply frequencies, junction inductances, and φ_zpf directly to get the full χ matrix — no EM solver needed. `Run on Binder ↗ <https://mybinder.org/v2/gh/zlatko-minev/pyEPR/master?filepath=_tutorial_notebooks%2FTutorial%206.%20EPR%20without%20HFSS%20%E2%80%94%20purely%20numerical%20workflow.ipynb>`__

   .. grid-item-card:: Tutorial 7 — EPR through PyAEDT (gRPC), no COM
      :link: _tutorial_notebooks/Tutorial 7. EPR through PyAEDT (gRPC) — no COM.html
      :link-type: url

      :bdg-warning:`Ansys HFSS required`

      Run the same EPR extraction through Ansys's official PyAEDT API entirely over gRPC — no COM. Attaches to a running AEDT session, feeds pyEPR's own diagonalizer, and matches the COM path digit-for-digit.

.. toctree::
   :hidden:
   :maxdepth: 1

   _tutorial_notebooks/Tutorial 1.  Startup example
   _tutorial_notebooks/Tutorial 2.  Field calculations - dielectric energy participation ratios (EPRs)
   _tutorial_notebooks/Tutorial 3.  toolbox_circuits
   _tutorial_notebooks/Tutorial 4. Parametric sweep options
   _tutorial_notebooks/Tutorial 5. Generic junction potential and fluxonium EPR
   _tutorial_notebooks/Tutorial 6. EPR without HFSS — purely numerical workflow
   _tutorial_notebooks/Tutorial 7. EPR through PyAEDT (gRPC) — no COM
