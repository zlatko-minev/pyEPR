About pyEPR
===========

.. |star this repo| image:: https://img.shields.io/github/stars/zlatko-minev/pyEPR?style=social
   :target: https://github.com/zlatko-minev/pyEPR/stargazers
.. |fork this repo| image:: https://img.shields.io/github/forks/zlatko-minev/pyEPR?style=social
   :target: https://github.com/zlatko-minev/pyEPR/fork

:bdg-success:`Open Source — BSD-3` :bdg-info:`Python 3.10+` |star this repo| |fork this repo|

.. contents:: On this page
   :local:
   :depth: 2


Overview
--------

**pyEPR** is an open-source Python library for the automated design and
quantization of superconducting quantum circuits.  It bridges classical
distributed microwave simulation and the quantum world, converting HFSS
(or equivalent) field solutions into a fully diagonalized quantum
Hamiltonian — frequencies, anharmonicities, dispersive shifts, and
cross-Kerr couplings — in a single automated pipeline.

pyEPR has two distinct capability layers:

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Layer
     - Description
   * - **EPR / quantum analysis**
     - Platform-independent post-processing of eigenmode field data.
       Works with data from any electromagnetic solver (HFSS, Palace,
       custom) or from a saved HDF5 file.
   * - **Ansys HFSS automation**
     - Python COM/DCOM wrapper for Ansys HFSS on Windows: geometry,
       boundary conditions, setup, sweep, and field extraction — all
       scriptable from Python.  Originally written before PyAEDT existed;
       the two tools are complementary.


The energy-participation ratio (EPR) method
-------------------------------------------

The EPR method, introduced in
`Minev et al., npj Quantum Information (2021) <https://arxiv.org/abs/2010.00620>`_,
provides a unified, systematic, and efficient approach for computing the
quantum Hamiltonian parameters of superconducting circuits.

The central idea is to characterize each non-linear element (Josephson
junction) by a single dimensionless number — the *energy participation ratio*
:math:`p_{mj}` — which quantifies what fraction of the electromagnetic energy
of mode :math:`m` is stored in junction :math:`j`.

**Participation ratio**

For eigenmode :math:`m` and junction :math:`j`:

.. math::

   p_{mj} \;=\; \frac{\text{inductive energy in junction } j}
                     {\text{total electromagnetic energy of mode } m}
         \;=\; \frac{\frac{1}{2} L_j^{-1} \langle \hat\phi_j^2 \rangle}
                    {E_m}

where :math:`L_j` is the Josephson inductance and :math:`E_m = \hbar\omega_m / 2`
is the zero-point energy of mode :math:`m` (with :math:`\hbar\omega_m` the
mode frequency).  The zero-point phase fluctuation across junction :math:`j`
in mode :math:`m` is:

.. math::

   \varphi_{mj}^\text{zpf} = \sqrt{p_{mj}\,\frac{\hbar\omega_m}{2E_J}}

where :math:`E_J = \hbar^2 / (4e^2 L_j)` is the Josephson energy.

**Hamiltonian parameters from EPR**

Once the participation ratios are extracted from the field simulation,
the leading-order quantum Hamiltonian parameters follow analytically
(to first order in :math:`p_{mj}` with a cosine junction potential
expanded to fourth order):

.. math::

   \hat H / \hbar \;=\;
   \sum_m \omega_m \hat a_m^\dagger \hat a_m
   - \frac{1}{2}\sum_{m} \alpha_m\, \hat a_m^\dagger \hat a_m^\dagger \hat a_m \hat a_m
   - \sum_{m < m'} \chi_{mm'}\, \hat a_m^\dagger \hat a_m \hat a_{m'}^\dagger \hat a_{m'}
   + \ldots

The self-Kerr (anharmonicity) of mode :math:`m` and the cross-Kerr coupling
between modes :math:`m` and :math:`m'` are:

.. math::

   \alpha_m = \sum_j p_{mj}^2\, \frac{e^2}{2C_j} \cdot \frac{1}{\hbar}
            = \sum_j p_{mj}^2\, E_{C,j}/\hbar

.. math::

   \chi_{mm'} = 2\sum_j p_{mj}\, p_{m'j}\, E_{C,j}/\hbar

where :math:`E_{C,j} = e^2 / (2C_j)` is the charging energy of junction
:math:`j`.  For the full derivation valid beyond the perturbative limit,
see `arXiv:2010.00620 <https://arxiv.org/abs/2010.00620>`_.

**Why EPR?**

* **Single simulation** — all Hamiltonian parameters are extracted from one
  eigenmode solve; no need for separate simulations per coupling.
* **Simulator-agnostic** — only eigenmode frequencies and field integrals are
  needed; the method works with any electromagnetic solver.
* **Systematic and scalable** — straightforwardly extends to many modes and
  junctions; participation ratios can be read off from standard energy plots.
* **Validated** — ten-percent to percent-level agreement with experiment
  over five orders of magnitude and across dozens of devices (3D cavities,
  transmons, fluxonium, flip-chip).


What pyEPR is **not**
---------------------

* It is **not** a geometry or layout tool — use
  `Qiskit Metal <https://github.com/Qiskit/qiskit-metal>`_ or
  `KLayout <https://www.klayout.de/>`_ for that.
* It is **not** a general AEDT scripting library — use
  `PyAEDT <https://github.com/ansys/pyaedt>`_ for mesh control, geometry
  parametrization, or driving other AEDT tools.
* It is **not** a circuit solver — it post-processes distributed-field
  solutions, not lumped-element netlists.


Relationship to other tools
----------------------------

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Tool
     - Relationship to pyEPR
   * - `PyAEDT <https://github.com/ansys/pyaedt>`_
     - Ansys's official cross-platform AEDT Python library.  Use for
       geometry, mesh, and solve scripting; pyEPR handles EPR quantization.
       The two are complementary.
   * - `Qiskit Metal <https://github.com/Qiskit/qiskit-metal>`_
     - IBM's chip-layout design tool.  Uses pyEPR internally for EPR
       analysis and Hamiltonian extraction.
   * - `Palace <https://github.com/awslabs/palace>`_
     - AWS open-source FEM solver (parallel, GPU-accelerated).  Does not
       yet have a direct pyEPR integration, but pyEPR's EPR analysis layer
       is solver-agnostic; see :ref:`without-hfss` for how to use pyEPR
       with non-HFSS field data.
   * - `QuTiP <https://qutip.org>`_
     - Used internally by pyEPR for numerical Hamiltonian diagonalization
       (``QuantumAnalysis``).  QuTiP ≥ 5.0 is required.


Citation
--------

If you use pyEPR in your research, please cite:

* **Method paper:** Z. K. Minev, Z. Leghtas, S. O. Mundhada, L. Christakis,
  I. M. Pop, and M. H. Devoret, "Energy-participation quantization of
  Josephson circuits," *npj Quantum Information* **7**, 131 (2021).
  `arXiv:2010.00620 <https://arxiv.org/abs/2010.00620>`_

* **Software:** Z. K. Minev et al., *pyEPR: The energy-participation-ratio
  (EPR) open-source framework for quantum device design* (2021).
  `DOI:10.5281/zenodo.4552482 <https://doi.org/10.5281/zenodo.4552482>`_

BibTeX entries are in ``pyEPR.bib`` at the root of the repository.


References
----------

* Z. K. Minev, Z. Leghtas *et al.*, npj Quantum Information **7**, 131 (2021)
  (`arXiv:2010.00620 <https://arxiv.org/abs/2010.00620>`_)
* Z. K. Minev, Ph.D. Dissertation, Yale University (2018), Chapter 4
  (`arXiv:1902.10355 <https://arxiv.org/abs/1902.10355>`_)
* Z. K. Minev, Z. Leghtas *et al.*, original pyEPR framework (2018)
  (`Zenodo <https://doi.org/10.5281/zenodo.4552482>`_)
