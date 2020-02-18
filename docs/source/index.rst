.. pyEPR documentation master file, created by
   sphinx-quickstart on Wed Jan 15 05:35:04 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*********************************************
Welcome to pyEPR 🍻!
*********************************************

Powerfull analysis and design of quantum microwave devices
***********************************************************
**Version**: |version| | **Date**: |today|


**pyEPR** is an open source, BSD-licensed library providing high-efficiency,
easy-to-use analysis functions and automation for the design of quantum chips based on superconducting quantum  circuits, both distributed and lumped.
pyEPR interfaces the classical distributed microwave analysis with that of quantum structures and Hamiltonians.
It is chiefly based on the `energy participation ratio <https://arxiv.org/abs/1902.10355>`_ approach; however, it has since v0.4 extended to cover a broad range of
design approaches. pyEPR stradels the analysis from Maxwell's to Schrodinger's equations, and converts the solutions of distributed microwve (typically eignmode simulations)
to a fully diagonalized spectrum of the energy levels, couplings, and key parameters of a many-body quantum Hamiltonian.

pyEPR contains both analytic and numeric solutions.


Contents
==================

.. toctree::
   :maxdepth: 3
   :caption: Contents:
   :numbered:

   about.rst
   installation.rst

.. toctree::
   :caption: API Reference:
   :glob:

   api/*



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

