.. pyEPR documentation master file

*********************************************
pyEPR — Energy-Participation-Ratio Framework
*********************************************

**Version**: |version| | **License**: BSD-3-Clause | `GitHub <https://github.com/zlatko-minev/pyEPR>`_ | `PyPI <https://pypi.org/project/pyEPR-quantum/>`_

.. image:: _static/read_me_0.png
   :width: 100%
   :alt: pyEPR overview
   :align: center

----

**pyEPR** is an open-source Python library for the automated design and
quantization of Josephson quantum circuits. It bridges classical distributed
microwave simulation (Ansys HFSS) and quantum circuit Hamiltonians using the
`energy-participation ratio (EPR) <https://arxiv.org/abs/2010.00620>`_ method.

.. grid:: 3
   :gutter: 3
   :margin: 4 4 0 0

   .. grid-item-card:: :octicon:`rocket;1.5em` Simulate
      :text-align: center

      Run an eigenmode simulation in Ansys HFSS with junction inductances.
      No manual circuit diagram needed — just the 3-D geometry.

   .. grid-item-card:: :octicon:`graph;1.5em` Extract
      :text-align: center

      Compute energy participation ratios *p*\ :sub:`mj` and zero-point
      phase fluctuations φ\ :sub:`zpf` for every mode and junction.

   .. grid-item-card:: :octicon:`beaker;1.5em` Diagonalize
      :text-align: center

      Numerically diagonalize the full Josephson Hamiltonian to get
      dressed frequencies, anharmonicities, and the χ matrix.

.. image:: _static/xmon-example.gif
   :width: 55%
   :alt: HFSS junction setup animation — Xmon qubit
   :align: center

----

.. card:: :octicon:`rocket;1em` New in v0.9.6 — cross-platform gRPC backend (Linux · macOS · Windows)
   :class-card: sd-border-2 sd-shadow-sm

   pyEPR now ships a second HFSS transport through Ansys's official
   `PyAEDT <https://github.com/ansys/pyaedt>`_ library — **no COM, no pywin32**,
   works on Linux and macOS for the first time.

   .. code-block:: bash

      pip install "pyEPR-quantum[pyaedt]"

   .. code-block:: python

      from pyEPR.ansys_pyaedt import PyaedtDistributedAnalysis
      eprd = PyaedtDistributedAnalysis(pinfo, aedt_version="2026.1")
      eprd.do_EPR_analysis()   # pure gRPC — identical physics, cross-platform

   Physics and results are identical to the COM path — validated digit-for-digit.
   Contributed by `Joey Yaker <https://github.com/joeyyaker>`_.
   See :doc:`pyaedt_backend` for full details.

----

Quick install
=============

.. tab-set::

   .. tab-item:: pip

      .. code-block:: bash

         pip install pyEPR-quantum

   .. tab-item:: conda

      .. code-block:: bash

         conda install -c conda-forge pyepr-quantum

   .. tab-item:: uv

      .. code-block:: bash

         uv pip install pyEPR-quantum

No Ansys licence? Start with :ref:`without-hfss` or open
`Tutorial 6 on Binder <https://mybinder.org/v2/gh/zlatko-minev/pyEPR/master?filepath=_tutorial_notebooks%2FTutorial%206.%20EPR%20without%20HFSS%20%E2%80%94%20purely%20numerical%20workflow.ipynb>`_.

----

Five-line quickstart
====================

.. code-block:: python

   import pyEPR as epr

   pinfo = epr.ProjectInfo(project_path=r'C:\sims', project_name='my_chip',
                            design_name='qubit_cavity')
   pinfo.junctions['jQ'] = {'Lj_variable':'Lj1', 'rect':'junc_rect',
                             'line':'junc_line', 'Cj_variable':'Cj1'}
   eprd = epr.DistributedAnalysis(pinfo)
   eprd.do_EPR_analysis()
   epra = epr.QuantumAnalysis(eprd.data_filename)
   epra.analyze_all_variations(cos_trunc=8, fock_trunc=15)

See the :doc:`examples_quick` page for more complete examples including the
no-HFSS numerical workflow.

----

Contents
========

.. toctree::
   :maxdepth: 2

   about.rst
   installation.rst
   key_classes_reference.rst
   hfss_setup.rst
   examples_quick.rst
   without_hfss.rst
   pyaedt_backend.rst
   tutorials.rst
   troubleshooting.rst

.. toctree::
   :caption: API Reference
   :glob:
   :hidden:

   api/*


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
