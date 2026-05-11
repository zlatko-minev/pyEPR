.. _install:

**************
Installation
**************

.. contents:: On this page
   :local:
   :depth: 2

Requirements
============

- **Python** 3.9 – 3.12
- **Operating system:** Windows, macOS, or Linux — see `Platform notes`_ below.

.. _install-quick:

Quick install
=============

.. tabs::

   .. tab:: uv *(recommended)*

      `uv <https://github.com/astral-sh/uv>`_ is a fast, modern Python package manager.
      Install it once, then:

      .. code-block:: bash

         uv pip install pyEPR-quantum

      Or, inside a uv project/virtual environment:

      .. code-block:: bash

         uv add pyEPR-quantum

   .. tab:: pip + venv

      .. code-block:: bash

         python -m venv .venv
         source .venv/bin/activate   # Windows: .venv\Scripts\activate
         pip install pyEPR-quantum

   .. tab:: conda

      ``pyEPR-quantum`` is available on the ``conda-forge`` channel:

      .. code-block:: bash

         conda create -n pyepr python=3.11
         conda activate pyepr
         conda install -c conda-forge pyepr-quantum

      .. note::

         The conda-forge package name is ``pyepr-quantum`` (lower-case).
         The PyPI name is ``pyEPR-quantum``.  Either way, you import it as
         ``import pyEPR as epr``.

.. _install-dev:

Development / editable install
================================

Clone the repository and install in editable mode so that your local changes
are picked up immediately:

.. tabs::

   .. tab:: uv *(recommended)*

      .. code-block:: bash

         git clone https://github.com/zlatko-minev/pyEPR.git
         cd pyEPR
         uv pip install -e ".[test]"

   .. tab:: pip + venv

      .. code-block:: bash

         git clone https://github.com/zlatko-minev/pyEPR.git
         cd pyEPR
         python -m venv .venv
         source .venv/bin/activate   # Windows: .venv\Scripts\activate
         pip install -e ".[test]"

   .. tab:: conda

      .. code-block:: bash

         git clone https://github.com/zlatko-minev/pyEPR.git
         cd pyEPR
         conda create -n pyepr python=3.11
         conda activate pyepr
         pip install -e ".[test]"

The ``[test]`` extra installs ``pytest`` and ``pytest-cov``.  Run the test
suite (no Ansys required):

.. code-block:: bash

   pytest

.. _install-platform:

Platform notes
==============

pyEPR has two distinct functional layers with different platform requirements:

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Feature
     - Windows
     - macOS
     - Linux
   * - EPR / quantum analysis
     - ✅ Full support
     - ✅ Full support
     - ✅ Full support
   * - Ansys HFSS COM interface
     - ✅ Full support
     - ⚠️ Limited (see below)
     - ⚠️ Limited (see below)

**EPR and quantum analysis** (``DistributedAnalysis``, ``QuantumAnalysis``) are
pure-Python and work identically on all platforms — no Ansys installation
required for post-processing existing data.

**Ansys HFSS COM interface** (``ansys.py``, ``HfssDesign``, etc.) was originally
written for Windows, where HFSS exposes a COM/DCOM automation interface via
``pythoncom`` / ``win32com``.  On macOS and Linux, Ansys ships a non-Windows
scripting interface (IronPython inside the AEDT product, or the ``ansys-aedt``
Python API); pyEPR's COM wrapper layer does not use these interfaces directly.
If you are running HFSS on a remote Windows machine or a Windows VM, you can
still drive it from macOS/Linux over the network COM bridge.

.. note::

   **PyAEDT** — Ansys's official open-source Python library
   (`github.com/ansys/pyaedt <https://github.com/ansys/pyaedt>`_) — provides
   a fully cross-platform scripting layer for AEDT.  pyEPR predates PyAEDT
   by several years and covers the quantum-circuit EPR analysis workflow that
   PyAEDT does not.  For new workflows that only need simulation control
   (geometry, mesh, solve) and not EPR quantization, PyAEDT is a natural
   complement or alternative for the COM layer.

.. _install-hfss-version:

Ansys HFSS version compatibility
=================================

AEDT 2021.2 renamed the strings returned by ``GetSolutionType()``.
AEDT 2024.1 changed the default design type created by ``InsertDesign``.
pyEPR ≥ 0.9.2 handles both transparently — always use the latest pyEPR.

For the full technical details see the :ref:`HFSS compatibility section
<install-hfss-version>` in the developer guide (``CLAUDE.md``).

.. _install-config:

Optional configuration
======================

After installation, you can customise pyEPR's behaviour by editing
``pyEPR/_config_user.py``.  The most useful settings are the default
data-save directory and logging verbosity.

To prevent git from tracking your local config changes:

.. code-block:: bash

   git update-index --skip-worktree pyEPR/_config_user.py
