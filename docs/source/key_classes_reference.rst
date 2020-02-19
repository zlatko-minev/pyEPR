Main classes
=================================

The first main class of pyEPR is :ref:`project-info`, which instansiates and stores the Ansys interfaces classes and user-defined parameters related to the design, such as junction names and properties.

The second main class of pyEPR is :ref:`distributed-analysis`, which performs the EPR analysis on the ansys eigenfield solutions from the fields. It saves the calculated energy participation ratios (EPRs) and realted convergences, and other paramete results. It does not calculate the Hamiltonian.
This is left for the third class.

The third main class of pyEPR is :ref:`quantum-analysis`, which uses the EPRs and other save quantities to create and diagonalizae the Hamiltonian.

.. _project-info:

ProjectInfo
-----------

.. autoclass:: pyEPR.project_info.ProjectInfo


.. _distributed-analysis:

DistributedAnalysis
----------------------

.. autoclass:: pyEPR.core_distributed_analysis.DistributedAnalysis

.. _quantum-analysis:

QuantumAnalysis
----------------------
.. autoclass:: pyEPR.core_quantum_analysis.QuantumAnalysis