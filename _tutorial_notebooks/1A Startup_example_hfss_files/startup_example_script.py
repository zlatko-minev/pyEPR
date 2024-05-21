# -*- coding: utf-8 -*-
"""
Example startup script to perform full quantization of a two qubit, one cavity Josephson circuit.
The results are saved, printed, and nicely plotted.

------~~~~!!!!------~~~~

Please also see the Jupyter notebook tutorials!

------~~~~!!!!------~~~~

@author: Zlatko
"""

from pyEPR import ProjectInfo, DistributedAnalysis, QuantumAnalysis

# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(
    "c:/sims",
    project_name="two_qubit_one_cavity",  # Project file name (string). "None" will get the current active one.
    design_name="Alice_Bob",  # Design name (string). "None" will get the current active one.
)

# 2a. Junctions. Specify junctions in HFSS model
project_info.junctions["jAlice"] = {
    "Lj_variable": "LJAlice",
    "rect": "qubitAlice",
    "line": "alice_line",
    "length": 0.0001,
}
project_info.junctions["jBob"] = {
    "Lj_variable": "LJBob",
    "rect": "qubitBob",
    "line": "bob_line",
    "length": 0.0001,
}

# 2b. Dissipative elements.
project_info.dissipative["dielectrics_bulk"] = [
    "si_substrate"
]  # supply names here, there are more options in project_info.dissipative.
project_info.dissipative["dielectric_surfaces"] = ["interface"]

# 3.  Run analysis
epr_hfss = DistributedAnalysis(project_info)
epr_hfss.do_EPR_analysis()

# 4.  Hamiltonian analysis
epr = QuantumAnalysis(epr_hfss.data_filename)
epr.analyze_all_variations(cos_trunc=8, fock_trunc=7)
epr.plot_hamiltonian_results()
