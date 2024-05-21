# -*- coding: utf-8 -*-
"""
My First pyEPR Script
"""

from pyEPR import *

# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(
    r"C:\zkm\my-first-pyEPR\\",
    project_name="HelloWorld-pyEPR",  # Project file name (string). "None" will get the current active one.
    design_name="MyFirstTest",  # Design name (string). "None" will get the current active one.
)

project_info.connect_to_project()


#
## 2a. Junctions. Specify junctions in HFSS model
# project_info.junctions['jAlice'] = {'Lj_variable':'LJAlice', 'rect':'qubitAlice', 'line': 'alice_line', 'length':0.0001}
# project_info.junctions['jBob']   = {'Lj_variable':'LJBob',   'rect':'qubitBob',   'line': 'bob_line',   'length':0.0001}
#
## 2b. Dissipative elements.
# project_info.dissipative['dielectrics_bulk']    = ['si_substrate']    # supply names here, there are more options in project_info.dissipative.
# project_info.dissipative['dielectric_surfaces'] = ['interface']
#
## 3.  Run analysis
# epr_hfss = DistributedAnalysis(project_info)
# epr_hfss.do_EPR_analysis()
#
## 4.  Hamiltonian analysis
# epr      = QuantumAnalysis(epr_hfss.data_filename)
# epr.analyze_all_variations(cos_trunc = 8, fock_trunc = 7)
# epr.plot_hamiltonian_results()
#
