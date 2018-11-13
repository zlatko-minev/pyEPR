# -*- coding: utf-8 -*-
"""
Example bare bones script to perform full quantization of a two qubit Josephson circuit.
@author: Zlatko
"""

from pyEPR import Project_Info, pyEPR_HFSS, pyEPR_Analysis

# 1.  Project and design. Open link to HFSS controls.
project_info = Project_Info('c:/sims',
                            project_name = 'two_qubit_one_cavity', # Project file name (string). "None" will get the current active one.
                            design_name  = 'Alice_Bob'             # Design name (string). "None" will get the current active one.
                            )

# 2a. Junctions. Specify junctions in HFSS model
project_info.junctions['jAlice'] = {'Lj_variable':'LJAlice', 'rect':'qubitAlice', 'line': 'alice_line', 'length':0.0001}
project_info.junctions['jBob']   = {'Lj_variable':'LJBob',   'rect':'qubitBob',   'line': 'bob_line',   'length':0.0001}

# 2b. Dissipative elements.
project_info.dissipative.dielectrics_bulk    = ['si_substrate']    # supply names here, there are more options in project_info.dissipative.
project_info.dissipative.dielectric_surfaces = ['interface']

# 3.  Run analysis
epr_hfss = pyEPR_HFSS(project_info)
epr_hfss.do_EPR_analysis()

# 4.  Hamiltonian analysis
epr      = pyEPR_Analysis(epr_hfss.data_filename)
epr.analyze_all_variations(cos_trunc = 8, fock_trunc = 7)
epr.plot_Hresults()