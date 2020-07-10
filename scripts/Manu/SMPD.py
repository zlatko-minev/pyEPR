# -*- coding: utf-8 -*-
"""
Script for analysing SMPD lambdaover2:
C:\GitHub\Quantrolib\drawings\Léo\SMPD_lambda2.py
"""

from pyEPR import *


# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(r'C:\HFSS_simu\\',
			     project_name = 'SMPD', # Project file name (string). "None" will get the current active one.
			     design_name  = 'fullcircuit_final'       # Design name (string). "None" will get the current active one.
			    )




# 2a. Junctions. Specify junctions in HFSS model
project_info.junctions['jtransmon'] = {'Lj_variable':'Jinduc', 'rect':'qubit_junction', 'line': 'qubit_junction_line', 'length':5e-6}
#

project_info.ports['Waste'] = {'rect':'WasP_connector_ohm', 'R': 50, 'line': 'WasP_connector_line'}
project_info.ports['Buffer'] = {'rect':'BufP_connector_ohm', 'R': 50,  'line': 'BufP_connector_line'}
project_info.ports['Qubit'] = {'rect':'qubit_connector_ohm', 'R': 50, 'line': 'qubit_connector_line'}
project_info.ports['Flux'] = {'rect':'BufR_squid_loop_lumped', 'R': 50, 'line': 'BufR_squid_connector_line'}



## 2b. Dissipative elements.
project_info.dissipative['dielectrics_bulk']    = ['silicon_substrate']    # supply names here, there are more options in project_info.dissipative.
project_info.dissipative['dielectric_surfaces'] = ['silicon_surface']

# 3.  Run analysis
epr_hfss = DistributedAnalysis(project_info)
epr_hfss.do_EPR_analysis()

# 4.  Hamiltonian analysis
epr      = QuantumAnalysis(epr_hfss.data_filename)
epr.analyze_all_variations(cos_trunc = 5, fock_trunc = 4)
epr.plot_hamiltonian_results()

