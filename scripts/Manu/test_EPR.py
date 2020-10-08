# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 19:31:33 2020

@author: Admin-local
"""

from pyEPR import *
from pyEPR.ansys import Optimetrics, HfssDesign
import numpy as np
import time
import matplotlib.pyplot as plt

################# 1.  Project and design. Open link to HFSS controls.
project_info = ProjectInfo(r'C:\HFSS_simu\\',
			     project_name = 'SMPD2', # Project file name (string). "None" will get the current active one.
			     design_name  = 'SMPD'       # Design name (string). "None" will get the current active one.
			    )


################# 2a. Junctions. Specify junctions in HFSS model
project_info.junctions['jtransmon'] = {'Lj_variable':'Jinduc', 'rect':'qubit_junction', 'line': 'qubit_junction_line', 'length':5e-6}
#




################ 2c. Define ports.

project_info.ports['Waste'] = {'rect':'WasP_connector_ohm', 'R': 50, 'line': 'WasP_connector_line'}
project_info.ports['Buffer'] = {'rect':'BufP_connector_ohm', 'R': 50,  'line': 'BufP_connector_line'}
project_info.ports['Qubit'] = {'rect':'qubit_connector_ohm', 'R': 50, 'line': 'qubit_connector_line'}
project_info.ports['Flux'] = {'rect':'BufR_squid_loop_lumped', 'R': 50, 'line': 'BufR_squid_connector_line'}

        ###load the optimetrics module from the ansys package
    
    ################# 3 - performs a pyEPR analysis on the new 'm' variations
    ###reload the list of the last variations
epr_hfss = DistributedAnalysis(project_info)
epr_hfss.do_EPR_analysis()

epr = QuantumAnalysis(epr_hfss.data_filename)
epr.analyze_all_variations(cos_trunc = 5, fock_trunc = 4)

    
 
    
        
var_list=['0']
var=0

### get the frequencies of the current variation
chi_dic=epr.results.get_chi_O1()
chis=np.abs(np.array(chi_dic[var_list[var]]))
### get the frequencies of the current variation
freq_dic=epr.results.get_frequencies_O1()
freq=np.abs(np.array(freq_dic[var_list[var]]))
### get the anharmonicity of the current variation
anharmonicity=np.diag(chis)
### get the Q of the current variation
total_Q_from_HFSS = np.array(epr.Qs)[:,var]
total_Q_from_couplings = 1/(1/np.array(epr.Qm_coupling[var_list[var]])).sum(1)
Q_couplings_adjusted=np.array([total_Q_from_HFSS/total_Q_from_couplings]).T*np.array(epr.Qm_coupling[str(var)])

Q_couplings_adjusted=np.array(epr.Qm_coupling[var_list[var]])
