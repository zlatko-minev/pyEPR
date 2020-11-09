# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:24:30 2020

@author: Jules
"""

from autotune import Autotune
from pyEPR.project_info import ProjectInfo
import numpy as np




##### 1. Inquire the project info
pinfo = ProjectInfo("D:/Jules/HFSS_Jules/Spice_project",
                    project_name='Spice_project', 
                    design_name='autotune_test_design', 
                    setup_name='Setup1')

# The names used here ('transmon',..) must be exactly the same as thoses in 
# targ_freqs etc....
pinfo.junctions['transmon'] = {'rect' : 'ind_JJ1ind_rect',
                                'line' : 'ind_JJ1_junction_line',
                                'Lj_variable' : 'Lj'}

#  The names used here ('memory',..) must be exactly the same as thoses in 
# target_freqs etc...
pinfo.resonators['memory'] = {'line' : 'mem_line'}
pinfo.resonators['readout'] = {'line' : 'ro_line'}

   

##### 2. Inquire the targets you wish to obtain (units : MHZ)

# Must be the same names as in the pinfo.resonators and pinfo.junctions

target_freqs = {'transmon' : 5000, 'memory' : 8000, 'readout' : 7500, }
target_Qs = {'readout' : 6000}
target_kers = {'transmon' : 200}
target_chis = {('transmon', 'memory') : 2}

# Eventually inquire the weights : 
    # weight_freqs = {'transmon' : 2, 'memory' : 1, 'readout' : 3}
    # weight_Qs    = target_Qs = {'readout' : 2, 'purcell' : 2}
    # weight_kers  = {'transmon' : 2}
    # weight_chis  = {('transmon', 'memory') : 1}


##### 3. Inquire the name of the HFSS variable that will be changed
var_name = ['tune_mem', 'tune_ro', 'Lj',
            'cpl_ro_50ohm', 'trm_length', 'cpl_mem_trm']

# If you want to test the functions, you can set an x array to change 
# manually the value of the variables :
x = np.array([[950, 2000, 11, 200, 450, 220]])

##### 4. Launch Autotune !
auto = Autotune(pinfo, target_freqs, target_Qs, 
                 target_kers, target_chis, var_name)


##### 5. Try to get sorted modes : 
# auto.get_sorted_modes(x)
    
##### 6. Optimize !
# Correct syntax of bounds for a bayesian optimization
# The name must be the same as the var_name
pbounds = {'tune_mem' : (300, 2700), 'tune_ro' : (650, 2800), 'Lj' : (5,20),
           'cpl_ro_50ohm' : (100, 400), 'trm_length' : (300, 600),
           'cpl_mem_trm' : (100,400)}
# auto.optimize(pbounds)



        
