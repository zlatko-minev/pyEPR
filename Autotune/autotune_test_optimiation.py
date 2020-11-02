# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:24:30 2020

@author: Jules
"""

from autotune import Autotune
from collections import namedtuple
from pyEPR.project_info import ProjectInfo
import numpy as np








pinfo = ProjectInfo("D:/Jules/HFSS_Jules/Spice_project",
                    project_name='Spice_project', 
                    design_name='mem_trm_ro_asym', 
                    setup_name='Setup1')
### The names used here ('transmon',..) must be exactly the same as thoses in 
# ### target_freqs etc...
pinfo.junctions['transmon'] = {'rect' : 'ind_JJind_rect',
                                'line' : 'ind_JJ_junction_line',
                                'Lj_variable' : 'Lj'}

# ### The names used here ('memory',..) must be exactly the same as thoses in 
### target_freqs etc...
pinfo.resonators['memory'] = {'line' : 'mem_line'}
pinfo.resonators['readout'] = {'line' : 'ro_line'}

   


### Must be the same names as in the pinfo.resonators and pinfo.junctions
kappa_pur = 200
Q_pur = 7500/kappa_pur


target_Qs = {}
target_freqs = {'memory' : 7500, 'readout' : 7500, 'transmon' : 6000}
target_kers = {'transmon' : 20}
target_chis = {('transmon', 'memory') : 0.2, ('transmon', 'readout') : 0.2}






var_name = ['Capa_U_mem']

# x = np.array([[7.18123456, 0.4778, 170, 170, 2500, 13, 520, 1.1],
#               [8.18123456, 0.4668, 180, 175, 2460, 13, 540, 1.1],
#               [9.18123456, 0.4668, 190, 175, 2460, 13, 540, 1.1],
#               [10.18123456, 0.4668, 200, 175, 2460, 13, 540, 1.1]])

x = np.array([[200]])

a = Autotune(pinfo, target_freqs, target_Qs, 
                 target_kers, target_chis, var_name)


pbounds = {'tune1' : (200, 2200), 'tune2' : (800, 2500), 'Xmon_gnd2' : (1,30)}



        
