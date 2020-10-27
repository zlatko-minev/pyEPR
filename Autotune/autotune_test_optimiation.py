# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:24:30 2020

@author: Jules
"""

from autotune import Autotune
from collections import namedtuple
from pyEPR.project_info import ProjectInfo
import numpy as np





# Target = namedtuple('Target', ['Freqs', 'Qs', 'Kers', 'Chis'])
# t = Target(Freqs = {'transmon' : 5500,  'readout' : 6500, 'purcell' : 6500, 'memory' : 8000},
#            Qs = {'readout' : 30000, 'purcell' : 300},
#            Kers = {'transmon' : 200},
#            Chis = {('transmon', 'memory') : 2, ('transmon', 'readout') : 0.2}
#            )


pinfo = ProjectInfo("D:/Jules/HFSS_Jules/jules_project",
                    project_name='jules_project', 
                    design_name='autotune', 
                    setup_name='Setup1')
### The names used here ('transmon',..) must be exactly the same as thoses in 
### target_freqs etc...
# pinfo.junctions['transmon'] = {'rect' : 'ind_JJind_rect',
#                                                       'line' : 'ind_JJ_junction_line',
#                                                       'Lj_variable' : 'Lj'}
### The names used here ('memory',..) must be exactly the same as thoses in 
### target_freqs etc...
pinfo.resonators['memory'] = {'line' : 'mem_line'}
pinfo.resonators['readout'] = {'line' : 'ro_line'}
pinfo.resonators['purcell'] = {'line' : 'pur_line'}

pinfo.ports['port50ohm'] = {'rect':'port_50_purcell_track',
                                             'line': 'port_50_purcell_line',
                                             'R': 50}


### Must be the same names as in the pinfo.resonators and pinfo.junctions
target_Qs = {'readout' : 30000, 'purcell' : 300}
target_freqs = {'readout' : 6500, 'purcell' : 6500, 'memory' : 8000}
# target_freqs = {'transmon' : 5500,  'readout' : 6500, 'purcell' : 6500, 'memory' : 8000}
# target_kers = {'transmon' : 200}
target_kers = {}
target_chis = {}
# target_chis = {('transmon', 'memory') : 2, ('transmon', 'readout') : 0.2}


# weight_xxxx must have the same keys as target_xxxx
# weight_freqs = {'readout' : 1, 'purcell' : 1, 'memory' : 3}
# weight_Qs    = target_Qs = {'readout' : 2, 'purcell' : 2}
# weight_kers  = {}
# weight_chis  = {}





var_name = ['Lj', 'trm_length', 'capa_mem_length', 'capa_ro_length',
              'shift_capa_mem_trm_Y', 'coupling_dist', 'connector_coupling',
              'tune_ro']

# x = np.array([[7.18123456, 0.4778, 170, 170, 2500, 13, 520, 1.1],
#               [8.18123456, 0.4668, 180, 175, 2460, 13, 540, 1.1],
#               [9.18123456, 0.4668, 190, 175, 2460, 13, 540, 1.1],
#               [10.18123456, 0.4668, 200, 175, 2460, 13, 540, 1.1]])

x = np.array([[7.18123456, 0.4778, 170, 170, 2500, 13, 520, 1.1],
              [8.18123456, 0.4668, 180, 175, 2460, 13, 540, 1.1]])

a = Autotune(pinfo, target_freqs, target_Qs, 
                 target_kers, target_chis, var_name)





# from bayes_opt import BayesianOptimization, UtilityFunction 
# from bayes_opt import SequentialDomainReductionTransformer


# xybounds = {'x' : (-2, 2), 'y' : (-1,1)}

# niter = 30

# def debug_cost (x, y):
#     return -(x+y-2)**2-x**2-0*3*np.cos(5*(x+y))




# bounds_transformer = SequentialDomainReductionTransformer()
# optimizer = BayesianOptimization(f = debug_cost, pbounds=xybounds, verbose=2, 
#             random_state=1, bounds_transformer = bounds_transformer)

# suggested_points = []
# targets = []

# for i in range(niter) : 
#     optimizer = BayesianOptimization(f = debug_cost, pbounds=xybounds, verbose=2, 
#             random_state=1)
#     if len(suggested_points) != 0 :
#         optimizer.register(suggested_points[-1], targets[-1])
#     utility = UtilityFunction(kind="ucb", kappa = 2.576, xi=0.0)
#     suggested_points.append(optimizer.suggest(utility))
#     next_point = tuple(suggested_points[-1].values())
#     targets.append(debug_cost(*next_point))
          
      
            
            
            
