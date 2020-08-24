# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 11:26:36 2020

@author: berdou
"""

import numpy as np
from HFSSdrawpy.utils import val
from HFSSdrawpy import Modeler, Body
from drawpylib.parameters import TRACK, GAP, RLC, MESH
import drawpylib.cpw_elements as elt
from pyEPR.ansys import Optimetrics, HfssDesign
from pyEPR.project_info import ProjectInfo
from pyEPR import DistributedAnalysis, QuantumAnalysis
from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer
import time



    
new_geometry = False
if new_geometry : 
    pm = Modeler('hfss')
    chip1 = Body(pm, 'chip1')
    
    chip_width= pm.set_variable("2mm")
    chip_length=pm.set_variable("2mm")
    chip_thickness= pm.set_variable("280um")
    
    
    #Transmon
    track_trm_ro = pm.set_variable('42um')
    gap_trm_ro = pm.set_variable('25um')
    Lj = pm.set_variable('12nH')
    trm_length = pm.set_variable('0.6mm')
    
    trm_posx = pm.set_variable('1mm')
    trm_posy = pm.set_variable('1mm')
    
    with chip1([trm_posx,trm_posy],[1,0]):  #+c.short_trm+c.trm_length/2
        elt.draw_half_transmon(chip1, track_trm_ro, gap_trm_ro, trm_length, '5um',Lj, name='transmon') #no port
     
    
    ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
    ground_plane.subtract(chip1.entities[GAP])
    ground_plane.unite(chip1.entities[TRACK])
    ground_plane.assign_perfect_E()
    
    #chip substrate
    chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness], material='silicon', name='substrate')
    chip1.box([0,0,0],[chip_width, chip_length, 3*chip_thickness], name='vaccuum')
    
    
    path = r'D:\Users\hqc\Documents\Jules\drawpy_camille_jules'
    pm.generate_gds(path, 'transmon_seul')


#TUNING
#def is_st = len(string_list)


def timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name


dir_path = r"D:\Users\hqc\Documents\Jules\HFSS"
project_info = ProjectInfo(r"D:\Users\hqc\Documents\Jules\\",
                               project_name  = 'Exemple_tomograpy',  # Name of the project file (string). "None" will get the current active one.
                               design_name   = 'HFSSDesign20',  # Name of the desgin file (string). "None" will get the current active one.
                               setup_name = 'Setup1')
project_info.junctions['jtransmon'] = {'rect':'ind_JJind_rect',  'line': 'ind_JJ_junction_line', 'Lj_variable':'Lj'}


target = np.array([5500, 200])  #freq and ker
Lj = 8.42440852233249
trm_length = 0.50111


tested_inputs = []
calculated_outputs = []

def cost_camille (Lj, trm_length):
    if [Lj, trm_length] in tested_inputs :
        freq = calculated_outputs[tested_inputs.index([Lj, trm_length])][0]
        chi = calculated_outputs[tested_inputs.index([Lj, trm_length])][1]
    else : 
        tested_inputs.append([Lj, trm_length])
        x0 = np.array([trm_length, Lj])
        epr_hfss = DistributedAnalysis(project_info)
        ###load the optimetrics module from the ansys package
        opti=Optimetrics( epr_hfss.design)
        
        ###load the optimetrics module from the ansys package
        parametric_name=timestamp_name('parametric')
        print(parametric_name)
        
        ###get the list of HFSS variable
        var=epr_hfss.design.get_variables()   
        
        ###list of variable to be optimized (should probably done outside)
        name=np.array(["trm_length","Lj"])
        
        ###dirty way to get the unit of the variable of interest (for some reason HFSS wants the same units than the project variable for the variation)
        units=[var[key][-2:] for key in name]
        
        ###dirty way to add the correct units to the 'x' array
        var_list=np.ones(x0.shape[0]).astype(str)
        for i in range(x0.shape[0]):
                var_list[i]=str(x0[i])+units[i]
        print(var_list)
                
        ###create and save the array of string with the correct format to import to HFSS
        arr=np.vstack([name,var_list])
        np.savetxt(dir_path+"\%s.txt"%parametric_name, arr, fmt='%s', delimiter='; ', newline='\n', header='', footer='', comments='# ', encoding=None)
        
        ###import the parametric setup with the list of variations (function I added to the ansys package)
        opti.import_setup(parametric_name,dir_path+"\%s.txt"%parametric_name)
        
        opti.solve_setup(parametric_name)
        
        # Reload the project with the new variation 
        
        epr_hfss = DistributedAnalysis(project_info)
        var_list=list(epr_hfss.variations)
        epr_hfss.do_EPR_analysis([var_list[-1]])
       
        epr = QuantumAnalysis(epr_hfss.data_filename, [var_list[-1]])
        epr.analyze_all_variations([var_list[-1]], cos_trunc = 5, fock_trunc = 4)
        
        ### Idea to check which variation was prevously done 
    #    string_list = ["Lj="+"'"+str(Lj)+'nH'+"'",
    #                   "trm_length="+"'"+str(trm_length)+'mm'+"'"]
    #    var_index = 0
    #    string = epr_hfss.get_variation_string
    #    for ii in range(len(var_list)): 
    #        if is_string_list_in_string(string_list, string(str(ii))) : 
    #            var_index = ii
                
        
        chi_dic = epr.results.get_chi_O1()
        chis = np.abs(np.array(chi_dic[var_list[-1]]))
        freq_dic = epr.results.get_frequencies_O1()
        freq = np.abs(np.array(freq_dic[var_list[-1]]))[0]
        chi = np.diag(chis)[0]    
        
        calculated_outputs.append([freq, chi])
        
    print('saved_freq', freq)
    print('saved_chi', chi)

    return(-(((freq-target[0])/target[0])**2 + ((chi-target[1])/target[1])**2))

pbounds = {'trm_length': (0.3, 0.95), 'Lj': (5.9, 15.9)} 
bounds_transformer = SequentialDomainReductionTransformer()
## Reduces the bounds ###
domred_optimizer = BayesianOptimization(f = cost_camille, pbounds=pbounds, verbose=2, 
                                          random_state=1, bounds_transformer = bounds_transformer)
 


domred_optimizer.maximize(10, 20)

##  Lj = 8.407
## trm_length = 0.4575