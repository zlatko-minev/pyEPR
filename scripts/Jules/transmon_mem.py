# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 11:26:36 2020

@author: berdou
"""

import numpy as np
from HFSSdrawpy.utils import val
from HFSSdrawpy.core.modeler import Modeler
from HFSSdrawpy.core.body import Body
from drawpylib.parameters import TRACK, GAP, RLC, MESH
import drawpylib.cpw_elements as elt
from pyEPR.project_info import ProjectInfo
from pyEPR.ansys import Optimetrics
from pyEPR import DistributedAnalysis, QuantumAnalysis
import time
from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer


new_geometry = False
if new_geometry : 
    pm = Modeler('hfss')
    relative = pm.set_variable('1mm')
    chip1 = Body(pm, 'chip1')
    
    is_bond = True
    track = pm.set_variable('20um')
    chip_width= pm.set_variable("6mm")
    chip_length=pm.set_variable("6mm")
    chip_thickness= pm.set_variable("280um")    
    pcb_thickness= pm.set_variable("320um")
    vaccuum_thickness= pm.set_variable(6*chip_thickness)
    
    
    #### Connectors
    if True:
        track= pm.set_variable("42um")
        gap = pm.set_variable("25um")
        pcb_track, pcb_gap, bond_length, bond_slope = '540um','210um','200um', '0.5'
        con = pm.set_variable('3.2mm')
        
    #    with chip1([0, chip_length/2-con], [1,0]):
    #        in_top_L_port, = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length)
    #    with chip1([0, chip_length/2+con], [1,0]):
    #        in_top_R_port, = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length)
    #    with chip1([chip_width, chip_length/2-con], [-1,0]):
    #        in_bot_L_port, = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length)
    #    with chip1([chip_width, chip_length/2+con], [-1,0]):
    #        in_bot_R_port, = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length)#pour faire une petite simu, on peut mettre un draw_end_cable à la place du port 
    #    with chip1([chip_width/2-con, 0], [0,1]):
    #        in_L, = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length)
    #    with chip1([chip_width/2-con, chip_length], [0,-1]):
    #        in_R, = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length)
    
    
    fillet =pm.set_variable('200um')
    # Capa buff mem
    dist_buff_mem = pm.set_variable('5um')
    track_capa = pm.set_variable('42um')
    gap_capa = pm.set_variable('25um')
    capa_length = pm.set_variable('1000um')
    #capa U mem trm
    mem_meander_length = pm.set_variable('450um')
    mem_meander_offset = pm.set_variable('-100um')
    shift_capa_mem_trm_Y = pm.set_variable('-2500um')
    shift_capa_mem_trm_X = pm.set_variable('-3500um')
    fillet_mem = pm.set_variable('200um')
    track_mem = pm.set_variable('42um')
    gap_mem = pm.set_variable('25um')
    capa_ro_length = pm.set_variable('136.7um')
    capa_U_length = pm.set_variable('158.7um')
    short_trm = pm.set_variable('5um')
    shift_const_X = pm.set_variable('800um')
    shift_const_Y = pm.set_variable('1900um')
    #Transmon
    track_trm_ro = pm.set_variable('42um')
    gap_trm_ro = pm.set_variable('25um')
    Lj = pm.set_variable('8.5nH')
    T_pad_width = pm.set_variable('0.042um')
    T_pad_spacing = pm.set_variable('0.025um')
    trm_length = pm.set_variable('0.465mm')
    T_fillet = pm.set_variable('0.005mm')
    #constrain port for readout 
    tune_ro = pm.set_variable('1.09mm')
    pos_constr_ro = pm.set_variable('3mm') 
    # U capas
    width_U_readout = 2*(gap_trm_ro + track_trm_ro)+2*short_trm + track_trm_ro+2*gap_trm_ro
    width_U_mem = 2*(gap_mem + track_mem)+2*short_trm + track_trm_ro+2*gap_trm_ro#p-e mettre 2 short-trm differents pr la mem et le readout (2 U differents)       
    #Purcell
    coupling_dist = pm.set_variable('13um')
    connector_coupling = pm.set_variable('540um')
    purcell_port_connect = pm.set_variable('500um')
    delta_ro_purc = pm.set_variable('-95um')
    #Position of the memory that will dict everyone's position
    pos_x_mem = pm.set_variable('0.5mm')
    pos_y_mem = pm.set_variable('3.75mm')
    
    with chip1([pos_x_mem, pos_y_mem],[0,-1]):
            T_portL, T_portR, T_port3,= elt.draw_T(chip1, track_capa, gap_capa)
            T_mesh = chip1.rect([-track_capa, - gap_capa], [2*track_capa, track_capa + 2*gap_capa], name = 'T_mesh', layer = MESH)
            T_mesh.assign_mesh_length(track_capa/3)
            with chip1([capa_length/2,0],[1,0]):
                capa_end_L, = elt.draw_end_cable(chip1, track_capa,gap_capa,typeEnd='open', name = 'capa_end_L')
            with chip1([-capa_length/2,0],[-1,0]):
                capa_end_R, = elt.draw_end_cable(chip1, track_capa, gap_capa, typeEnd = 'open', name = 'capa_end_R')
            chip1. draw_cable(T_portL, capa_end_L,
                         is_bond=is_bond, fillet=fillet, name ='capa_L')
            chip1.draw_cable(T_portR, capa_end_R, 
                             is_bond= is_bond, fillet= fillet, name = 'capa_R')
            
            #capa U mem trm
            with chip1([-shift_capa_mem_trm_Y,-shift_capa_mem_trm_X],[-1,0]):
                Xmon_mem, = elt.draw_Xmon_end_cable(chip1, track_mem, gap_mem, width_U_mem, capa_U_length, name ='U_coupler' )
            #constrain 
                with chip1([-shift_capa_mem_trm_Y,shift_const_X],[0,1]): #[shift_const_X,-shift_capa_mem_trm_Y
                    constrain_mem, = elt.create_port(chip1, name = 'constrain_mem')    
              # Transmon
                with chip1([2*gap_mem+track_mem+trm_length/2+gap_trm_ro+short_trm,0],[0,1]):  #+c.short_trm+c.trm_length/2
                    elt.draw_half_transmon(chip1, track_trm_ro, gap_trm_ro, trm_length, '5um',Lj, name='transmon') #no port
           
    # draw mem cable
    chip1.draw_cable(T_port3, constrain_mem, Xmon_mem, is_bond = is_bond,
                                 fillet = fillet_mem, name = 'mem')
    
    
    # Meshing 
    ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
    ground_plane.subtract(chip1.entities[GAP])
    ground_plane.unite(chip1.entities[TRACK])
    ground_plane.assign_perfect_E()
    
    #chip substrate
    chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness], material='silicon', name='substrate')
    chip1.box([0,0,0],[chip_width, chip_length, 6*chip_thickness], name='vaccuum')
    
    
    path = r'D:\Users\hqc\Documents\Jules\Python\drawpy_camille_jules\K_layout'
    pm.generate_gds(path, 'trm_mem')
    
 




    #######################
##########    TUNING     ##########
    ########################
    ### Before tuning, make sure that you have ran one converging HFSS simulation with the right 
    ### setup parameters 

def timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name


dir_path = r"D:\Users\hqc\Documents\Jules\HFSS"
project_info = ProjectInfo(r"D:\Users\hqc\Documents\Jules\\",
                               project_name  = 'Exemple_tomograpy',  
                               design_name   = 'HFSSDesign8', 
                               setup_name = 'Setup1')
project_info.junctions['jtransmon'] = {'rect':'ind_JJind_rect',  'line': 'ind_JJ_junction_line', 'Lj_variable':'Lj'}



target = np.array([5500, 8000,200, 0.3])  #freq_tr / freq_mem / ker_tr / chi_tr_mem
# good guesses (already optimized)
Lj = 8.20073 
trm_length = 0.437058092535451 #mm
shift_capa_mem_trm_Y = -2328.027 #um
capa_U_length = 136.7304630287759 #um


trm_mem_inputs = []
trm_mem_outputs = []

def cost2_camille (Lj, trm_length, shift_capa_mem_trm_Y, capa_U_length):
    x0 = [trm_length, Lj, shift_capa_mem_trm_Y, capa_U_length]
    if x0 in trm_mem_inputs :
        freq_mem = trm_mem_outputs[trm_mem_inputs.index(x0)][0]
        freq_trm = trm_mem_outputs[trm_mem_inputs.index(x0)][1]   
        chi_trm = trm_mem_outputs[trm_mem_inputs.index(x0)][2]
        chi_mem_trm = trm_mem_outputs[trm_mem_inputs.index(x0)][3]
    
    else : 
        trm_mem_inputs.append(x0)
        x0 = np.array(x0)
        epr_hfss = DistributedAnalysis(project_info)
        
        ###OPTIMETRICS
        opti = Optimetrics( epr_hfss.design)
        parametric_name=timestamp_name('parametric')
        print(parametric_name)
        
        ###VARIABLES
        var=epr_hfss.design.get_variables()   
        name=np.array(["trm_length","Lj","shift_capa_mem_trm_Y", "capa_U_length"])
        
        ###UNITS
        units=[var[key][-2:] for key in name]    
        var_list=np.ones(x0.shape[0]).astype(str)
        for i in range(x0.shape[0]):
                var_list[i]=str(x0[i])+units[i]
        print(var_list)
                
        ###SAVE 
        arr=np.vstack([name,var_list])
        np.savetxt(dir_path+"\%s.txt"%parametric_name, arr, fmt='%s', delimiter='; ', newline='\n', header='', footer='', comments='# ', encoding=None)
           
        ###Import to HFSS and solve
        opti.import_setup(parametric_name,dir_path+"\%s.txt"%parametric_name)
        opti.solve_setup(parametric_name)
            
        ### RELOAD 
        epr_hfss = DistributedAnalysis(project_info)
        var_list=list(epr_hfss.variations)
        ### 
        epr_hfss.do_EPR_analysis([var_list[-1]])
    
        epr = QuantumAnalysis(epr_hfss.data_filename,[var_list[-1]])
        epr.analyze_all_variations([var_list[-1]],cos_trunc = 5, fock_trunc = 4)
        
        freqs=np.array(epr.get_frequencies()).T   
        chi_dic = epr.results.get_chi_O1()
        chis = np.abs(np.array(chi_dic[var_list[-1]]))
        freq_dic = epr.results.get_frequencies_O1()
        freq = np.abs(np.array(freq_dic[var_list[-1]]))
        chi = np.diag(chis)   
        
        ### sorting modes  
        ### define the qubit as the mode with the largest anharmanocity
        index={}
        index['trm']=np.argsort(chi)[-1]
        index['mem']=np.argsort(chi)[-2]
            
        
        chi_mem_trm = chis[index['mem']][index['trm']]
        chi_trm = chi[index['trm']]
        freq_trm = freq[index['trm']]
        freq_mem = freq[index['mem']]
        trm_mem_outputs.append([freq_mem, freq_trm, chi_trm, chi_mem_trm])
  
    print('freq_mem', freq_mem)
    print('freq_trm', freq_trm)
    print('chi_trm', chi_trm)
    print('chi_mem_trm', chi_mem_trm)
    return (-(((freq_trm-target[0])/target[0])**2 + ((freq_mem-target[1])/target[1])**2
              +((chi_trm-target[2])/target[2])**2 + ((chi_mem_trm-target[3])/target[3])**2))

    


##### BAYESIAN OPTIMIZATION ###
######The bounds are given by the previous codes (transmon alone and memory alone)
pbounds = {'Lj': (7, 10), 'trm_length': (0.35, 0.52), 'shift_capa_mem_trm_Y': (-2500, -2200),
           'capa_U_length' : (155, 170) } 
bounds_transformer = SequentialDomainReductionTransformer()
### Reduces the bounds ###
domred_optimizer = BayesianOptimization(f = cost2_camille,pbounds=pbounds, verbose=2, 
                                          random_state=1, bounds_transformer = bounds_transformer)

#domred_optimizer.maximize(10,50)

#When the optimizer is running he converges at f_trm = 5000 GHZ and the other values (freq_mem, 
#chis..) are extremely close rome the target