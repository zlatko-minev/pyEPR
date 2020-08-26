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
from pyEPR.ansys import Optimetrics
import time
from pyEPR.project_info import ProjectInfo
from pyEPR import DistributedAnalysis, QuantumAnalysis
from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer

new_geometry = False
if new_geometry : 
    pm = Modeler('hfss')
    relative = pm.set_variable('1mm')
    chip1 = Body(pm, 'chip1')
    
    is_bond = False
    track = pm.set_variable('20um')
    chip_width= pm.set_variable("3.5mm")
    chip_length=pm.set_variable("4.5mm")
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
    shift_capa_mem_trm_Y = pm.set_variable('-1500um')
    shift_capa_mem_trm_X = pm.set_variable('-3500um')
    fillet_mem = pm.set_variable('200um')
    track_mem = pm.set_variable('42um')
    gap_mem = pm.set_variable('25um')
    capa_ro_length = pm.set_variable('136.7um')
    capa_U_length = pm.set_variable('200um')
    short_trm = pm.set_variable('5um')
    shift_const_X = pm.set_variable('800um')
    shift_const_Y = pm.set_variable('1900um')
    #Transmon
    track_trm_ro = pm.set_variable('42um')
    gap_trm_ro = pm.set_variable('25um')
    Lj = pm.set_variable('12nH')
    T_pad_width = pm.set_variable('0.042um')
    T_pad_spacing = pm.set_variable('0.025um')
    trm_length = pm.set_variable('0.6mm')
    T_fillet = pm.set_variable('0.005mm')
    #constrain port for readout 
    tune_ro = pm.set_variable('1.257mm')
    pos_constr_ro = pm.set_variable('3mm') #distance from transmon to curve of readout cable
    width_U_readout = 2*(gap_trm_ro + track_trm_ro)+2*short_trm + track_trm_ro+2*gap_trm_ro
    
    pos_x_mem = pm.set_variable('-2mm')
    pos_y_mem = pm.set_variable('-1mm')
    with chip1([pos_x_mem, pos_y_mem],[0,-1]):
        with chip1([-shift_capa_mem_trm_Y,-shift_capa_mem_trm_X],[-1,0]):   
              with chip1([2*gap_mem+track_mem+trm_length/2+gap_trm_ro+short_trm,0],[0,1]):  #+c.short_trm+c.trm_length/2
                   #purcell filter
                    coupling_dist = pm.set_variable('10um')
                    connector_coupling = pm.set_variable('1000um')
                    purcell_port_connect = pm.set_variable('500um')
                    delta_ro_purc = pm.set_variable('-200um')

                        #we could just put a end with an open but we do this so that the Purcell res is perfeclty symmetric to the readout res                            
                    with chip1(['-0.25mm',-pos_constr_ro,],[1,0]): #we go to the const ro 
                        with chip1([0,-(track_trm_ro+2*gap_trm_ro +coupling_dist)],[1,0]):#purcell const
                            constrain_pur, = elt.create_port(chip1, widths=None, name = 'constrain_pur')
                            with chip1(['0.25mm', -(pos_constr_ro -(trm_length/2+gap_trm_ro+short_trm+gap_trm_ro*2+track_trm_ro)+delta_ro_purc)],[0,-1]) : 
                                end_pur_U, = elt.draw_Xmon_end_cable(chip1,track_trm_ro, gap_trm_ro, width_U_readout, capa_ro_length,
                                                                     name = 'end_pur_U')
                            with chip1(['-0.25mm',-tune_ro],[0,-1]):
                                end_pur, =  elt.draw_end_cable(chip1,track_trm_ro, gap_trm_ro, typeEnd = 'short', name = 'end_pur' )
#                           Connector from purcell to connector 
                                with chip1([-connector_coupling, -track_trm_ro/2], [0,-1]):
                                        contact_con, = elt.create_port(chip1, widths=[track_trm_ro, track_trm_ro+2*gap_trm_ro], name = 'constrain_pur')
                                        contact_pur_mesh = chip1.rect([-track_trm_ro, -2*gap_trm_ro], [2*track_trm_ro, track_trm_ro + 2*gap_trm_ro],
                                                                      name = 'contact_pur_mesh', layer = MESH)
                                        contact_pur_mesh.assign_mesh_length(track_trm_ro/3)
                                        with chip1([purcell_port_connect, 0], [1,0]) : 
                                            port_50_purcell, = elt.draw_end_cable(chip1, track_trm_ro, gap_trm_ro, typeEnd = 'RLC', R = '50ohm', name = 'port_50_purcell')
    ##draw the cables
    #draw at the end because then otherwise compilation errors with the different level of indentations of 
    
    
  
    #draw Purcell cable
    chip1.draw_cable(end_pur_U, constrain_pur, end_pur, is_bond = is_bond, fillet= fillet, name='purcell')
    #draw cable from purcell to connector
    chip1.draw_cable(contact_con,port_50_purcell, reverse_adaptor=True, name='cable_connect_purcell')#reverse permet d'adapter la taille du port etxterieur plutot que celui du cable 
    
    
    ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
    ground_plane.subtract(chip1.entities[GAP])
    ground_plane.unite(chip1.entities[TRACK])
    ground_plane.assign_perfect_E()
    
    #chip substrate
    chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness], material='silicon', name='substrate')
    chip1.box([0,0,0],[chip_width, chip_length, 6*chip_thickness], name='vaccuum')
    
    
    path = r'D:\Users\hqc\Documents\Jules\Python\drawpy_camille_jules\K_layout'
    pm.generate_gds(path, 'purcell_seul')
    
     
    
    
    
    
    
    
def timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name


dir_path = r"D:\Users\hqc\Documents\Jules\HFSS"
project_info = ProjectInfo(r"D:\Users\hqc\Documents\Jules\\",
                               project_name  = 'Exemple_tomograpy',  # Name of the project file (string). "None" will get the current active one.
                               design_name   = 'HFSSDesign6',  # Name of the desgin file (string). "None" will get the current active one.
                               setup_name = 'Setup1')
#project_info.junctions['jtransmon'] = {'rect':'ind_JJind_rect',  'line': 'ind_JJ_junction_line',
#                      'Lj_variable':'Lj'}
project_info.ports['Readout'] = {'rect':'port_50_purcell_track',  'line': 'port_50_purcell_line',
                  'R': 50}
 
target = np.array([6500, 200])#frq and kappa
delta_ro_purc = -90 #um
connector_coupling = 540 #um


pur_inputs = []
pur_outputs = []

def cost4_camille (delta_ro_purc, connector_coupling): 
    x0 = [delta_ro_purc, connector_coupling]
    if x0 in pur_inputs : 
       freq_pur = pur_outputs[pur_inputs.index(x0)][0]
       kappa_pur = pur_outputs[pur_inputs.index(x0)][1]
    else :     
        pur_inputs.append(x0)
        x0 = np.array(x0)
        
        epr_hfss = DistributedAnalysis(project_info)
  
        opti=Optimetrics(epr_hfss.design)
        
        ###load the optimetrics module from the ansys package
        parametric_name=timestamp_name('parametric')
        print(parametric_name)
        
        ###get the list of HFSS variable
        var=epr_hfss.design.get_variables()   
        
        ###list of variable to be optimized (should probably done outside)
        name=np.array(["delta_ro_purc", "connector_coupling"])
        
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
        opti.edit_setup(parametric_name, CopyMesh=False)
        
        opti.solve_setup(parametric_name)
        
    # Reload the project with the new variation 
    # DO NOT use quantum analysis since there is no junction !
    
        epr_hfss = DistributedAnalysis(project_info)
        var_list=epr_hfss.variations
        epr_hfss.update_ansys_info()
        epr_hfss.pinfo.save()
        epr_hfss.set_variation(var_list[-1])
        freqs, Qs = epr_hfss.get_freqs_bare_pd(var_list[-1], frame=False)
        
        freq_pur = freqs[0]
        Q_pur = Qs[0]
        kappa_pur = freq_pur/Q_pur
        
        pur_outputs.append([freq_pur, kappa_pur])
    
    print('fre_pur', freq_pur)
    print('kappa_pur', kappa_pur)
    
    return(-(((freq_pur-target[0])/target[0])**2 + ((kappa_pur-target[1])/target[1])**2))
    
#    
#pbounds = {'delta_ro_purc' : (-200, 0), 'connector_coupling' :(100, 1000)} 
#bounds_transformer = SequentialDomainReductionTransformer()
#domred_optimizer = BayesianOptimization(f = cost4_camille,pbounds=pbounds, verbose=2, 
#                                          random_state=1, bounds_transformer = bounds_transformer)
#
#domred_optimizer.maximize(10, 20)