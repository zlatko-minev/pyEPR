# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 16:50:15 2020

@author: hqc
"""
#################################################
### Tuning : cat(memory) + transmon + readout + purcell ###
#################################################
import numpy as np
# Import PyEPR stuff
from pyEPR import *
from pyEPR.ansys import Optimetrics, HfssDesign
import time
from pyEPR.project_info import ProjectInfo
from pyEPR import DistributedAnalysis, QuantumAnalysis
# Import optimization stuff
import scipy.optimize as sp
from bayes_opt import BayesianOptimization
from bayes_opt import SequentialDomainReductionTransformer
# Import Drawpylib stuff
from HFSSdrawpy.utils import val
from HFSSdrawpy.core.modeler import Modeler
from HFSSdrawpy.core.body import Body
from drawpylib.parameters import TRACK, GAP, RLC, MESH
import drawpylib.cpw_elements as elt

new_geometry = True
tuning = False
# This code is to run in 5 steps, each step tunes some parts of the circuits
# 1. Transmon alone 
# 2. Memory (cat) + Transmon
# 3. Readout + transmon  
# 4. Purcell alone
# 5. Purcell + Readout

transmon = False 
memory_transmon = False
readout_transmon = False
purcell = True
purcell_readout = False


### DRAWING THE CIRCUIT###
if new_geometry :        
    pm = Modeler('hfss')
    relative = pm.set_variable('1mm')
    chip1 = Body(pm, 'chip1')
    
    is_bond = False
    track = pm.set_variable('20um')
    chip_width= pm.set_variable("9.74mm")
    chip_length=pm.set_variable("10.74mm")
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
    #        in_bot_R_port, = elt.draw_connector(chip1, pcb_track, pcb_gap, bond_length) 
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
    capa_ro_length = pm.set_variable('200um')
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
    tune_ro = pm.set_variable('1.7mm')
    pos_constr_ro = pm.set_variable('3mm') #distance from transmon to curve of readout cable
    # U capas
    width_U_readout = 2*(gap_trm_ro + track_trm_ro)+2*short_trm + track_trm_ro+2*gap_trm_ro
    width_U_mem = 2*(gap_mem + track_mem)+2*short_trm + track_trm_ro+2*gap_trm_ro
    #Purcell
    coupling_dist = pm.set_variable('10um')
    connector_coupling = pm.set_variable('1000um')
    purcell_port_connect = pm.set_variable('500um')
    delta_ro_purc = pm.set_variable('-200um')
    #Position of the memory that will change everyone's position
    pos_x_start = pm.set_variable('3mm')#Pos of the T (memory)
    pos_y_start = pm.set_variable('5mm')
    
    
    if transmon :
        pos_x_start = pm.set_variable('1mm')
        pos_y_start = pm.set_variable('1mm')
        chip_width= pm.set_variable("2mm")
        chip_length=pm.set_variable("2mm")
        with chip1([pos_x_start, pos_y_start],[-1,0]):
            # Draw transmon 
            elt.draw_half_transmon(chip1, track_trm_ro, gap_trm_ro, trm_length, '5um',Lj, name='transmon')
            
        ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
        ground_plane.subtract(chip1.entities[GAP])
        ground_plane.unite(chip1.entities[TRACK])
        ground_plane.assign_perfect_E()
            
        #chip substrate
        chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness], 
                  material='silicon', name='substrate')
        chip1.box([0,0,0],[chip_width, chip_length, 6*chip_thickness], name='vaccuum')
            
        path = r'D:\Users\hqc\Documents\Jules\Python\drawpy_camille_jules\K_layout'
        pm.generate_gds(path, 'trm_seul')
    
    if memory_transmon : 
        pos_x_start = pm.set_variable('0.5mm')
        pos_y_start = pm.set_variable('3.75mm')
        chip_width= pm.set_variable("6mm")
        chip_length=pm.set_variable("6mm")        
        with chip1([pos_x_start, pos_y_start],[0,-1]):
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
            width_U_mem = 2*(gap_mem + track_mem)+2*short_trm + track_trm_ro+2*gap_trm_ro#p-e mettre 2 short-trm differents pr la mem et le readout (2 U differents)
            with chip1([-shift_capa_mem_trm_Y, shift_capa_mem_trm_X],[1,0]):
                Xmon_mem, = elt.draw_Xmon_end_cable(chip1, track_mem, gap_mem, width_U_mem, capa_U_length, name ='U_coupler' )
            #constrain 
                with chip1([shift_capa_mem_trm_Y,-shift_const_X],[0,-1]): #[shift_const_X,-shift_capa_mem_trm_Y
                    constrain_mem, = elt.create_port(chip1, name = 'constrain_mem')    
              # Transmon
                with chip1([2*gap_mem+track_mem+trm_length/2+gap_trm_ro+short_trm,0],[0,1]):  #+c.short_trm+c.trm_length/2
                    elt.draw_half_transmon(chip1, track_trm_ro, gap_trm_ro, trm_length, '5um',Lj, name='transmon') #no port
             
        # draw mem cable
        chip1.draw_cable(T_port3, constrain_mem, Xmon_mem, is_bond = is_bond,
                                     fillet = fillet_mem, name = 'mem')
        
        ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
        ground_plane.subtract(chip1.entities[GAP])
        ground_plane.unite(chip1.entities[TRACK])
        ground_plane.assign_perfect_E()
    
        #chip substrate
        chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness], material='silicon', name='substrate')
        chip1.box([0,0,0],[chip_width, chip_length, 6*chip_thickness], name='vaccuum')
        
        path = r'D:\Users\hqc\Documents\Jules\Python\drawpy_camille_jules'
        pm.generate_gds(path, 'trm_mem')
            
    if readout_transmon :
          pos_x_start = pm.set_variable('-2mm')
          pos_y_start = pm.set_variable('2.5mm')
          chip_width= pm.set_variable("4mm")
          chip_length=pm.set_variable("6mm")         
          with chip1([pos_x_start, pos_y_start],[0,-1]):
              with chip1([-shift_capa_mem_trm_Y,-shift_capa_mem_trm_X],[-1,0]):
            # Transmon
                    with chip1([2*gap_mem+track_mem+trm_length/2+gap_trm_ro+short_trm,0],[0,1]):  #+c.short_trm+c.trm_length/2
                        elt.draw_half_transmon(chip1, track_trm_ro, gap_trm_ro, trm_length, '5um',Lj, name='transmon') #no port
             # read out capa 
                        with chip1([0,-(trm_length/2+gap_trm_ro+short_trm+gap_trm_ro*2+track_trm_ro)],[0,1]):
                            capa_ro, = elt.draw_Xmon_end_cable(chip1, track_trm_ro, gap_trm_ro,
                                                               width_U_readout, capa_ro_length,
                                                               name = 'readout_capa')
            # constrain port for readout
                        with chip1(['-0.25mm',-pos_constr_ro,],[1,0]):# '-0.25' correspond à la hauteur (par rapport au trm) du port de contrainte qui permet le couplage avec le purcell, à tuner 
                            constraint_ro, = elt.create_port(chip1, widths=None, 
                                                      name = 'constrain_ro')
            #end of cable ro
                            with chip1(['-0.25mm',tune_ro],[0,1]):
                                end_ro, = elt.draw_end_cable(chip1,track_trm_ro, gap_trm_ro, typeEnd = 'short', name = 'end_ro' )
                
                # draw ro cable
          chip1.draw_cable( capa_ro, constraint_ro, end_ro,is_bond= is_bond,
                             fillet = fillet, name ='ro')
            
          ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
          ground_plane.subtract(chip1.entities[GAP])
          ground_plane.unite(chip1.entities[TRACK])
          ground_plane.assign_perfect_E()
            
          #chip substrate
          chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness], material='silicon', name='substrate')
          chip1.box([0,0,0],[chip_width, chip_length, 6*chip_thickness], name='vaccuum')
            
          path = r'D:\Users\hqc\Documents\Jules\Python\drawpy_camille_jules\K_layout'
          pm.generate_gds(path, 'trm_readout')
     
    if purcell : 
        pos_x_start = pm.set_variable('-2.5mm')
        pos_y_start = pm.set_variable('-1mm')
        chip_width= pm.set_variable("3mm")
        chip_length=pm.set_variable("4.5mm")         
        with chip1([pos_x_start, pos_y_start],[0,-1]):
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
        
    if purcell_readout : 
         pos_x_start = pm.set_variable('-2mm')
         pos_y_start = pm.set_variable('3.75mm')
         chip_width= pm.set_variable("6mm")
         chip_length=pm.set_variable("6mm")         
         with chip1([pos_x_start, pos_y_start],[0,-1]):
            with chip1([-shift_capa_mem_trm_Y,-shift_capa_mem_trm_X],[-1,0]):
                with chip1([2*gap_mem+track_mem+trm_length/2+gap_trm_ro+short_trm,0],[0,1]):  #+c.short_trm+c.trm_length/2
             # read out capa 
                    with chip1([0,-(trm_length/2+gap_trm_ro+short_trm+gap_trm_ro*2+track_trm_ro)],[0,1]):
                        capa_ro, = elt.draw_Xmon_end_cable(chip1, track_trm_ro, gap_trm_ro,
                                width_U_readout, capa_ro_length,
                                name = 'readout_capa')
            # constrain port for readout
                    with chip1(['-0.25mm',-pos_constr_ro,],[1,0]):# '-0.25' correspond à la hauteur (par rapport au trm) du port de contrainte qui permet le couplage avec le purcell, à tuner 
                        constraint_ro, = elt.create_port(chip1, widths=None, 
                                                      name = 'constrain_ro')
            #end of cable ro
                        with chip1(['-0.25mm',tune_ro],[0,1]):
                            end_ro, = elt.draw_end_cable(chip1,track_trm_ro, gap_trm_ro, typeEnd = 'short', name = 'end_ro' )
          
                    
                   #purcell filter
                        coupling_dist = pm.set_variable('10um')
    #                    tune_pur = pm.set_variable('6.3mm')
                        connector_coupling = pm.set_variable('1000um')
                        purcell_port_connect = pm.set_variable('500um')
                        delta_ro_purc = pm.set_variable('-200um')
    #                with chip1([0,-tune_pur],[0,-1]):
    #                    end_pur_U, = elt.draw_Xmon_end_cable(chip1,track_trm_ro, gap_trm_ro, width_U_readout, capa_ro_length,
    #                            name = 'end_pur_U')
                        #we could just put a end with an open but we do this so that the Purcell res is perfeclty symmetric to the readout res                            
                    with chip1(['-0.25mm',-pos_constr_ro,],[1,0]): #we go to the const ro 
                        with chip1([0,-(track_trm_ro+2*gap_trm_ro +coupling_dist)],[1,0]):#purcell const
                                constrain_pur, = elt.create_port(chip1, widths=None, name = 'constrain_pur')
                                with chip1(['0.25mm', -(pos_constr_ro -(trm_length/2+gap_trm_ro+short_trm+gap_trm_ro*2+track_trm_ro)+delta_ro_purc)],[0,-1]) : 
                                    end_pur_U, = elt.draw_Xmon_end_cable(chip1,track_trm_ro, gap_trm_ro, width_U_readout, capa_ro_length,
                                                                         name = 'end_pur_U')
                                with chip1(['-0.25mm',-tune_ro],[0,-1]):
                                    end_pur, =  elt.draw_end_cable(chip1,track_trm_ro, gap_trm_ro, typeEnd = 'short', name = 'end_pur' )
                               #Connector from purcell to connector 
                                    with chip1([-connector_coupling, -track_trm_ro/2], [0,-1]):
                                            contact_con, = elt.create_port(chip1, widths=[track_trm_ro, track_trm_ro+2*gap_trm_ro], name = 'constrain_pur')
                                            contact_pur_mesh = chip1.rect([-track_trm_ro, -2*gap_trm_ro], [2*track_trm_ro, track_trm_ro + 2*gap_trm_ro],
                                                                          name = 'contact_pur_mesh', layer = MESH)
                                            contact_pur_mesh.assign_mesh_length(track_trm_ro/3)
                                            with chip1([purcell_port_connect, 0], [1,0]) :            
                                                port_50_purcell, = elt.draw_end_cable(chip1, track_trm_ro, gap_trm_ro, typeEnd = 'RLC', R='50ohm', name = 'port_50_purcell')
  
        # draw ro cable
            chip1.draw_cable( capa_ro, constraint_ro, end_ro,is_bond= is_bond,
                         fillet = fillet, name ='ro')
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
            pm.generate_gds(path, 'purcell_readout')

        
    #######################
#########    TUNING     ##########
    ########################
    
if tuning :    
    def timestamp_name(name):
        secondsSinceEpoch = time.time()
        timeObj = time.localtime(secondsSinceEpoch)
        timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,
                                       timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
        return timestamp+'_'+name
    
    
    dir_path = r"D:\Users\hqc\Documents\Jules\HFSS"
    project_info = ProjectInfo(r"D:\Users\hqc\Documents\Jules\\",
                               project_name  = 'Exemple_tomograpy',  # Name of the project file 
                               design_name   = 'HFSSDesign',  # Name of the desgin file 
                               setup_name = 'Setup1')
    
    
    if transmon : 
        target = np.array([5500, 200])  #freq_tr / ker_tr 
        project_info.junctions['jtransmon'] = {'rect':'ind_JJind_rect',
                          'line': 'ind_JJ_junction_line', 'Lj_variable':'Lj'}
    if memory_transmon : 
        target = np.array([5500, 8000,200, 0.3])  #freq_tr / freq_mem / ker_tr / chi_tr_mem
        project_info.junctions['jtransmon'] = {'rect':'ind_JJind_rect',
                          'line': 'ind_JJ_junction_line', 'Lj_variable':'Lj'}
    if readout_transmon : 
        target = np.array([5500, 6500,200, 3])  #freq_tr / freq_ro / ker_tr / chi_tr_ro
        project_info.junctions['jtransmon'] = {'rect':'ind_JJind_rect',
                          'line': 'ind_JJ_junction_line', 'Lj_variable':'Lj'}
    if purcell : 
        target = np.array([6500,200])  # freq_ro / kappa_purc
        project_info.ports['Readout'] = {'rect':'port_50_purcell_track',  'line': 'port_50_purcell_line',
                      'R': 50}
    if purcell_readout : 
        target = np.array([6500, 6500,3, 200])  #freq_ro / freq_purc / kappa_ro / kappa_pur
        project_info.ports['Readout'] = {'rect':'port_50_purcell_track',  'line': 'port_50_purcell_line',
                      'R': 50} 
    
    def cost_camille (Lj, trm_length, shift_capa_mem_trm_Y, capa_U_length, tune_ro, capa_ro_length,
                      delta_ro_purc, connector_coupling, coupling_dist):
        if transmon : 
            x0 = np.array([Lj, trm_length])
        if memory_transmon : 
            x0 = np.array([Lj, trm_length, shift_capa_mem_trm_Y, capa_U_length])
        if readout_transmon : 
            x0 = np.array([Lj, trm_length, tune_ro, capa_ro_length])
        if purcell : 
            x0 = np.array([delta_ro_purc, connector_coupling])
        if purcell_readout : 
            x0 = np.array([coupling_dist, connector_coupling, tune_ro, delta_ro_purc])
        print('x0 = ', x0)
        
        epr_hfss = DistributedAnalysis(project_info)
        
        ###OPTIMETRICS
        opti=Optimetrics( epr_hfss.design)
                parametric_name=timestamp_name('parametric')
        print(parametric_name)
        
        ###HFSS variables
        var=epr_hfss.design.get_variables()   
        
        ###list of variable to be optimized (should probably done outside)
        if transmon : 
            name=np.array(["trm_length","Lj"])
        if memory_transmon : 
            name=np.array(["trm_length","Lj","shift_capa_mem_trm_Y", "capa_U_length"])
        if readout_transmon : 
            name=np.array(["trm_length","Lj","tune_ro", "capa_ro_length"])
        if purcell : 
            name=np.array(["delta_ro_purc", "connector_coupling"])
        if purcell_readout : 
            name=np.array(["coupling_dist","connector_coupling","tune_ro", "delta_ro_purc"])
        
        ###UNITS 
        units=[var[key][-2:] for key in name]
        var_list=np.ones(x0.shape[0]).astype(str)
        for i in range(x0.shape[0]):
                var_list[i]=str(x0[i])+units[i]
        print(var_list)
                
        ###create and save the array of string with the correct format to import to HFSS
        arr=np.vstack([name,var_list])
        np.savetxt(dir_path+"\%s.txt"%parametric_name, arr, fmt='%s', delimiter='; ',
                   newline='\n', header='', footer='', comments='# ', encoding=None)
        
        ###import the parametric setup with the list of variations (function I added to the ansys package)
        opti.import_setup(parametric_name,dir_path+"\%s.txt"%parametric_name)
        
        opti.solve_setup(parametric_name)
            
        # Reload the project with the new variation 
        epr_hfss = DistributedAnalysis(project_info)
        
        if transmon or memory_transmon or readout_transmon : 
            var_list=list(epr_hfss.variations)
            epr_hfss.do_EPR_analysis([var_list[-1]])
            epr = QuantumAnalysis(epr_hfss.data_filename,[var_list[-1]])
            epr.analyze_all_variations([var_list[-1]],cos_trunc = 5, fock_trunc = 4)
            
            freqs=np.array(epr.get_frequencies()).T   
            chi_dic = epr.results.get_chi_O1()
            chis = np.abs(np.array(chi_dic[var_list[-1]]))
            freq_dic = epr.results.get_frequencies_O1()
            freq = np.abs(np.array(freq_dic[var_list[-1]]))
            chi = np.diag(chis)   
        
        if purcell : 
            var_list=epr_hfss.variations
            epr_hfss.update_ansys_info()
            epr_hfss.pinfo.save()
            epr_hfss.set_variation(var_list[-1])
            freqs, Qs = epr_hfss.get_freqs_bare_pd(var_list[-1], frame=False)
        
        if purcell_readout : 
            var_list=epr_hfss.variations
            epr_hfss.update_ansys_info()
            epr_hfss.pinfo.save()
            epr_hfss.set_variation(var_list[-1])
            freqs, Qs = epr_hfss.get_freqs_bare_pd(var_list[-1], frame=False)
            Q = [Qs[0], Qs[1]]
            ### define the purcell as the mode with the largest losses = smallest Q
            index={}
            index['pur']=np.argsort(Q)[-2]
            index['ro']=np.argsort(Q)[-1]
                
    
            freq_pur = freqs[index['pur']]*1000
            freq_ro = freqs[index['ro']]*1000
            kappa_pur = freq_pur/Qs[index['pur']]
            kappa_ro = freq_ro/Qs[index['ro']]
        
        # GIVE THE FREQUENCIES
        if transmon : 
            freq_trm = freq[0]
            chi_trm = chi[0]
            print('freq_trm =', freq_trm)
            print('chi_trm =', chi_trm)
            return(-(((freq_tr-target[0])/target[0])**2 + ((ker_ter-target[1])/target[1])**2))
    
        if memory_transmon : 
            ### define the transmon as the mode with the largest anharmanocity###
            index={}
            index['trm']=np.argsort(chi)[-1]
            index['mem']=np.argsort(chi)[-2]
            
            chi_mem_trm = chis[index['mem']][index['trm']]
            chi_trm = chi[index['trm']]
            freq_trm = freq[index['trm']]
            freq_mem = freq[index['mem']]
        
            print('freq_trm =', freq_trm)
            print('freq_mem =', freq_mem)
            print('chi_trm =', chi_trm) 
            print('chi_mem_trm =', chi_mem_trm)
            return (-(((freq_trm-target[0])/target[0])**2 + ((freq_mem-target[1])/target[1])**2
                      +((chi_trm-target[2])/target[2])**2 + ((chi_mem_trm-target[3])/target[3])**2))
            
        if readout_transmon : 
            ### define the transmon as the mode with the largest anharmanocity###
            index={}
            index['trm']=np.argsort(chi)[-1]
            index['ro']=np.argsort(chi)[-2]
            
            chi_ro_trm = chis[index['ro']][index['trm']]
            chi_trm = chi[index['trm']]
            freq_trm = freq[index['trm']]
            freq_ro = freq[index['ro']]
        
            print('freq_ro =', freq_ro)
            print('freq_trm =', freq_trm)
            print('chi_trm =', chi_trm)
            print('chi_ro_trm =', chi_ro_trm)
            return (-(((freq_trm-target[0])/target[0])**2 + ((freq_ro-target[1])/target[1])**2
                  +((chi_trm-target[2])/target[2])**2 + ((chi_ro_trm-target[3])/target[3])**2))   
            
        if purcell : 
            total_Q_from_HFSS = np.array(epr.Qs)[:,var]
            
            freq_pur = freq[0]
            kappa_pur = freq_pur/total_Q_from_HFSS[0]
            
            print('freq_pur =', freq_pur)
            print('kappa_pur =', kappa_pur)
            return(-(((freq_pur-target[0])/target[0])**2 + ((kappa_pur-target[1])/target[1])**2))
    
        if purcell_readout : 
            total_Q_from_HFSS = np.array(epr.Qs)[:,var]
            
            ### define the purcell as the mode with the largest losses (lowest Q)
            index={}
            index['pur']=np.argsort(total_Q_from_HFSS)[-2]
            index['ro']=np.argsort(total_Q_from_HFSS)[-1]
            
            freq_pur = freq[index['pur']]
            freq_ro = freq[index['ro']]
            kappa_pur = freq_pur/total_Q_from_HFSS[index['pur']]
            kappa_ro = freq_ro/total_Q_from_HFSS[index['ro']]

            print('freq_pur =', freq_pur)
            print('freq_ro =', freq_ro)
            print('kappa_pur =', kappa_pur)
            print('kappa_ro =', kappa_ro)
            return (-(((freq_ro-target[0])/target[0])**2 + ((freq_pur-target[1])/target[1])**2
                      +((kappa_ro-target[2])/target[2])**2 + ((kappa_pur-target[3])/target[3])**2))        
        
       
        Lj, trm_length, shift_capa_mem_trm_Y, capa_U_length, tune_ro, capa_ro_length,
                      delta_ro_purc, connector_coupling, coupling_dist
        
        
##### BAYESIAN OPTIMIZATION ###
######The bounds are given by the previous calculations

                      
    if transmon : 
        pbounds = {'Lj': (8.4, 8.55), 'trm_length': (0.4, 0.5), 'shift_capa_mem_trm_Y': (-2328, -2328),
                   'capa_U_length' : (155, 165), "tune_ro" : (1, 1.4), "capa_ro_length" :(110, 200),
               "delta_ro_purc" : (-200, 30), "connector_coupling" : (100, 1000), "coupling_dist" :(7, 15)} 
#not finished yet...  
    bounds_transformer = SequentialDomainReductionTransformer()
    domred_optimizer = BayesianOptimization(f = cost2_camille,pbounds=pbounds, verbose=2, 
                                              random_state=1, bounds_transformer = bounds_transformer)
#    
    #domred_optimizer.maximize(10, 20)    