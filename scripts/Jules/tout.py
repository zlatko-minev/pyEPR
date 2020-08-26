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
    pos_constr_ro = pm.set_variable('3mm') #distance from transmon to curve of readout cable
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
    pos_y_mem = pm.set_variable('3.25mm')
    
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
                                            contact_pur_mesh.assign_mesh_length (track_trm_ro/3)
                                            with chip1([purcell_port_connect, 0], [1,0]) : 
                                                    port_50_purcell, = elt.draw_end_cable(chip1, track_trm_ro, gap_trm_ro, typeEnd = 'RLC', R='50ohm', name = 'port_50_purcell')
    ##draw the cables
    #draw at the end because then otherwise compilation errors with the different level of indentations of 
    
    # draw mem cable
    chip1.draw_cable(T_port3, constrain_mem, Xmon_mem, is_bond = is_bond,
                                 fillet = fillet_mem, name = 'mem')
    
    
    # draw ro cable
    chip1.draw_cable( capa_ro, constraint_ro, end_ro,is_bond= is_bond,
                     fillet = fillet, name ='ro')
    #draw Purcell cable
    chip1.draw_cable(end_pur_U, constrain_pur, end_pur, is_bond = is_bond, fillet= fillet, name='purcell')
    #draw cable from purcell to connector
    chip1.draw_cable(contact_con,port_50_purcell, reverse_adaptor=True, name='cable_connect_purcell')
    #reverse permet d'adapter la taille du port etxterieur plutot que celui du cable 
    
    # Meshing 
    ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
    ground_plane.subtract(chip1.entities[GAP])
    ground_plane.unite(chip1.entities[TRACK])
    ground_plane.assign_perfect_E()
    
    #chip substrate
    chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness], material='silicon', name='substrate')
    chip1.box([0,0,0],[chip_width, chip_length, 6*chip_thickness], name='vaccuum')
    
    
    path = r'D:\Users\hqc\Documents\Jules\Python\drawpy_camille_jules\K_layout'
    pm.generate_gds(path, 'tout')
    
    
    
    
    
#TUNING
def timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name


dir_path = r"D:\Users\hqc\Documents\Jules\HFSS"
project_info = ProjectInfo(r"D:\Users\hqc\Documents\Jules\\",
                               project_name  = 'Exemple_tomograpy',  # Name of the project file (string). "None" will get the current active one.
                               design_name   = 'HFSSDesign10',  # Name of the desgin file (string). "None" will get the current active one.
                               setup_name = 'Setup1')
project_info.junctions['jtransmon'] = {'rect':'ind_JJind_rect',  'line': 'ind_JJ_junction_line',
                      'Lj_variable':'Lj'}
project_info.ports['Readout'] = {'rect':'port_50_purcell_track',  'line': 'port_50_purcell_line',
                  'R': 50}


Lj = 8.30073 
trm_length = 0.4653858092535451 #mm
shift_capa_mem_trm_Y = -2544.027 #um
capa_U_length = 158.7304630287759 #um

coupling_dist = 13 # um 
connector_coupling = 540 #um
tune_ro = 1.050284085297964 #mm
delta_ro_purc = -85 # um

capa_ro_length = 136.73585031623904 # um



x0 = np.array([trm_length, Lj, shift_capa_mem_trm_Y, capa_U_length, coupling_dist,
               connector_coupling, tune_ro, delta_ro_purc, capa_ro_length])
epr_hfss = DistributedAnalysis(project_info)
###load the optimetrics module from the ansys package
opti=Optimetrics(epr_hfss.design)

###load the optimetrics module from the ansys package
parametric_name=timestamp_name('parametric')
print(parametric_name)

###get the list of HFSS variable
var=epr_hfss.design.get_variables()   

###list of variable to be optimized (should probably done outside)
name=np.array(["trm_length","Lj","shift_capa_mem_trm_Y", "capa_U_length","coupling_dist",
               "connector_coupling", "tune_ro", "delta_ro_purc", "capa_ro_length"])

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

epr = QuantumAnalysis(epr_hfss.data_filename,[var_list[-1]])
epr.analyze_all_variations([var_list[-1]],cos_trunc = 5, fock_trunc = 4)

chi_dic = epr.results.get_chi_O1()
chis = np.abs(np.array(chi_dic[var_list[-1]]))
freq_dic = epr.results.get_frequencies_O1()
freq = list(np.abs(np.array(freq_dic[var_list[-1]])))
chi = np.diag(chis)   
Qs = epr_hfss.get_freqs_bare_pd(var_list[-1], frame=False)[1]
Q = [Qs[0], Qs[1], Qs[2], Qs[3]]


 ### sorting modes  
### define the purcell as the mode with the largest losses = smallest Q
index={}
index['pur']=np.argsort(Q)[0]
index['ro']=np.argsort(Q)[1]
 ### define the qubit as the mode with the largest anharmanocity
index['trm']=np.argsort(chi)[-1]


freq_trm = freq[index['trm']]
freq_ro = freq[index['ro']]
freq_pur = freq[index['pur']]

# The memory is the last frequency of the list
known_freqs = [freq_trm, freq_ro, freq_pur]
index['mem'] = freq.index(list(set(freq)-set(known_freqs)))

freq_mem = freq[index['mem']]



chi_trm = chi[index['trm']]
chi_ro_trm = chis[index['ro']][index['trm']]
chi_mem_trm = chis[index['mem']][index['trm']]


kappa_pur = freq_pur/Qs[index['pur']]
kappa_ro = freq_ro/Qs[index['ro']]


print('freq_mem =', freq_mem)
print('freq_trm =', freq_trm)
print('freq_ro =', freq_ro)
print('freq_pur =', freq_pur)
print('chi_trm =', chi_trm)
print('chi_mem_trm =', chi_mem_trm)
print('chi_ro_trm =', chi_ro_trm)
print('kappa_pur =', kappa_pur)
print('kappa_ro =', kappa_ro)










        

