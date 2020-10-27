# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 11:26:36 2020

@author: berdou
"""

import numpy as np
from HFSSdrawpy.utils import val
from HFSSdrawpy.core.modeler import Modeler
from HFSSdrawpy.core.body import Body
from drawpylib.parameters import TRACK, GAP, RLC, MESH, DEFAULT
import drawpylib.cpw_elements as elt



pm = Modeler('hfss')
relative = pm.set_variable('1mm')
chip1 = Body(pm, 'chip1')

chip_width= pm.set_variable("5.25mm")
chip_length=pm.set_variable("7.25mm")

chip_thickness= pm.set_variable("280um")    
pcb_thickness= pm.set_variable("320um")
chip_height= pm.set_variable(6*chip_thickness)
is_bond = True

track= pm.set_variable("25um")
gap = pm.set_variable("15um")
fillet =pm.set_variable('200um')

# Capa buff mem
track_capa = pm.set_variable('25um')
gap_capa = pm.set_variable('15um')
capa_length = pm.set_variable('1000um')
#capa U mem trm
shift_capa_mem_trm_Y = pm.set_variable('2500um')
shift_capa_mem_trm_X = pm.set_variable('3500um')
fillet_mem = pm.set_variable('200um')
track_mem = pm.set_variable('25um')
gap_mem = pm.set_variable('15um')
capa_ro_length = pm.set_variable('136.7um')
capa_mem_length = pm.set_variable('158.7um')
short_trm = pm.set_variable('5um')
shift_const_X = pm.set_variable('800um')
shift_const_Y = pm.set_variable('1900um')
#Transmon
track_trm_ro = pm.set_variable('25um')
gap_trm_ro = pm.set_variable('15um')
Lj = pm.set_variable('8.5nH')
trm_length = pm.set_variable('0.465mm')
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
pos_x_mem = pm.set_variable('0.4mm')
pos_y_mem = pm.set_variable('3.1mm')

res_line_length = pm.set_variable('30um')

with chip1([pos_x_mem, pos_y_mem],[0,-1]):
          T_portL, T_portR, T_port3,= elt.draw_T(chip1, track_capa, gap_capa)
          T_mesh = chip1.rect([-track_capa, - gap_capa],
            [2*track_capa, track_capa + 2*gap_capa], name = 'T_mesh', layer = MESH)
          T_mesh.assign_mesh_length(track_capa/3)
          with chip1([capa_length/2,0],[1,0]):
              capa_end_L, = elt.draw_end_cable(chip1, track_capa,gap_capa
                                               ,typeEnd='open', name = 'capa_end_L')
          with chip1([-capa_length/2,0],[-1,0]):
              capa_end_R, = elt.draw_end_cable(chip1, track_capa, gap_capa,
                                               typeEnd = 'open', name = 'capa_end_R')
          chip1. draw_cable(T_portL, capa_end_L,
                       is_bond=is_bond, fillet=fillet, name ='capa_L')
          chip1.draw_cable(T_portR, capa_end_R, 
                           is_bond= is_bond, fillet= fillet, name = 'capa_R')
          
          #capa U mem trm
          with chip1([shift_capa_mem_trm_Y,shift_capa_mem_trm_X],[-1,0]):
              Xmon_mem, = elt.draw_Xmon_end_cable(chip1, track_mem, gap_mem,
                                width_U_mem, capa_mem_length, name ='U_coupler' )
              # memory line
              with chip1([-track ,0],[-1,0]):
                  points =[('0um', '0um'), (res_line_length, '0um')]
                  mem_line = chip1.polyline(points, closed=False,
                                            name='mem_line', layer = DEFAULT)
                                             
          #constrain 
              with chip1([shift_capa_mem_trm_Y,shift_const_X],[0,1]): #[shift_const_X,shift_capa_mem_trm_Y
                  constrain_mem, = elt.create_port(chip1, name = 'constrain_mem')    
            # Transmon
              with chip1([2*gap_mem+track_mem+trm_length/2+gap_trm_ro+short_trm,0],[0,1]):  #+c.short_trm+c.trm_length/2
                  elt.draw_half_transmon(chip1, track_trm_ro, gap_trm_ro,
                                         trm_length, '5um',Lj, name='transmon') #no port
           # read out capa 
                  with chip1([0,-(trm_length/2+gap_trm_ro+short_trm+gap_trm_ro*2+track_trm_ro)],[0,1]):
                      capa_ro, = elt.draw_Xmon_end_cable(chip1, track_trm_ro, gap_trm_ro,
                              width_U_readout, capa_ro_length, name = 'readout_capa')
                 # Readout line
                      with chip1([-track ,0],[-1,0]):
                          points =[('0um','0um'),(res_line_length,'0um')]
                          ro_line = chip1.polyline(points, closed=False,
                                                   name='ro_line', layer = DEFAULT)
                          # rect_test = chip1.rect(['10um', '20um'], ['50um', '50um'], name = 'contact_pur_mesh', layer = MESH)
                
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

                                 # Purcel line line
                                  with chip1([-track ,0],[-1,0]):
                                      points =[('0um','0um'),(res_line_length,'0um')]
                                      pur_line = chip1.polyline(points, closed=False,
                                                   name='pur_line', layer = DEFAULT)
                          # rect_test = chip1.rect(['10um', '20um'], ['50um', '50um'], name = 'contact_pur_mesh', layer = MESH)
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
chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness],
          material='silicon', name='substrate')
chip1.box([0,0,0],[chip_width, chip_length, chip_height], name='vaccuum')


path = r'D:/Jules/Autotune/gds_files'
pm.generate_gds(path, 'autotun_test')
  
 

 
