import os
from HFSSdrawpy import Modeler, Body
import drawpylib.cpw_elements as elt
from HFSSdrawpy.utils import val
from quantrolib.parameters import TRACK, GAP, RLC, MESH, MASK, DEFAULT, ELEC, \
    eps
import numpy as np
pm = Modeler('hfss')
pm.is_hfss=True
#############################################################################################
############################## Draw the box #################################################
#############################################################################################

copper_box = Body(pm, name='copper_box')

# create the copper box, values given by Eric
box_width = pm.set_variable('6mm')
box_length = pm.set_variable('34mm')
box_height = pm.set_variable('20mm')

box0=copper_box.box([-0.5*box_width,-0.5*box_length,-0.5*box_height], [box_width, box_length, box_height], material='vacuum')

# create the ports
################ 2 Pins parameters and drawing
connect_inner_radius = pm.set_variable('0.25mm')
connect_outer_radius = pm.set_variable('0.6mm')
connect_length = pm.set_variable('2mm')
connect_penetrationlength1 = pm.set_variable('3mm')

### Pin position 1
connect_positionX1 = -0.5*box_width-connect_length
connect_positionY1 = pm.set_variable('5mm')
connect_positionZ1 = pm.set_variable('5mm')

cylinder_outer1=copper_box.cylinder([connect_positionX1, connect_positionY1, connect_positionZ1],connect_outer_radius,connect_length,"X",name='connect_outer1')
cylinder_outer1.assign_material('vacuum')
copper_box.unite([box0, cylinder_outer1])
box0.assign_perfect_E()

#pin 1
cylinder_inner1=copper_box.cylinder([connect_positionX1, connect_positionY1, connect_positionZ1], connect_inner_radius, connect_length+connect_penetrationlength1,"X",name='connect_inner1')
cylinder_inner1.assign_material('vacuum')
cylinder_inner1.assign_perfect_E()

#### Add ports add the end of the pins
connect_outer_port1=copper_box.disk([connect_positionX1, connect_positionY1, connect_positionZ1], connect_outer_radius,"X",name='connect_outer_port1')
connect_inner_port1=copper_box.disk([connect_positionX1, connect_positionY1, connect_positionZ1], connect_inner_radius,"X",name='connect_inner_port1')
connect_outer_port1.subtract([connect_inner_port1])

start1 = [connect_positionX1, connect_positionY1+connect_inner_radius, connect_positionZ1]
end1 = [connect_positionX1, connect_positionY1+connect_outer_radius, connect_positionZ1]

connect_outer_port1.assign_lumped_RLC([end1,start1],["50ohm", 0, 0])

#####################################################################################
###################### Draw the chip ###################################################
#####################################################################################

substrate = Body(pm, name='substrate')

# Parameters for silicon chip: fixed parameters
silicon_width = pm.set_variable('3mm')
silicon_length = pm.set_variable('8mm')
silicon_thickness = pm.set_variable('300um')

silicon_chip=substrate.box([-silicon_width/2, -silicon_length/2, 0], [silicon_width, silicon_length, -silicon_thickness], material='silicon_mk')



######################################################################################
######################### RESONATOR DESIGN ##########################################
#######################################################################################

reso = Body(pm, name='reso')



Jinduc = pm.set_variable( '10nH')
pad_length=pm.set_variable('0.095mm')
pad_width=pm.set_variable('0.7mm')
pad_spacing = pm.set_variable('0.07mm')


cutout_size = [3*pad_length,2*pad_width]
pad_size = [pad_length, pad_width]
Jwidth = '10um'
Jlength = '0.5um'
gap='12um'
track='21um'
with reso([0, 0], [1, 0]):
    elt.draw_IBM_transmon(reso, cutout_size=cutout_size,
                          pad_spacing=pad_spacing,
                          pad_size=pad_size,
                          Jwidth=Jwidth, Jlength=Jlength,
                          track=track,
                          gap=gap,
                          Jinduc=Jinduc,
                          nport=0, fillet='0.015mm', name='qubit')    

for elmt in reso.entities[TRACK]:
    elmt.assign_perfect_E()


name_list=np.array(["connect_penetrationlength1","pad_length","pad_width","Jinduc","box_height","pad_spacing"])
unit_list=np.array(["mm","mm","mm","nH","mm","mm"])
x0=np.array([ 0.87991289,  0.35229587,  0.35229587,  9.2995847 , 24.99623459,  0.15      ])
for value,unit,name in zip(x0,unit_list,name_list):
    pm.set_variable(str(value)+unit, name)