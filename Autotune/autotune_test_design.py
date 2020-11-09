"""
Created on Tue Sep  1 11:19:53 2020

@author: jules

"""

from HFSSdrawpy.core.modeler import Modeler
from HFSSdrawpy.core.body import Body
from drawpylib.parameters import TRACK, GAP, RLC, MESH, DEFAULT
import drawpylib.cpw_elements as elt



pm = Modeler('gds')
relative = pm.set_variable('1mm')
chip1 = Body(pm, 'chip1')

chip_width= pm.set_variable("3500um")
chip_length=pm.set_variable("3500um")
chip_thickness= pm.set_variable("280um")
pcb_thickness= pm.set_variable("320um")
vaccuum_thickness= pm.set_variable('700um')
fillet = pm.set_variable('0mm')

track = pm.set_variable('25um')
gap = pm.set_variable('15um')

delta = pm.set_variable('250um')
cpl_ro_trm = pm.set_variable('150um')
Xmon_gnd = pm.set_variable('5um')
Xmon_width = 2*Xmon_gnd + 4*gap + 3*track
width_J = pm.set_variable('5um')
gnd_cpl_50ohm = pm.set_variable('5um')

lines_lenght = pm.set_variable('30um')
port_len = pm.set_variable('450um')


# Tunable parameters in the example
tune_mem = pm.set_variable('1000um')
tune_ro = pm.set_variable('2000um')
Lj = pm.set_variable('10nH')
cpl_ro_50ohm = pm.set_variable('200um')
trm_length = pm.set_variable('500um')
cpl_mem_trm = pm.set_variable('200um')






#Memory
with chip1([chip_width/10, tune_mem], [0,-1]):
    first_port1, = elt.draw_end_cable(chip1, track, gap, typeEnd='open',
                              name = 'first_port1')
with chip1([2*chip_width/10, chip_length-delta], [-1,0]):
    first_cst_port1, = elt.create_port(chip1, widths=None, name = 'first_cst_port1')

with chip1([4*chip_width/10, delta], [-1,0]):
    second_cst_port1, = elt.create_port(chip1, widths=None, name = 'second_cst_port1')

with chip1([5*chip_width/10, chip_length/3], [0,1]):
    Xmon_end_mem, = elt.draw_Xmon_end_cable(chip1, track, gap, Xmon_width, cpl_mem_trm ,
                                      name = 'Xmon_end_mem')
    # memory line
    with chip1([-track ,0],[-1,0]):
        points =[('0um', '0um'), (lines_lenght, '0um')]
        mem_line = chip1.polyline(points, closed=False,
                                  name='mem_line', layer = DEFAULT)
    # Transmon
    with chip1([3*gap + track + Xmon_gnd , 0], [-1,0]):
        port_trm1, = elt.draw_end_cable(chip1, track, gap, typeEnd='open',
                                    name = 'port_trm1')
        with chip1([-trm_length/2, track/2 + gap/2], [0,1]):
            ind1 = elt.draw_ind_inline(chip1, track, gap, gap, width_J, 2*Lj,
                            name= 'ind_JJ1')#, premesh=True):
        with chip1([-trm_length/2, -track/2 - gap/2], [0,1]):
            ind2 = elt.draw_ind_inline(chip1, track, gap, gap, width_J,2*Lj,
                            name= 'ind_JJ2')#, premesh=True):
        with chip1([-trm_length, 0], [-1,0]):
            port_trm2, = elt.draw_end_cable(chip1, track, gap, typeEnd='open',
                                            name = 'port_trm2')
    
        with chip1([-(trm_length + 3*gap + track + Xmon_gnd), 0], [1,0]):
            Xmon_end_ro, = elt.draw_Xmon_end_cable(chip1, track, gap, Xmon_width,
                                                    cpl_ro_trm ,name = 'Xmon_end_ro')
            # readout line
            with chip1([-track ,0],[-1,0]):
                points =[('0um', '0um'), (lines_lenght, '0um')]
                ro_line = chip1.polyline(points, closed=False,
                                      name='ro_line', layer = DEFAULT)
# Readout           
with chip1([6*chip_width/10, chip_length-delta], [1,0]):
    second_cst_port2, = elt.create_port(chip1, widths=None,
                                        name = 'second_cst_port2')
with chip1([8*chip_width/10, delta], [1,0]):
    first_cst_port2, = elt.create_port(chip1, widths=None,
                                       name = 'first_cst_port2')
with chip1([9*chip_width/10, tune_ro], [0,1]):
    capa_cpl, = elt.draw_Xmon_end_cable(chip1, track, gap, Xmon_width, cpl_ro_50ohm ,
                                      name = 'capa_cpl')
        
    # 50 ohm port
    with chip1([3*gap + track + gnd_cpl_50ohm, 0], [-1,0]):
        first_50port, = elt.draw_end_cable(chip1, track, gap, typeEnd='open',
                             name = 'first_50port')
    with chip1([port_len, 0], [1,0]):
        Rport, = elt.draw_end_cable(chip1, track, gap, typeEnd='RLC', R='50ohm',
                             name = 'Rport')
 
    
chip1.draw_cable(first_port1, first_cst_port1, second_cst_port1, Xmon_end_mem,
                  name = 'mem', is_bond = False)
chip1.draw_cable(port_trm1, port_trm2,  name = 'trm', is_bond = False)
chip1.draw_cable(capa_cpl, first_cst_port2, second_cst_port2, Xmon_end_ro,
                  name = 'ro', is_bond = False)
chip1.draw_cable(first_50port, Rport, name = 'res_port', is_bond = False)



ground_plane = chip1.rect([0, 0], [chip_width, chip_length], layer=TRACK)
ground_plane.subtract(chip1.entities[GAP])
ground_plane.unite(chip1.entities[TRACK])
ground_plane.assign_perfect_E()

#chip substrate
chip1.box([0,0,-chip_thickness],[chip_width, chip_length, chip_thickness],
          material='silicon', name='substrate')
chip1.box([0,0,-chip_thickness],[chip_width, chip_length, -pcb_thickness],
          material='Rogers TMM 10i (tm)', name='pcb')
chip1.box([0,0,0],[chip_width, chip_length, vaccuum_thickness], name='vaccuum')


path = r'D:/Jules/Spice project/klayout'
pm.generate_gds(path, 'mem_trm_ro_50')


