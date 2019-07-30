"""
Created on Fri Oct 30 14:21:45 2015

Configuration file 
"""
from . import AttrDict

#------------------------------------------------------------
# Directories 
'Folder to save result data in'
root_dir = r'D:\zkm-data-pyEPR'

'''
pyEPR_default_options
	Default options to use

method_calc_P_mj  :  
	'line_voltage' : Uses the line voltage integral
	'J_surf_mag'   : takes the avg. Jsurf over the rect. Make sure you have seeded 
					 lots of tets here. I recommend starting with 4 across smallest dimension.
					 Multi-junction calculation of energy participation ratio matrix based on <I_J>. Current is integrated average of J_surf by default: (zkm 3/29/16)
            		 Will calculate the Pj matrix for the selected modes for the given junctions junc_rect array & length of juuncs
'''
options_hfss = AttrDict(dict(
	method_calc_P_mj = 'line_voltage', # 'line_voltage' or 'J_surf_mag'
	save_mesh_stats  = True,
))



class Dissipation_params:
	''' Loss properties of various materials and surfaces '''
	#TOOD: Turn into a dictionary 

	# bulk dielectric:
	# refs:  https://arxiv.org/abs/1308.1743
	#        http://arxiv.org/pdf/1509.01854.pdf
	tan_delta_sapp = 1e-6  # tan(delta) for bulk surface
	epsi           = 10    # dielectric

	# surface dielectric:
	# ref: http://arxiv.org/pdf/1509.01854.pdf
	th             = 3e-9  # surface dielectric (dirt) thickness
	eps_r          = 10    # surface dielectric (dirt) constant
	tan_delta_surf = 1e-3  # surface dielectric (dirt) loss tangent, tan(delta)

	# thin-film surface loss:
	# ref:  https://arxiv.org/abs/1308.1743
	surface_Rs     = 250e-9# Ohms

	# seams current loss:
	# ref: http://arxiv.org/pdf/1509.01119.pdf
	gseam          = 1.0e3 # per Ohm meter: seam conductance


import matplotlib.pyplot as plt
class Plotting_Options:
    default_color_map = plt.cm.viridis

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Internal
ipython        = None
