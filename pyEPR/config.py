"""
Created on Fri Oct 30 14:21:45 2015

@author: Zaki, Zlatko
"""

#------------------------------------------------------------
# Directories 
root_dir = r'C:\data\\'


class Dissipation_params:
	''' Loss properties of various materials and surfaces '''

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


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Internal
ipython        = None