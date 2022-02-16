"""
User configuration file.

The dictionary of options specified here overwrites the pyEPR default
config defined in _config_default.py

Do not edit `_config_default.py` directly. Rather, overwrite attributes here

@author: Your name goes here
"""

from . import Dict

config = Dict(

    # Folder to save result data to.
    # PLEASE CHANGE THIS
    root_dir=r'C:\data-pyEPR',

    # Loss properties of various materials and surfaces
    dissipation=Dict(

        ##################################################
        # Bulk dielectric
        # refs: https://arxiv.org/abs/1308.1743
        #       http://arxiv.org/pdf/1509.01854.pdf
        tan_delta_sapp=1e-6,  # tan(delta) for bulk surface
        epsi=10,    # dielectric

        ##################################################
        # Surface dielectric
        # ref: http://arxiv.org/pdf/1509.01854.pdf

        # Surface dielectric (dirt) thickness
        # units: meters
        th=3e-9,

        # Surface dielectric (dirt) constant
        # units: relative permittivity
        eps_r=10,

        # Surface dielectric (dirt) loss tangent
        # units: unitless, since this is tan(delta)
        tan_delta_surf=1e-3,

        ##################################################
        # Thin-film surface loss
        # units:  Ohms
        # ref:    https://arxiv.org/abs/1308.1743
        surface_Rs=250e-9,

        ##################################################
        # Seam current loss
        # units: per Ohm meter; i.e., seam conductance
        # ref:   http://arxiv.org/pdf/1509.01119.pdf
        gseam=1.0e3,
    ),

    ansys=Dict(
        # method_calc_P_mj sets the method used to calculate the participation ratio in eigenmode.
        # Valid values:
        # 	'line_voltage' : Uses the line voltage integral
        # 	'J_surf_mag'   : takes the avg. Jsurf over the rect. Make sure you have seeded
        # 					lots of tets here. I recommend starting with 4 across smallest dimension.
        # 					Multi-junction calculation of energy participation ratio matrix based on <I_J>.
        # 					Current is integrated average of J_surf by default: (zkm 3/29/16)
        # 					Will calculate the Pj matrix for the selected modes for the given junctions
        # 					junc_rect array & length of junctions
        method_calc_P_mj='line_voltage',

    ),

    plotting=Dict(
        # Default color map for plotting. Better if made into a string name
        # taken from matplotlib.cm
        default_color_map='viridis',  # pylint: disable=no-member
    ),
)


__all__ = ['config']
