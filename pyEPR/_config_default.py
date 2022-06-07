"""
 --- DO NOT MODIFY THIS FILE ---

Default configuration file for pyEPR

This file is NOT meant for users to modify.
Rather, a user should update any config settings they want
in a dictionary called CONFIG in a file called config.py

@author: Zlatko Minev and the pyEPR team
@date: Created on Fri Oct 30 14:21:45 2015
"""

import collections.abc
from . import Dict

# If we are reloading the package, then config will already be defined, then do not overwrite it.
__config_defined__ = 'config' in locals()


config = Dict( # pylint: disable=invalid-name

    # Folder to save result data to.
    root_dir=r'C:\data-pyEPR',
    save_format=r'%Y-%m-%d %H-%M-%S',

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

        # To save or not the mesh statistics from an HFSS run
        save_mesh_stats=True,
    ),

    epr = Dict(

        # Define the participation renormalization method
        # False : no extra renormalization to enforce
        #         can be more problematic for large pj, when sim isn't well converged
        # True or 1 : use enforcement of U_J_total to be U_mode-U_H
        #         can be more problematic for small pj, when sim isn't well converged
        # 2     : use enforcement of U_J_total to be U_mode-U_H (i.e., 1)
        #         only when the total participation is above a certain threshold
        #         preferred method.
        renorm_pj = 2,
    ),

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

    plotting=Dict(
        # Default color map for plotting. Better if made into a string name
        # taken from matplotlib.cm
        default_color_map='viridis',  # pylint: disable=no-member
    ),

    # Not to be used by the user. Just internal
    internal=Dict(

        # Are we using ipython
        ipython=None,

        # Error message for loading packages
        error_msg_missing_import="""\N{face with head-bandage}
        If you need a part of pyEPR that uses this package,
        then please install it. Then add it to the system path (if needed).
        See online setup instructions at
            http://www.github.com/zlatko-minev/pyEPR""",

        # Warn on missing import
        warn_missing_import=False,
    ),

    # Logging
    log=Dict(

        # '%(name)s - %(levelname)s - %(message)s\n   ::%(pathname)s:%(lineno)d: %(funcName)s\n')
        format='%(levelname)s %(asctime)s [%(funcName)s]: %(message)s',

        datefmt='%I:%M%p', #'%I:%M%p %Ss'

        level='INFO'
    )

)


def is_using_ipython():
    """Check if we're in IPython.

    Returns:
        bool -- True if ran in IPython
    """
    try:
        __IPYTHON__  # pylint: disable=undefined-variable, pointless-statement
        return True
    except NameError:
        return False


def update_recursive(d:collections.abc.Mapping, u:collections.abc.Mapping):
    """Recursive update of dictionaries.

    Arguments:
        d {collections.abc.Mapping} -- dict to overwrite
        u {collections.abc.Mapping} -- dict used to update

    Returns:
        same as d; Updated d
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_recursive(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def get_config():
    """Returns the config pointer.

    If the config is not yet loaded, it will load the default config and then
    update it with the _config_user.config dictionary.

    Else, it will just return the pointer to the above-updated config, which the
    user could have modified. The modifications will be kept.

    Returns:
        Dict : the config dictionary
    """
    if __config_defined__:
        #print('Config is already defined.') # not sure we ever make it here
        return config

    else:
        # Config is only loaded for the first time, set it up.
        #print('First time load of config')

        # Update with user config
        from . import _config_user
        _config = update_recursive(config, _config_user.config)

        # Add to config any bootup params
        config.internal.ipython = is_using_ipython()

        return config


__all__ = ['get_config']
