# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 19:30:12 2017

@author: Zlatko K. Minev
"""

from __future__ import division, print_function, absolute_import

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from .config import Plotting_Options
from matplotlib.colors import rgb2hex


# ==============================================================================
# Plotting
# ==============================================================================

def mpl_dpi(dpi=200):
    mpl.rcParams['figure.dpi'] = dpi
    mpl.rcParams['savefig.dpi'] = dpi

def legend_translucent(ax, values=[], loc=0, alpha=0.5, leg_kw={}):
    '''
    values = [ ["%.2f" %k for k in RES] ]

    Also, you can use the following:
    leg_kw = dict(fancybox   =True, fontsize = 9,
                  framealpha =0.5,  ncol     = 1)
    blah.plot().legend(**leg_kw )
    '''
    if ax.get_legend_handles_labels() == ([], []):
        return None

    leg = ax.legend(*values, loc=loc, fancybox=True, **leg_kw)
    leg.get_frame().set_alpha(alpha)
    return leg


def cmap_discrete(n, cmap_kw={}):
    ''' Discrete colormap.
        cmap_kw = dict(colormap = plt.cm.gist_earth, start = 0.05, stop = .95)

        cmap_kw
        -----------------------
        helix = True, Allows us to instead call helix from here
    '''
    if cmap_kw.pop('helix', False):
        return cmap_discrete_CubeHelix(n, helix_kw=cmap_kw)

    cmap_KW = dict(colormap=Plotting_Options.default_color_map,
                   start=0.05, stop=.95)
    cmap_KW.update(cmap_kw)

    return get_color_cycle(n+1, **cmap_KW)


def get_color_cycle(n, colormap=Plotting_Options.default_color_map, start=0., stop=1., format='hex'):
    pts = np.linspace(start, stop, n)
    if format == 'hex':
        colors = [rgb2hex(colormap(pt)) for pt in pts]
    return colors


def cmap_discrete_CubeHelix(n, helix_kw={}):
    '''
        https://github.com/jiffyclub/palettable/blob/master/demo/Cubehelix%20Demo.ipynb
        cube.show_discrete_image()
    '''
    from palettable import cubehelix
    helix_KW = dict(start_hue=240., end_hue=-300., min_sat=1., max_sat=2.5,
                    min_light=0.3, max_light=0.8, gamma=.9)
    helix_KW.update(helix_kw)
    cube = cubehelix.Cubehelix.make(n=n, **helix_KW)
    return cube.mpl_colors


def xarr_heatmap(fg, title=None, kwheat={}, fmt=('%.3f', '%.2f'), fig=None):
    ''' Needs seaborn and xarray'''
    fig = plt.figure() if fig == None else fig
    df = fg.to_pandas()
    # format indecies
    df.index = [float(fmt[0] % x) for x in df.index]
    df.columns = [float(fmt[1] % x) for x in df.columns]
    import seaborn as sns
    ax = sns.heatmap(df, annot=True, **kwheat)
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel(fg.dims[1])
    ax.set_ylabel(fg.dims[0])


__all__ = ['legend_translucent', 'cmap_discrete',
           'get_color_cycle', 'xarr_heatmap']

"""
Jupyter widgets:
--------------------------
Not seeing widgets: https://github.com/tqdm/tqdm/issues/451

    conda update tqdm
    # This might aleady work, will require a lot of updates, if not then do:
    conda install nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    jupyter nbextension enable --py widgetsnbextension

"""
