# -*- coding: utf-8 -*-
"""
Created on Sat Feb 04 09:32:46 2017

@author: Zlatko K. Minev, pyEPR ream
"""
from __future__ import division, print_function, absolute_import   # Python 2.7 and 3 compatibility
import platform         # Which OS we run
import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt

# Constants
from collections import OrderedDict
from ..calcs.constants import Planck, elementary_charge, epsilon_0, pi, π, ħ, ϕ0, e_el
from .. import Dict

# ==============================================================================
# Utility functions
# ==============================================================================

def combinekw(kw1, kw2):
    ''' Copy kw1,  update with kw2, return result '''
    kw = kw1.copy()
    kw.update(kw2)
    return kw


def isint(value):
    try:
        int(value)
        return True
    except ValueError:
        return False


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def floor_10(x):
    ''' round to nearest lower power of 10 c'''
    return 10.**(np.floor(np.log10(x)))


def fact(n):
    ''' Factorial '''
    if n <= 1:
        return 1
    return n * fact(n-1)


def nck(n, k):
    ''' choose '''
    return fact(n)/(fact(k)*fact(n-k))


def get_above_diagonal(M):
    ''' extract the values that are above the diagonal.
        Assumes square matrix
    '''
    return M[np.triu_indices(M.shape[0], k=1)]


def df_find_index(s: pd.Series, find, degree=2, ax=False):
    """
    Given a Pandas Series such as of freq with index Lj,
    find the Lj that would give the right frequency
    """
    max_ = max(s.index.values)
    min_ = min(s.index.values)
    if find <= max_ and find >= min_:
        # interpolate
        z = pd.Series(list(s.index.values)+[np.NaN], index=list(s) + [find])
        z = z.sort_index()
        z = z.interpolate()
        return z[find], z
    else:
        print('extrapolating')
        z = pd.Series(list(s.index.values), index=list(s))
        p = df_extrapolate(z, degree=degree, ax=False)
        value = p(find)
        return value, p


def df_interpolate_value(s: pd.Series, find, ax=False, method='index'):
    """
    Given a Pandas Series such as of freq with index Lj,
    find the freq that would correspond to Lj given a value not in the index
    """
    z = pd.Series(list(s) + [np.NaN], index=list(s.index.values)+[find])
    z = z.sort_index()
    z = z.interpolate(method=method)
    return z[find], z


def df_extrapolate(s, degree=2, ax=False, rng_scale=2.):
    """
    For a pandas series

    Returns np.poly1d
    """
    z = np.polyfit(s.index.values, s.values, degree)
    p = np.poly1d(z)

    if ax:
        if ax is True:
            ax = plt.gca()

        max_ = max(s.index.values)
        min_ = min(s.index.values)
        rng = max_ - min_
        xp = np.linspace(min_-rng_scale*rng, max_+rng_scale*rng, 100)
        ys = p(xp)
        ax.plot(xp, ys)
        s.plot(marker='o', ax=ax)

    return p


def df_regress_value(s: pd.Series, index, degree=2, ax=False, method='index',
                     rng_scale=2.):
    """
    for pandas series.
    calls either  df_interpolate_value or df_extrapolate
    """
    max_ = max(s.index.values)
    min_ = min(s.index.values)
    if index > max_ or index < min_:
        # print('extrapolate')
        p = df_extrapolate(s, degree=degree, ax=ax, rng_scale=rng_scale)
        value = p(index)
    else:
        value = df_interpolate_value(s, index, method=method)[0]
    return value


def series_of_1D_dict_to_multi_df(Uj_ind: pd.Series):
    df = pd.DataFrame(dict([(k, v) for k, v in Uj_ind.items()])).transpose()
    df.index.set_names(Uj_ind.index.names, inplace=True)
    return df


def sort_df_col(df):
    '''         sort by numerical int order    '''
    return df.sort_index(axis=1)

    # Buggy code, doesn't handles ints as inputs or floats as inputs
    col_names = df.columns
    if np.all(col_names.map(isint)):
        return df[col_names.astype(int).sort_values().astype(str)]
    elif np.all(col_names.map(isfloat)):
        # raises error in some cases
        return df[col_names.astype(float).sort_values().astype(str)]
    else:
        return df


def sort_Series_idx(sr):
    '''         sort by numerical int order    '''
    idx_names = sr.index
    if np.all(idx_names.map(isint)):
        return sr[idx_names.astype(int).sort_values().astype(str)]
    if np.all(idx_names.map(isfloat)):
        return sr[idx_names.astype(float).sort_values().astype(str)]
    else:
        return sr


def get_instance_vars(obj, Forbidden=[]):
    VARS = {}
    for v in dir(obj):
        if not (v.startswith('__') or v.startswith('_')):
            if not callable(getattr(obj, v)):
                # Added for using addict.Dict which is not callable. 
                if not isinstance(getattr(obj, v), Dict):
                    if not (v in Forbidden):
                        VARS[v] = getattr(obj, v)
    return VARS


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used. See StackExchange"""
    def newFunc(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)  # turn off filter
        warnings.warn("Call to deprecated function {}.".format(
            func.__name__), category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc


def info_str_platform():
    return '''

 System platform information:

    system   : %s
    node     : %s
    release  : %s
    machine  : %s
    processor: %s
    summary  : %s
    version  : %s

 Python platform information:

    version  : %s (implem: %s)
    compiler : %s

    ''' % (
        platform.system(),
        platform.node(),
        platform.release(),
        platform.machine(),
        platform.processor(),
        platform.platform(),
        platform.version(),
        platform.python_version(), platform.python_implementation(),
        platform.python_compiler())


# ==============================================================================
# Matrix
# ==============================================================================


def print_matrix(M, frmt="{:7.2f}", append_row=""):
    M = np.mat(M)
    for row in np.array(M.tolist()):
        print(' ', end='')
        for chi in row:
            print(frmt.format(chi), end='')
        print(append_row+"\n", end='')


def divide_diagonal_by_2(CHI0, div_fact=2.):
    CHI = CHI0.copy()
    CHI[np.diag_indices_from(CHI)] /= div_fact
    return CHI


def print_NoNewLine(text):
    print((text), end='')


def print_color(text, style=0, fg=24, bg=43, newline=True):
    '''For newer, see pc (or Print_colors)
       style 0..8;   fg  30..38;  bg  40..48
     '''
    format = ';'.join([str(style), str(fg), str(bg)])
    s = '\x1b[%sm %s \x1b[0m' % (format, text)
    if newline:
        print(s)
    else:
        print(s, end='')


class Print_colors:
    '''Colors class:reset all colors with colors.reset; two
    sub classes fg for foreground
    and bg for background; use as colors.subclass.colorname.
    i.e. colors.fg.red or colors.bg.green also, the generic bold, disable,
    underline, reverse, strike through,
    and invisible work with the main class i.e. colors.bold
    https://www.geeksforgeeks.org/print-colors-python-terminal/

    Example use:
    ..codeblock python
        print(colors.bg.green, "adgd", colors.fg.red, "dsgdsg")
        print(colors.bg.lightgrey, "dsgsd", colors.fg.red, "sdgsd")
    '''
    reset = '\033[0m'
    bold = '\033[01m'
    disable = '\033[02m'
    underline = '\033[04m'
    reverse = '\033[07m'
    strikethrough = '\033[09m'
    invisible = '\033[08m'

    class fg:
        black = '\033[30m'
        red = '\033[31m'
        green = '\033[32m'
        orange = '\033[33m'
        blue = '\033[34m'
        purple = '\033[35m'
        cyan = '\033[36m'
        lightgrey = '\033[37m'
        darkgrey = '\033[90m'
        lightred = '\033[91m'
        lightgreen = '\033[92m'
        yellow = '\033[93m'
        lightblue = '\033[94m'
        pink = '\033[95m'
        lightcyan = '\033[96m'

    class bg:
        black = '\033[40m'
        red = '\033[41m'
        green = '\033[42m'
        orange = '\033[43m'
        blue = '\033[44m'
        purple = '\033[45m'
        cyan = '\033[46m'
        lightgrey = '\033[47m'


pc = Print_colors

# ==============================================================================
# %%         Dataframe
# ==============================================================================


def DataFrame_col_diff(PS, indx=0):
    ''' check weather the columns of a dataframe are equal,
        returns a T/F series of the row index that specifies which rows are different
        USE:
            PS[DataFrame_col_diff(PS)]
    '''
    R = []
    for i in range(PS.shape[1]-1):
        R += [PS.iloc[:, i] == PS.iloc[:, i+1]]
    if len(R) == 1:
        return np.logical_not(R[0])
    else:
        return np.logical_not(np.logical_and.reduce(R))  # pylint: disable=no-member


def DataFrame_display_side_by_side(*args, do_display=True):
    '''
    from pyEPR.toolbox.pythonic import display_dfs
    https://stackoverflow.com/questions/38783027/jupyter-notebook-display-two-pandas-tables-side-by-side
    '''
    from IPython.display import display_html
    html_str = ''
    for df in args:
        html_str += df.to_html()
    text = html_str.replace('table', 'table style="display:inline"')
    if do_display:
        display_html(text, raw=True)
    return text


display_dfs = DataFrame_display_side_by_side


def xarray_unravel_levels(arr, names, my_convert=lambda x: x):
    ''' Takes in nested dict of dict of dataframes
        names : names of lists; you dont have to include the last two dataframe columns & rows, but you can to override them
        requires  xarray
    '''
    import xarray  # pylint: disable=import-error
    if type(arr) == pd.DataFrame:
        return xarray.DataArray(arr, dims=None if len(names) == 0 else names)
    elif type(arr) in [OrderedDict, dict]:
        return xarray.concat([xarray_unravel_levels(item, names[1:]) for k, item in arr.items()], pd.Index(arr.keys(), name=names[0]))
    elif type(arr) == xarray.DataArray:
        return arr
    else:
        return my_convert(arr)


def robust_percentile(calc_data, ROBUST_PERCENTILE=2.):
    '''
        analysis helper function
    '''
    vmin = np.percentile(calc_data, ROBUST_PERCENTILE)
    vmax = np.percentile(calc_data, 100 - ROBUST_PERCENTILE)
    return vmin, vmax


__all__ = ['fact', 'nck', 'combinekw',
           'divide_diagonal_by_2', 'df_find_index',
           'sort_df_col', 'sort_Series_idx',
           'print_matrix', 'print_NoNewLine',
           'DataFrame_col_diff', 'xarray_unravel_levels', 'robust_percentile']
