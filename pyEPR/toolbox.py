# -*- coding: utf-8 -*-
"""
Created on Sat Feb 04 09:32:46 2017

@author: Minev
"""
from __future__ import division, print_function, absolute_import   # Python 2.7 and 3 compatibility
import warnings
import numpy as np
import pandas as pd


### Constants
from collections import OrderedDict
from scipy.constants import hbar, Planck, e as e_el, epsilon_0, pi

fluxQ = hbar / (2*e_el) # reduced flux quantum

#==============================================================================
# Utility functions
#==============================================================================

def combinekw(kw1, kw2):
    ''' Copy kw1,  update with kw2, return result '''
    kw = kw1.copy()
    kw.update(kw2)
    return kw

def fact(n):
    ''' Factorial '''
    if n <= 1:
        return 1
    return n * fact(n-1)

def nck(n, k):
    ''' choose '''
    return fact(n)/(fact(k)*fact(n-k))

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


def sort_df_col(df):
    '''         sort by numerical int order    '''
    col_names = df.columns
    if np.all(col_names.map(isint)):
        return df[col_names.astype(int).sort_values().astype(str)]
    else:
        return df

def sort_Series_idx(sr):
    '''         sort by numerical int order    '''
    idx_names = sr.index
    if np.all(idx_names.map(isint)):
        return sr[idx_names.astype(int).sort_values().astype(str)]
    else:
        return sr



def get_instance_vars(obj, Forbidden=[]):
    VARS = {}
    for v in dir(obj):
        if not ((v.startswith('__')) or (v.startswith('_'))):
            print(v)
            if not callable(getattr(obj,v)):
                if not (v in Forbidden):
                    VARS[v] = getattr(obj,v)
    return VARS


def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used. See StackExchange"""
    def newFunc(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning) #turn off filter
        warnings.warn("Call to deprecated function {}.".format(func.__name__), category=DeprecationWarning, stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning) #reset filter
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc

#==============================================================================
# Matrix
#==============================================================================

def print_matrix(M, frmt = "{:7.2f}", append_row = ""):
    M = np.mat(M)
    for row in np.array(M.tolist()):
        print(' ', end='')
        for chi in row:
            print( frmt.format(chi), end='')
        print( append_row+"\n", end='')

def divide_diagonal_by_2(CHI0, div_fact = 2.):
    CHI = CHI0.copy();
    CHI[np.diag_indices_from(CHI)] /= div_fact
    return CHI

def print_NoNewLine(text):
    print((text), end='')

def print_color(text, style = 0, fg=24, bg = 43, newline = True):
    '''style 0..8;   fg  30..38;  bg  40..48'''
    format = ';'.join([str(style), str(fg), str(bg)])
    s = '\x1b[%sm %s \x1b[0m' % (format, text)
    if newline:
        print(s)
    else:
        print(s, end='')


#==============================================================================
#%%         Dataframe
#==============================================================================

def DataFrame_col_diff(PS, indx=0):
    ''' check weather the columns of a dataframe are equal,
        returns a T/F series of the row index that specifies which rows are differnt
        USE:
            PS[DataFrame_col_diff(PS)]
    '''
    R = []
    for i in range(PS.shape[1]-1):
        R += [ PS.iloc[:,i] == PS.iloc[:,i+1] ]
    if len(R) == 1:
        return np.logical_not(R[0])
    else:
        return np.logical_not(np.logical_and.reduce(R))


def xarray_unravel_levels(arr, names, my_convert = lambda x: x):
    ''' Takes in nested dict of dict of dataframes
        names : names of lists; you dont have to include the last two dataframe columns & rows, but you can to override them
        requires  xarray
    '''
    import xarray
    if type(arr) == pd.DataFrame:
        return xarray.DataArray(arr, dims = None if len(names)==0 else names)
    elif type(arr) in  [OrderedDict, dict]:
        return xarray.concat([xarray_unravel_levels(item, names[1:]) for k, item in arr.items()], pd.Index(arr.keys(), name=names[0]) )
    elif type(arr) == xarray.DataArray:
        return arr
    else:
        return my_convert(arr)

def robust_percentile(calc_data, ROBUST_PERCENTILE = 2.):
    '''
        analysis helper function
    '''
    vmin = np.percentile(calc_data, ROBUST_PERCENTILE)
    vmax = np.percentile(calc_data, 100 - ROBUST_PERCENTILE)
    return vmin, vmax




__all__  = ['hbar', 'e_el', 'epsilon_0', 'pi', 'fluxQ',
            'fact', 'nck', 'combinekw',
            'divide_diagonal_by_2',
            'sort_df_col', 'sort_Series_idx',
            'print_matrix', 'print_NoNewLine',
            'DataFrame_col_diff', 'xarray_unravel_levels', 'robust_percentile']

