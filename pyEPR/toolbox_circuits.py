#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:14:08 2019

Unit and variable conversions.

@author: Zlatko Minev
"""

from __future__ import division, print_function, absolute_import   # Python 2.7 and 3 compatibility
import numpy as np
import pandas as pd
#from collections import OrderedDict

from numpy import sqrt
from .toolbox import fluxQ, Planck, hbar, e_el, pi
#from scipy.constants import hbar, Planck, e as e_el, epsilon_0, pi

class Convert(object):
    '''
        Static container class for conversions of units and variables. 
        
        TEST CONVERSION:
        ```python
            from pyEPR.toolbox_circuits import Convert
        
            Lj_nH, Cs_fF = 11, 60
            Convert.transmon_print_all_params(Lj_nH, Cs_fF);

        ```
    '''
    
    @staticmethod
    def Ejj_from_Lj(Lj, units_Lj='nH'):
        ''' 
            Josephson Junction energy from Josephson inductance. 
            Returns in MHz
                
            $E_j = \phi_0^2 / L_J$
        '''
        Lj = 1.0*Lj # to float
        if units_Lj == 'nH':
            return (fluxQ**2.)/(Lj*1.0E-3*Planck)
        else:
            raise(NotImplementedError()) 

    @staticmethod
    def Ljj_from_Ejj(Ej, units_Ej='MHz'):
        ''' 
            Josephson Junction ind from Josephson energy in MHZ. 
            Returns in SI units of Henries
                
            $E_j = \phi_0^2 / L_J$
        '''
        Ej = 1.0*Ej # to float
        if units_Ej == 'MHz':
            return (fluxQ**2.)/(Ej*1.0E6*Planck)
        else:
            raise(NotImplementedError()) 
            
    @staticmethod
    def Ec_from_Cs(Cs, units_Lj='fF'):
        ''' 
            Charging energy 4Ec n^2, where n=Q/2e
            Returns in MHz
                
            $E_{C}=\frac{e^{2}}{2C}J$
        '''
        Cs = 1.0*Cs # to float
        if units_Lj == 'fF':
            return (e_el**2.)/(2.*Cs*1.0E-9*Planck)
        else:
            raise(NotImplementedError()) 
            
    @staticmethod
    def C_from_Ec(Ec, units_Ec='MHz'):
        ''' 
            Charging energy 4Ec n^2, where n=Q/2e
            
            Returns in SI units, in Farads.
                
            $E_{C}=\frac{e^{2}}{2C}J$
        '''
        Ec = 1.0*Ec # to float
        if units_Ec == 'MHz':
            return (e_el**2.)/(2.*Ec*1.0E6*Planck)
        else:
            raise(NotImplementedError()) 
            
    @staticmethod
    def ZPF_from_LC(L, C):
        '''
            Input units assumed to be identical
            
            Returns Phi ZPF in and Q_ZPF in NOT reduced units, but SI 
        '''
        Z = sqrt(L/C)
        return ( sqrt(hbar*Z/2.), sqrt(hbar/(2.*Z)) )  # Phi , Q 
    
    @staticmethod
    def Omega_from_LC(L, C):
        '''
            Calculate the resonant *angular* frequency
        '''
        return sqrt(1./(L*C))
    
    @staticmethod
    def transmon_get_all_params(Ej_MHz, Ec_MHz):
        '''Convinince func'''
        Ej, Ec         = Ej_MHz, Ec_MHz
        Lj_H, Cs_F     = Convert.Ljj_from_Ejj(Ej), Convert.C_from_Ec(Ec)  # SI units 
        Phi_ZPF, Q_ZPF = Convert.ZPF_from_LC(Lj_H, Cs_F)
        Omega_MHz      = sqrt(1./(Lj_H*Cs_F)) * 1E-6 # MHz
        f_MHz          = Omega_MHz / (2*pi)*1E-3
        Z_Ohms         = sqrt(Lj_H/Cs_F)
        phi_ZPF        = Phi_ZPF/fluxQ
        n_ZPF          = Q_ZPF / (2*e_el)
        return {'Ej_MHz'   : Ej_MHz,    'Ec_MHz': Ec_MHz,
                'Lj_H'     : Lj_H,      'Cs_F'  : Cs_F,
                'Lj_nH'    : Lj_H*1E9,  'Cs_fF' : Cs_F*1E15,
                'Phi_ZPF'  : Phi_ZPF,   'Q_ZPF' : Q_ZPF,
                'phi_ZPF'  : phi_ZPF,   'n_ZPF' : n_ZPF,
                'Omega_MHz': Omega_MHz,
                'f_MHz'    : f_MHz,
                'Z_Ohms'   : Z_Ohms,
               }

    @staticmethod
    def transmon_print_all_params(Lj_nH, Cs_fF):
        # Parameters - duplicates with transmon_get_all_params
        Ej, Ec     = Convert.Ejj_from_Lj(Lj_nH), Convert.Ec_from_Cs(Cs_fF) # MHz
        Lj_H, Cs_F = Convert.Ljj_from_Ejj(Ej), Convert.C_from_Ec(Ec)  # SI units 
        Phi_ZPF, Q_ZPF = Convert.ZPF_from_LC(Lj_H, Cs_F)
        Omega_MHz  = sqrt(1./(Lj_H*Cs_F)) * 1E-6 # MHz
        
        ### Print
        text = r"""
        \begin{align}
            L_J               &=%.1f \mathrm{\ nH}       &  C_\Sigma &=%.1f \mathrm{\ fF}   \\
            E_J               &=%.2f \mathrm{\ GHz}      &  E_C      &=%.0f \mathrm{\ MHz}  \\
            \omega_0  &=2\pi\times %.2f \mathrm{\ GHz}   &  Z_0 &= %.0f \mathrm{\ \Omega}   \\
            \phi_\mathrm{ZPF} &= %.2f \ \ \phi_0         &  n_\mathrm{ZPF} &=%.2f \ \ (2e)  \\ 
        \end{align}
        """ %(Lj_H*1E9, Cs_F*1E15, Ej/1E3, Ec,\
             Omega_MHz / (2*pi)*1E-3, sqrt(Lj_H/Cs_F),
             Phi_ZPF/fluxQ,   Q_ZPF / (2*e_el)) 
        
        from IPython.display import display, Math
        display(Math(text))

        return text



