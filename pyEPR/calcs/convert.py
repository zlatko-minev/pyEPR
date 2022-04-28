#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:14:08 2019

Unit and variable conversions.

@author: Zlatko Minev
"""

from __future__ import (absolute_import,  # Python 2.7 and 3 compatibility
                        division, print_function)

import numpy as np
import pandas as pd
from numpy import sqrt

from .basic import CalcsBasic
from .constants import (Planck, e_el, elementary_charge, fluxQ, hbar, pi, ħ, π,
                        ϕ0)


class Convert():
    '''
        Static container class for conversions of units and variables.

        TEST CONVERSION:

        .. code-block:: python

           from pyEPR.toolbox.conversions import Convert
           
           Lj_nH, Cs_fF = 11, 60
           Convert.transmon_print_all_params(Lj_nH, Cs_fF);
    '''
    # Known SI prefixed
    _prefix = {'y': -24,  # yocto
               'z': -21,  # zepto
               'a': -18,  # atto
               'f': -15,  # femto
               'p': -12,  # pico
               'n': -9,   # nano
               'u': -6,   # micro
               'm': -3,   # mili
               'c': -2,   # centi
               'd': -1,   # deci
               ' ': 0,
               'k': 3,    # kilo
               'M': 6,    # mega
               'G': 9,    # giga
               'T': 12,   # tera
               'P': 15,   # peta
               'E': 18,   # exa
               'Z': 21,   # zetta
               'Y': 24,   # yotta
               }

    # Known SI units
    _SI_units = ['H',  # Henries
                 'F',  # Farads
                 'Hz',  # Hertz
                 'Ohm',  # Ohms
                 'Ω',  # Ohms
                 'Wb'  # Webers
                 'J',  # Joules
                 'A'   # Amps
                 ]

    @staticmethod
    def toSI(number, from_units: str):
        r"""
        Convert from SI unit prefix to regular SI units
        If the from_units is ' ' or not in the prefix list,
        then the unit is assumed to be
        """
        if from_units in Convert._SI_units:
            from_units = ' '
        # else: we assume that the first letter is a prefix
        return number*(10**Convert._prefix.get(from_units[0]))

    @staticmethod
    def fromSI(number, from_units: str):
        r"""Convert a number with SI units, such as fF to F.

        Arguments:
            number {[numeric]} -- number
            from_units {str} -- string

        Returns:
            numeric number, with units expanded
        """
        if from_units in Convert._SI_units:
            from_units = ' '
        # else: we assume that the first letter is a prefix
        return number*(10**(-Convert._prefix.get(from_units[0])))

    @staticmethod
    def _convert_num(out_func, in_num, in_units, out_units):
        in_num = 1.0*in_num  # to float
        # convert units of input number
        in_num = Convert.toSI(in_num, in_units)
        out_num = out_func(in_num)  # Assume func processes all in SI units
        out_num = Convert.fromSI(out_num, out_units)
        return out_num

    @staticmethod
    def Ej_from_Lj(Lj, units_in='nH', units_out='MHz'):
        r'''
        Josephson Junction energy from Josephson inductance.
        Returns in MHz

        :math:`E_j = \phi_0^2 / L_J`
        '''
        return Convert._convert_num(
            # Plank to go from Joules to Hz
            lambda _Lj: Planck**-1 * (ϕ0**2)/_Lj,
            Lj, units_in, units_out)

    @staticmethod
    def Lj_from_Ej(Ej, units_in='MHz', units_out='nH'):
        r'''
        Josephson Junction ind from Josephson energy in MHZ.
        Returns in units of nano Henries by default

        :math:`E_j = \phi_0^2 / L_J`
        '''
        return Convert._convert_num(
            lambda _x: (ϕ0**2.)/(_x*Planck),  # Plank to go from Joules to Hz
            Ej, units_in, units_out)

    @staticmethod
    def Ic_from_Lj(Lj, units_in='nH', units_out='nA'):
        r'''
        Josephson Junction crit. curr from Josephson inductance.

        :math:`E_j = \phi_0^2 / L_J = \phi_0 I_C`
        '''
        return Convert._convert_num(
            lambda _x: ϕ0/_x,  # Plank to go from Joules to Hz
            Lj, units_in, units_out)

    @staticmethod
    def Lj_from_Ic(Lj, units_in='nA', units_out='nH'):
        r'''
        Josephson Junction crit. curr from Josephson inductance.

        :math:`E_j = \phi_0^2 / L_J = \phi_0 I_C`
        '''
        return Convert._convert_num(
            lambda _x: ϕ0/_x,  # Plank to go from Joules to Hz
            Lj, units_in, units_out)

    @staticmethod
    def Ec_from_Cs(Cs,  units_in='fF', units_out='MHz'):
        r'''
        Charging energy :math:`4E_c n^2`, where :math:`n=Q/2e`
        Returns in MHz

        :math:`E_{C}=\frac{e^{2}}{2C}J`
        '''
        return Convert._convert_num(
            # Plank to go from Joules to Hz
            lambda _x: Planck**-1 * (e_el**2.)/(2.*_x),
            Cs, units_in, units_out)

    @staticmethod
    def Cs_from_Ec(Ec, units_in='MHz', units_out='fF'):
        r'''
        Charging energy :math:`4E_c n^2`, where :math:`n=Q/2e`

        Returns in SI units, in Farads.

        :math:`E_{C}=\frac{e^{2}}{2C}J`
        '''
        return Convert._convert_num(
            # Plank to go from Joules to Hz
            lambda _x: (e_el**2.)/(2.*_x*Planck),
            Ec, units_in, units_out)

    @staticmethod
    def ZPF_from_LC(L, C):
        r'''
        Input units assumed to be identical

        Returns Phi ZPF in and Q_ZPF in NOT reduced units, but SI
        '''
        Z = sqrt(L/C)
        return (sqrt(hbar*Z/2.), sqrt(hbar/(2.*Z)))  # Phi , Q

    @staticmethod
    def ZPF_from_EPR(hfss_freqs, hfss_epr_, hfss_signs, hfss_Ljs,
                     Lj_units_in='H', to_df=False):
        r"""
        Parameters:
            Can be either Pandas or numpy arrays.

            hfss_freqs : HFSS Freqs. (standard units: GHz, but these will cancel with Ejs) (list/Series)
            hfss_epr : EPR ratio matrix, dim = M x J (2D array/DataFrame)
            hfss_signs : Sign matrix, dim = M x J  (2D array/DataFrame)
            hfss_Ljs : Assumed in Henries (see Lj_units_in). (list/Series)

            Lj_units_in : Default 'H' for Henries. Can change here.

        Returns:
            M x J matrix of reduced ZPF; i.e., scaled by reduced flux quantum.
            type: np.array
            and a tuple of matrices.

        Example use:
            ϕzpf, (Ωm, Ej, Pmj, Smj) = Convert.ZPF_from_EPR(hfss_freqs, hfss_epr, hfss_signs, hfss_Ljs, to_df=True)
        """

        hfss_freqs, hfss_epr, hfss_signs, hfss_Ljs = map(
            np.array, (hfss_freqs, hfss_epr_, hfss_signs, hfss_Ljs))

        Ωd = np.diagflat(hfss_freqs)
        Ej = Convert.Ej_from_Lj(
            hfss_Ljs, units_in=Lj_units_in, units_out='GHz')
        Ej = np.diagflat(Ej)

        ϕzpfs = CalcsBasic.epr_to_zpf(hfss_epr, hfss_signs, Ωd, Ej)

        if to_df:
            ϕzpfs = pd.DataFrame(
                ϕzpfs, columns='ϕ'+hfss_epr_.columns.values, index=hfss_epr_.index)

        return ϕzpfs, (Ωd, Ej, hfss_epr, hfss_signs)

    @staticmethod
    def Omega_from_LC(L, C):
        r'''
        Calculate the resonant *angular* frequency
        '''
        return sqrt(1./(L*C))
