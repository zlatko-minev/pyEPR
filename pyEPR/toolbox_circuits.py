#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit and variable conversions.

@author: Zlatko K. Minev
"""

# Python 2.7 and 3 compatibility
from __future__ import division, print_function, absolute_import
import numpy as np
import pandas as pd
#from collections import OrderedDict

from numpy import sqrt
from .toolbox import ϕ0, fluxQ, Planck, hbar, e_el, pi, ħ, elementary_charge, π,\
    divide_diagonal_by_2
#from scipy.constants import hbar, Planck, e as e_el, epsilon_0, pi


from pyEPR.toolbox_circuits import *
import pint
ureg = pint.UnitRegistry()
ureg.enable_contexts('spectroscopy')
ureg.ϕ0 = ϕ0
#_c_jj = pint.Context('jj')
#_c_jj.add_transformation('([length] ** 2 * [mass] / [current] ** 2 / [time] ** 2)', '[current]', lambda ureg, x, **kw: ureg.ϕ0 / x)
#ureg.add_context(_c_jj)
#ureg.enable_contexts('jj')

class Convert(object):
    '''
        Static container class for conversions of units and variables.

        TEST CONVERSION:
        ```python
            from pyEPR.toolbox_circuits import Convert

            Lj_nH, Cs_fF = 11, 60
            Convert.transmon_print_all_params(Lj_nH, Cs_fF);

        ```

        TODO:  Superseed with pint
        with ureg.context('sp'):
            print(ureg.parse_expression('1 J').to('MHz'))
        See
            https://pint.readthedocs.io/en/0.9/contexts.html
            https://github.com/hgrecco/pint/blob/master/pint/default_en.txt
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
        """
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
        '''
            Josephson Junction energy from Josephson inductance.
            Returns in MHz

            $E_j = \phi_0^2 / L_J$
        '''
        return Convert._convert_num(
            # Plank to go from Joules to Hz
            lambda _Lj: Planck**-1 * (ϕ0**2)/_Lj,
            Lj, units_in, units_out)

    @staticmethod
    def Lj_from_Ej(Ej, units_in='MHz', units_out='nH'):
        '''
            Josephson Junction ind from Josephson energy in MHZ.
            Returns in units of nano Henries by default

            $E_j = \phi_0^2 / L_J$
        '''
        return Convert._convert_num(
            lambda _x: (ϕ0**2.)/(_x*Planck),  # Plank to go from Joules to Hz
            Ej, units_in, units_out)

    @staticmethod
    def Ic_from_Lj(Lj, units_in='nH', units_out='nA'):
        '''
            Josephson Junction crit. curr from Josephson inductance.

            $E_j = \phi_0^2 / L_J = \phi_0 I_C $
        '''
        return Convert._convert_num(
            lambda _x: ϕ0/_x,  # Plank to go from Joules to Hz
            Lj, units_in, units_out)

    @staticmethod
    def Lj_from_Ic(Lj, units_in='nA', units_out='nH'):
        '''
            Josephson Junction crit. curr from Josephson inductance.

            $E_j = \phi_0^2 / L_J = \phi_0 I_C $
        '''
        return Convert._convert_num(
            lambda _x: ϕ0/_x,  # Plank to go from Joules to Hz
            Lj, units_in, units_out)

    @staticmethod
    def Ec_from_Cs(Cs,  units_in='fF', units_out='MHz'):
        '''
            Charging energy 4Ec n^2, where n=Q/2e
            Returns in MHz

            $E_{C}=\frac{e^{2}}{2C}J$
        '''
        return Convert._convert_num(
            # Plank to go from Joules to Hz
            lambda _x: Planck**-1 * (e_el**2.)/(2.*_x),
            Cs, units_in, units_out)

    @staticmethod
    def Cs_from_Ec(Ec, units_in='MHz', units_out='fF'):
        '''
            Charging energy 4Ec n^2, where n=Q/2e

            Returns in SI units, in Farads.

            $E_{C}=\frac{e^{2}}{2C}J$
        '''
        return Convert._convert_num(
            # Plank to go from Joules to Hz
            lambda _x: (e_el**2.)/(2.*_x*Planck),
            Ec, units_in, units_out)

    @staticmethod
    def ZPF_from_LC(L, C):
        '''
            Input units assumed to be identical

            Returns Phi ZPF in and Q_ZPF in NOT reduced units, but SI
        '''
        Z = sqrt(L/C)
        return (sqrt(hbar*Z/2.), sqrt(hbar/(2.*Z)))  # Phi , Q

    @staticmethod
    def ZPF_from_EPR(hfss_freqs, hfss_epr_, hfss_signs, hfss_Ljs,
                     Lj_units_in='H', to_df=False):
        """
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
            and a tuple of matricies.

        Example use:
            ϕzpf, (Ωm, Ej, Pmj, Smj) = Convert.ZPF_from_EPR(hfss_freqs, hfss_epr, hfss_signs, hfss_Ljs, to_df=True)
        """

        hfss_freqs, hfss_epr, hfss_signs, hfss_Ljs = map(
            np.array, (hfss_freqs, hfss_epr_, hfss_signs, hfss_Ljs))

        Ωd = np.diagflat(hfss_freqs)
        Ej = Convert.Ej_from_Lj(
            hfss_Ljs, units_in=Lj_units_in, units_out='GHz')
        Ej = np.diagflat(Ej)

        ϕzpfs = Calcs_basic.epr_to_zpf(hfss_epr, hfss_signs, Ωd, Ej)

        if to_df:
            ϕzpfs = pd.DataFrame(
                ϕzpfs, columns='ϕ'+hfss_epr_.columns.values, index=hfss_epr_.index)

        return ϕzpfs, (Ωd, Ej, hfss_epr, hfss_signs)

    @staticmethod
    def Omega_from_LC(L, C):
        '''
            Calculate the resonant *angular* frequency
        '''
        return sqrt(1./(L*C))


class Calcs_basic(object):

    @staticmethod
    def dispersiveH_params_PT_O1(Pmj, Ωm, Ej):
        """
        First order PT on the 4th power of the JJ cosine.

        Pmj : Matrix MxJ
        Ωm : GHz Matrix MxM
        Ej : GHz Matrix JxJ

        returns f_O1, χ_O1
        χ_O1 has diagonal divided by 2 so as to give true anharmonicity.

        Example use:
        ..codeblock python
            f_O1, χ_O1 = Calc_basic.dispersiveH_params_PT_O1(Pmj, Ωm, Ej) # PT_01: Calculate 1st order PT results
        """
        from numpy.linalg import inv

        Pmj, Ωm, Ej = map(np.array, (Pmj, Ωm, Ej))

        assert Ωm.shape[0] == Ωm.shape[1]
        assert Ej.shape[0] == Ej.shape[1]
        assert Ωm.shape[1] == Pmj.shape[0]
        assert Pmj.shape[1] == Ej.shape[0]

        f_0 = np.diag(Ωm)

        χ_O1 = 0.25 * Ωm @ Pmj @ inv(Ej) @ Pmj.T @ Ωm * 1000.  # GHz to MHz

        f_O1 = f_0 - 0.5*np.ndarray.flatten(np.array(χ_O1.sum(1))) / \
            1000.  # 1st order PT expect freq to be dressed down by alpha

        # Make the diagonals alpha
        χ_O1 = divide_diagonal_by_2(χ_O1)

        return f_O1, χ_O1

    def epr_to_zpf(Pmj, SJ, Ω, EJ):
        '''
        INPUTS:
            All as matrices (numpy arrays)
            :Pnj: MxJ energy-participatuion-ratio matrix, p_mj
            :SJ: MxJ sign matrix, s_mj
            :Ω: MxM diagonal matrix of frequencies (GHz, not radians, diagonal)
            :EJ: JxJ diagonal matrix matrix of Josephson energies (in same units as Om)

        RETURNS:
            reduced zpf  (in units of $\phi_0$)
        '''
        (Pmj, SJ, Ω, EJ) = map(np.array, (Pmj, SJ, Ω, EJ))

        assert (Pmj > 0).any(), "ND -- p_{mj} are not all > 0; \n %s" % (Pmj)

        ''' technically, there the equation is hbar omega / 2J, but here we assume
        that the hbar is absrobed in the units of omega, and omega and Ej have the same units.
        PHI=np.zeros((3,3))
        for m in range(3):
            for j in range(3):
                PHI[m,j] = SJ[m,j]*sqrt(PJ[m,j]*Om[m,m]/(2.*EJ[j,j]))
        '''
        return SJ * sqrt(0.5 * Ω @ Pmj @ np.linalg.inv(EJ))

    @staticmethod
    def transmon_get_all_params(Ej_MHz, Ec_MHz):
        """
            Linear harmonic oscillator approximation of transmon.
            Convinince func
        """
        Ej, Ec = Ej_MHz, Ec_MHz
        Lj_H, Cs_F = Convert.Lj_from_Ej(
            Ej, 'MHz', 'H'), Convert.Cs_from_Ec(Ec, 'MHz', 'F')  # SI units
        Phi_ZPF, Q_ZPF = Convert.ZPF_from_LC(Lj_H, Cs_F)
        Omega_MHz = sqrt(1./(Lj_H*Cs_F)) * 1E-6  # MHz
        f_MHz = Omega_MHz / (2*pi)*1E-3
        Z_Ohms = sqrt(Lj_H/Cs_F)
        phi_ZPF = Phi_ZPF/fluxQ
        n_ZPF = Q_ZPF / (2*e_el)
        return {'Ej_MHz': Ej_MHz,    'Ec_MHz': Ec_MHz,
                'Lj_H': Lj_H,      'Cs_F': Cs_F,
                'Lj_nH': Lj_H*1E9,  'Cs_fF': Cs_F*1E15,
                'Phi_ZPF': Phi_ZPF,   'Q_ZPF': Q_ZPF,
                'phi_ZPF': phi_ZPF,   'n_ZPF': n_ZPF,
                'Omega_MHz': Omega_MHz,
                'f_MHz': f_MHz,
                'Z_Ohms': Z_Ohms,
                }

    @staticmethod
    def transmon_print_all_params(Lj_nH, Cs_fF):
        """
            Linear harmonic oscillator approximation of transmon.
            Convinince func
        """
        # Parameters - duplicates with transmon_get_all_params
        Ej, Ec = Convert.Ej_from_Lj(Lj_nH, 'nH', 'MHz'), Convert.Ec_from_Cs(
            Cs_fF, 'fF', 'MHz')  # MHz
        Lj_H, Cs_F = Convert.Lj_from_Ej(Ej, 'MHz', 'H'),     Convert.Cs_from_Ec(
            Ec, 'MHz', 'F')     # SI units
        Phi_ZPF, Q_ZPF = Convert.ZPF_from_LC(Lj_H, Cs_F)
        Omega_MHz = sqrt(1./(Lj_H*Cs_F)) * 1E-6  # MHz

        # Print
        text = r"""
        \begin{align}
            L_J               &=%.1f \mathrm{\ nH}       &  C_\Sigma &=%.1f \mathrm{\ fF}   \\
            E_J               &=%.2f \mathrm{\ GHz}      &  E_C      &=%.0f \mathrm{\ MHz}  \\
            \omega_0  &=2\pi\times %.2f \mathrm{\ GHz}   &  Z_0 &= %.0f \mathrm{\ \Omega}   \\
            \phi_\mathrm{ZPF} &= %.2f \ \ \phi_0         &  n_\mathrm{ZPF} &=%.2f \ \ (2e)  \\
        \end{align}
        """ % (Lj_H*1E9, Cs_F*1E15, Ej/1E3, Ec,
               Omega_MHz / (2*pi)*1E-3, sqrt(Lj_H/Cs_F),
               Phi_ZPF/fluxQ,   Q_ZPF / (2*e_el))

        from IPython.display import display, Math
        display(Math(text))

        return text


class Matrix_Operations(object):

    @staticmethod
    def to_cos(op_cos_arg):
        return 0.5*((1j*op_cos_arg).expm() + (-1j*op_cos_arg).expm())
