"""
Transmon calculations

"""
import math

import numpy as np
from numpy import pi, sqrt, exp
from numpy.linalg import inv

from .constants import e_el, fluxQ
from .convert import Convert
from ..toolbox.pythonic import divide_diagonal_by_2



class CalcsTransmon():
    """
    Common calculations and parameter reporting used for transmon qubits.
    """

    @staticmethod
    def dispersiveH_params_PT_O1(Pmj, Ωm, Ej):
        """
        First order PT on the 4th power of the JJ cosine.

        This function applied to an unfrustrated Josephson junction.

        Pmj : Matrix MxJ
        Ωm : GHz Matrix MxM
        Ej : GHz Matrix JxJ

        returns f_O1, χ_O1
        χ_O1 has diagonal divided by 2 so as to give true anharmonicity.

        Example use:
        ..codeblock python
            # PT_01: Calculate 1st order PT results
            f_O1, χ_O1 = Calc_basic.dispersiveH_params_PT_O1(Pmj, Ωm, Ej)
        """

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

    @staticmethod
    def transmon_get_all_params(Ej_MHz, Ec_MHz):
        """
        Linear harmonic oscillator approximation of transmon.
        Convenience func
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
        Convenience func
        """
        # Parameters - duplicates with transmon_get_all_params
        Ej, Ec = Convert.Ej_from_Lj(Lj_nH, 'nH', 'MHz'), Convert.Ec_from_Cs(
            Cs_fF, 'fF', 'MHz')  # MHz
        Lj_H, Cs_F = Convert.Lj_from_Ej(Ej, 'MHz', 'H'), Convert.Cs_from_Ec(
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
               Phi_ZPF/fluxQ, Q_ZPF / (2*e_el))

        from IPython.display import display, Math
        display(Math(text))

        return text

    @staticmethod
    def charge_dispersion_approx(m, Ec, Ej):
        """
        Use Eq. (2.5) of Koch's paper.
        """

        return sqrt(2./pi) * Ec * (-1.)**(m) * 2.**(4.*m+5.) * exp(-sqrt(8*Ej/Ec)) * (Ej/(2*Ec))**(m/2.+3./4.)\
            / math.factorial(m)
