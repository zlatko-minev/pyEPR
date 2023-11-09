"""
Basic calculations that apply in general .
"""

import numpy as np
from numpy import sqrt
from .. import logger

class CalcsBasic():

    @staticmethod
    def epr_to_zpf(Pmj, SJ, Ω, EJ):
        r'''
        Arguments, All as matrices (numpy arrays):
            :Pnj: MxJ energy-participation-ratio matrix, p_mj
            :SJ: MxJ sign matrix, s_mj
            :Ω: MxM diagonal matrix of frequencies (GHz, not radians, diagonal)
            :EJ: JxJ diagonal matrix matrix of Josephson energies (in same units as Om)

        RETURNS:
            reduced zpf  (in units of :math:`\phi_0`)
        '''
        (Pmj, SJ, Ω, EJ) = map(np.array, (Pmj, SJ, Ω, EJ))

        if (Pmj < 0).any():
            print('BAD!')
            logger.error(f"""The simulation is not converged!!! \N{nauseated face}
            Some of the energy participations are less than zero.
            This happens when some participations are tiny 10^-8 or less
            or when not enough passes have been taken. The Pmj matrix is
            {Pmj}""")

        # Technically, there the equation is hbar omega / 2J, but here we assume
        # that the hbar is absorbed in the units of omega, and omega and Ej have the same units.
        # PHI=np.zeros((3,3))
        # for m in range(3):
        #     for j in range(3):
        #         PHI[m,j] = SJ[m,j]*sqrt(PJ[m,j]*Om[m,m]/(2.*EJ[j,j]))

        return SJ * sqrt(0.5 * Ω @ Pmj @ np.linalg.inv(EJ))

    @staticmethod
    def epr_cap_to_nzpf(Pmj_cap, SJ, Ω, Ec):
        """
        Experimental. To be tested
        """
        (Pmj, SJ, Ω, EJ) = map(np.array, (Pmj_cap, SJ, Ω, Ec))
        return SJ * sqrt(Ω @ Pmj @ np.linalg.inv(Ec) /(4*4))
