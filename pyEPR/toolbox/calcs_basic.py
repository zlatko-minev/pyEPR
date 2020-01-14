"""
Basic calculations that apply in general .
"""

import numpy as np
from numpy import sqrt
#from numpy.linalg import inv


class CalcsBasic():

    @staticmethod
    def epr_to_zpf(Pmj, SJ, Ω, EJ):
        r'''
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
