"""
Hamiltonian and Matrix Operations.
"""
import pyEPR.calcs.quantum as qop
import numpy as np
import scipy.linalg as sla

from ..toolbox.pythonic import fact


class MatrixOps(object):

    @staticmethod
    def cos(op_cos_arg: np.ndarray):
        """
        Make cosine opertor matrix from arguemnt  op_cos_arg

            op_cos_arg (np.ndarray) : argumetn of the cosine
        """

        return 0.5*(sla.expm(1j*op_cos_arg) + sla.expm(-1j*op_cos_arg))

    @staticmethod
    def cos_approx(x, cos_trunc=5):
        """
        Create a Taylor series matrix approximation of the cosine, up to some order.
        """
        return sum((-1)**i * x**(2*i) / float(fact(2*i)) for i in range(2, cos_trunc + 1))

    @staticmethod
    def dot(ais, bis):
        """
        Dot product
        """
        return sum(ai*bi for ai, bi in zip(ais, bis))


class HamOps(object):

    @staticmethod
    def fock_state_on(d: dict, fock_trunc: int, N_modes: int):
        ''' d={mode number: # of photons} In the bare eigen basis
        '''
        # give me the value d[i]  or 0 if d[i] does not exist
        return qop.lkron([qop.basis(fock_trunc, d.get(i, 0)) for i in range(N_modes)])

    @staticmethod
    def closest_state_to(s: np.ndarray, energyMHz, evecs):
        """
        Returns the enery of the closest state to s
        """
        distance = lambda s2: np.linalg.norm(s.T.conj() * s2[1])
        return max(zip(energyMHz, evecs), key=distance)

    @staticmethod
    def closest_state_to_idx(s: np.ndarray, evecs):
        """
        Returns the index
        """
        distance = lambda s2: np.linalg.norm(s.T.conj() * s2[1])
        return max(zip(range(len(evecs)), evecs), key=distance)

    @staticmethod
    def identify_Fock_levels(fock_trunc: int, evecs,
                             N_modes=2,
                             Fock_max=4):
        """
        Return quantum numbers in terms of the undiagonalized eigenbasis.
        """
        #  to do: need to turn Fock_max into arb algo on each mdoe

        def fock_state_on(d):
            return HamOps.fock_state_on(d, fock_trunc, N_modes)

        def closest_state_to_idx(s):
            return HamOps.closest_state_to_idx(s, evecs)

        FOCKr = {}
        for d1 in range(Fock_max):
            for d2 in range(Fock_max):
                d = {0: d1, 1: d2}
                FOCKr[closest_state_to_idx(fock_state_on(d))[0]] = d
        return FOCKr
