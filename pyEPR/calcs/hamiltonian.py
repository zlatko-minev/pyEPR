"""
Hamiltonian and Matrix Operations.
Hamiltonian operations heavily draw on qutip package.
This package must be installed for them to work.
"""
try:
    import qutip
    from qutip import Qobj  # basis, tensor,
except (ImportError, ModuleNotFoundError):
    Qobj=None
    pass

from ..toolbox.pythonic import fact


class MatrixOps(object):

    @staticmethod
    def cos(op_cos_arg: Qobj):
        """
        Make cosine operator matrix from argument  op_cos_arg

            op_cos_arg (qutip.Qobj) : argument of the cosine
        """

        return 0.5*((1j*op_cos_arg).expm() + (-1j*op_cos_arg).expm())

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
        return qutip.tensor(*[qutip.basis(fock_trunc, d.get(i, 0))
                              for i in range(N_modes)])

    @staticmethod
    def closest_state_to(s: Qobj, energyMHz, evecs):
        """
        Returns the energy of the closest state to s
        """
        def distance(s2):
            return (s.dag() * s2[1]).norm()
        return max(zip(energyMHz, evecs), key=distance)

    @staticmethod
    def closest_state_to_idx(s: Qobj, evecs):
        """
        Returns the index
        """
        def distance(s2):
            return (s.dag() * s2[1]).norm()
        return max(zip(range(len(evecs)), evecs), key=distance)

    @staticmethod
    def identify_Fock_levels(fock_trunc: int, evecs,
                             N_modes=2,
                             Fock_max=4):
        """
        Return quantum numbers in terms of the undiagonalized eigenbasis.
        """
        #  to do: need to turn Fock_max into arb algo on each mode

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
