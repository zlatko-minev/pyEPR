"""
Hamiltonian and Matrix Operations.
Hamiltonian operations heavily draw on qutip package.
This package must be installded for them to work.
"""
try:
    import qutip
    from qutip import QObj #  basis, tensor,
except (ImportError, ModuleNotFoundError):
    pass

from ..toolbox.pythonic import fact

class MatrixOps(object):

    @staticmethod
    def cos(op_cos_arg:QObj):
        """
        Make cosine opertor matrix from arguemnt  op_cos_arg

            op_cos_arg (qutip.QObj) : argumetn of the cosine
        """

        return 0.5*( (1j*op_cos_arg).expm() + (-1j*op_cos_arg).expm() )

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
    def fock_state_on(d, fock_trunc, N_modes):
        ''' d={mode number: # of photons} In the bare eigen basis
        '''
        return qutip.tensor(*[qutip.basis(fock_trunc, d.get(i, 0)) for i in range(N_modes)])  # give me the value d[i]  or 0 if d[i] does not exist

    @staticmethod
    def closest_state_to(s, energyMHz, evecs):
        def distance(s2):
            return (s.dag() * s2[1]).norm()
        return max(zip(energyMHz, evecs), key=distance)

    @staticmethod
    def closest_state_to_idx(s, evecs):
        def distance(s2):
            return (s.dag() * s2[1]).norm()
        return max(zip(range(len(evecs)), evecs), key=distance)

    @staticmethod
    def identify_Fock_levels(fock_trunc, #= my_params['Nmax']
                         evecs,
                         N_modes = 2,
                         Fock_max = 4):
        fock_state_on        = lambda d: HamOps.fock_state_on(d, fock_trunc, N_modes)
        closest_state_to_idx = lambda s: HamOps.closest_state_to_idx(s, evecs)

        FOCKr = {} #TODO:# need to turn into arb algo
        for d1 in range(Fock_max):
            for d2 in range(Fock_max):
                d = {0:d1, 1:d2}
                FOCKr[closest_state_to_idx(fock_state_on(d))[0]] = d
        return FOCKr
