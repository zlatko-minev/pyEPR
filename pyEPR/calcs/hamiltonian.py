"""
Hamiltonian and Matrix Operations.
Hamiltonian operations heavily draw on qutip package.
This package must be installed for them to work.
"""

try:
    import qutip
    from qutip import Qobj  # basis, tensor,
except (ImportError, ModuleNotFoundError):
    Qobj = None
    pass

from ..toolbox.pythonic import fact


class MatrixOps(object):
    @staticmethod
    def cos(op_cos_arg: Qobj):
        """Exact cosine operator via matrix exponential: (e^{iφ} + e^{-iφ}) / 2.

        Parameters
        ----------
        op_cos_arg : qutip.Qobj
            Phase operator φ (Hermitian).

        Returns
        -------
        qutip.Qobj
            cos(φ) as a matrix in Fock space.
        """
        return 0.5 * ((1j * op_cos_arg).expm() + (-1j * op_cos_arg).expm())

    @staticmethod
    def cos_full_correction(op_cos_arg: Qobj):
        """Exact EPR nonlinear correction: cos(φ) - I + φ²/2 (no truncation).

        This is the infinite-order equivalent of :func:`cos_approx`.  In the
        EPR Hamiltonian the linear eigenfrequencies already account for the
        harmonic (φ²/2) Josephson energy, so the nonlinear potential that must
        be subtracted is the remainder:

        .. math::

            H_{\\mathrm{nl}} = -E_J \\bigl[\\cos(\\varphi) - 1 + \\tfrac{\\varphi^2}{2}\\bigr]
                             = -E_J \\sum_{n\\ge 2} \\frac{(-1)^n \\varphi^{2n}}{(2n)!}

        For weakly anharmonic circuits (transmon, φ_zpf ≲ 0.3) the truncated
        :func:`cos_approx` is equivalent and faster.  For strongly anharmonic
        circuits (fluxonium, φ_zpf ≳ 1) use this function to avoid systematic
        errors from the truncated series.

        Parameters
        ----------
        op_cos_arg : qutip.Qobj
            Phase operator φ (Hermitian), in the full Fock-space tensor product.

        Returns
        -------
        qutip.Qobj
            cos(φ) - I + φ²/2 as a matrix.

        See Also
        --------
        cos_approx : Truncated Taylor-series version (faster for small φ_zpf).

        References
        ----------
        arXiv:2411.15039 — EPR analysis for very anharmonic superconducting circuits.
        """
        import qutip
        cos_phi = MatrixOps.cos(op_cos_arg)
        # Build identity in the same tensor-product Hilbert space as op_cos_arg
        I = qutip.qeye(op_cos_arg.dims[0])
        return cos_phi - I + op_cos_arg**2 / 2

    @staticmethod
    def apply_scalar_function(op: Qobj, func) -> Qobj:
        """Evaluate a real scalar function on a Hermitian operator via eigendecomposition.

        For a Hermitian operator H with real eigenvalues λ_i and eigenvectors :math:`|i\\rangle`:

        .. math::

            f(H) = \\sum_i f(\\lambda_i) |i\\rangle\\langle i|

        This lets you evaluate *any* analytic scalar function—not just polynomials or
        exponentials—as an operator, which is needed for generic junction potentials.

        Parameters
        ----------
        op : qutip.Qobj
            Hermitian operator (e.g., the phase operator φ in Fock space).
        func : callable
            Real scalar function ``func(float) -> float``.

        Returns
        -------
        qutip.Qobj
            ``func(op)`` in the same Hilbert space.

        Examples
        --------
        >>> import qutip, numpy as np
        >>> from pyEPR.calcs.hamiltonian import MatrixOps
        >>> a = qutip.destroy(8)
        >>> phi = 0.3 * (a + a.dag())
        >>> cos_phi = MatrixOps.apply_scalar_function(phi, np.cos)
        """
        import numpy as np
        import qutip
        mat = op.full()
        evals, evecs = np.linalg.eigh(mat)          # Hermitian → real eigenvalues
        f_evals = np.vectorize(func)(evals)
        result = (evecs * f_evals) @ evecs.conj().T  # V @ diag(f) @ V†
        return qutip.Qobj(result, dims=op.dims)

    @staticmethod
    def cos_approx(x, cos_trunc=5):
        """
        Create a Taylor series matrix approximation of the cosine, up to some order.
        """
        return sum(
            (-1) ** i * x ** (2 * i) / float(fact(2 * i))
            for i in range(2, cos_trunc + 1)
        )

    @staticmethod
    def dot(ais, bis):
        """
        Dot product
        """
        return sum(ai * bi for ai, bi in zip(ais, bis))


class HamOps(object):
    @staticmethod
    def fock_state_on(d: dict, fock_trunc: int, N_modes: int):
        """d={mode number: # of photons} In the bare eigen basis"""
        # give me the value d[i]  or 0 if d[i] does not exist
        return qutip.tensor(
            *[qutip.basis(fock_trunc, d.get(i, 0)) for i in range(N_modes)]
        )

    @staticmethod
    def closest_state_to(s: Qobj, energyMHz, evecs):
        """
        Returns the energy of the closest state to s
        """

        def distance(s2):
            inner = s.dag() * s2[1]
            # qutip 4 returns a 1x1 Qobj (.norm()); qutip 5 returns a complex scalar (abs())
            return inner.norm() if hasattr(inner, "norm") else abs(inner)

        return max(zip(energyMHz, evecs), key=distance)

    @staticmethod
    def closest_state_to_idx(s: Qobj, evecs):
        """
        Returns the index
        """

        def distance(s2):
            inner = s.dag() * s2[1]
            # qutip 4 returns a 1x1 Qobj (.norm()); qutip 5 returns a complex scalar (abs())
            return inner.norm() if hasattr(inner, "norm") else abs(inner)

        return max(zip(range(len(evecs)), evecs), key=distance)

    @staticmethod
    def identify_Fock_levels(fock_trunc: int, evecs, N_modes=2, Fock_max=4):
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
