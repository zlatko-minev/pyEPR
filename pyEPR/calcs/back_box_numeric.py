"""
Numerical diagonalization of quantum Hamiltonian and parameter
extraction.

@author: Phil Reinhold, Zlatko Minev, Lysander Christakis

Original code on black_box_hamiltonian and make_dispersive functions by Phil Reinhold
Revisions and updates by Zlatko Minev & Lysander Christakis
"""
# pylint: disable=invalid-name


from __future__ import print_function

from functools import reduce

import numpy as np

from .constants import Planck as h
from .constants import fluxQ, hbar
from .hamiltonian import MatrixOps

try:
    import qutip
    from qutip import basis, tensor
except (ImportError, ModuleNotFoundError):
    pass

__all__ = [
    "epr_numerical_diagonalization",
    "make_dispersive",
    "black_box_hamiltonian",
    "black_box_hamiltonian_nq",
    "cos_full_correction",
    "make_nonlinear_potential",
]

dot = MatrixOps.dot
cos_approx = MatrixOps.cos_approx
cos_full_correction = MatrixOps.cos_full_correction


def make_nonlinear_potential(V, phi_min=0.0, _eps=1e-5):
    r"""Build an EPR-compatible nonlinear potential from a user-supplied scalar function.

    **Background — EPR Hamiltonian split**

    In the EPR framework the Hamiltonian is decomposed as

    .. math::

        H = H_{\rm lin} - \sum_j \frac{E_{J,j}}{h}\,
            \underbrace{\Bigl[V(\varphi_{j,0} + \delta\varphi_j)
            - V(\varphi_{j,0})
            - \tfrac12 V''(\varphi_{j,0})\,\delta\varphi_j^2
            \Bigr] / |V''(\varphi_{j,0})|}_{\displaystyle \mathrm{nl}(\delta\varphi_j)}

    :math:`H_{\rm lin}` is built from the HFSS eigenmode frequencies, which already
    encode the **harmonic** (quadratic) Josephson energy via the linearised inductance
    :math:`L_j`.  The nonlinear correction ``nl(δφ)`` must therefore subtract exactly
    the constant and quadratic parts of *V* that are already captured.

    The function returned by :func:`make_nonlinear_potential` does this automatically:
    it evaluates *V* as an operator (via eigendecomposition of the phase operator),
    subtracts the constant :math:`V(\varphi_0)` and the quadratic term
    :math:`\tfrac12 V''(\varphi_0)\,\delta\varphi^2`, and normalises by
    :math:`|V''(\varphi_0)|` for consistency with the effective Josephson energy
    :math:`E_{J,\rm eff} = (\Phi_0/2\pi)^2 / L_j`.

    **Convention**: *V* is the normalised potential appearing in :math:`U_J = -E_J V(\varphi)`.
    For a standard Josephson junction :math:`V(\varphi) = \cos\varphi`, giving
    :math:`V''(0) = -1` and recovering the default :func:`cos_full_correction`.

    Parameters
    ----------
    V : callable
        Scalar potential function ``V(phi: float) -> float``.

        *phi* is the **reduced phase in the absolute frame** (not a fluctuation).
        Common choices:

        * ``np.cos`` — standard Josephson junction (reproduces :func:`cos_full_correction`)
        * ``lambda phi: np.cos(phi - phi_ext)`` — flux-biased junction at *phi_ext*
        * ``lambda phi: d_J * np.cos(phi - phi_ext)`` — asymmetric SQUID (scaled)

    phi_min : float, optional
        Phase at the potential minimum (equilibrium / bias point) in the same frame
        as *V*.  The EPR phase-fluctuation operator δφ is centred here.
        Default 0 (unbiased junction at its minimum).
    _eps : float, optional
        Finite-difference step for estimating :math:`V''(\varphi_0)`.  Default 1e-5.

    Returns
    -------
    callable
        ``nl(phi_op: qutip.Qobj) -> qutip.Qobj``, suitable as the
        ``non_linear_potential`` argument of :func:`epr_numerical_diagonalization`
        or :func:`black_box_hamiltonian`.

    Raises
    ------
    ValueError
        If :math:`|V''(\varphi_0)| < 10^{-12}`, indicating *phi_min* is not a
        proper quadratic extremum of *V*.

    Examples
    --------
    **Standard junction** — equivalent to the built-in default:

    >>> import numpy as np
    >>> nl = make_nonlinear_potential(np.cos)   # V(phi) = cos(phi), phi_min=0
    >>> # nl(phi_op) ≡ cos_full_correction(phi_op)

    **Flux-biased junction** at half a flux quantum (φ_ext = π/2):

    >>> phi_ext = np.pi / 2
    >>> nl = make_nonlinear_potential(lambda phi: np.cos(phi - phi_ext),
    ...                               phi_min=phi_ext)
    >>> f_ND, chi_ND = epr_numerical_diagonalization(
    ...     freqs, Ljs, phi_zpf, fock_trunc=20,
    ...     non_linear_potential=nl,
    ... )

    **Asymmetric SQUID** tuned to flux bias *phi_ext*, asymmetry *d = (Ej1-Ej2)/(Ej1+Ej2)*:

    >>> import numpy as np
    >>> phi_ext = 0.4  # external flux in reduced units
    >>> d = 0.1       # junction asymmetry
    >>> # Effective potential (see arXiv:2010.00620, §IV)
    >>> def V_squid(phi):
    ...     return np.cos(phi) * np.cos(phi_ext) + d * np.sin(phi) * np.sin(phi_ext)
    >>> phi_min = np.arctan(-d * np.tan(phi_ext))   # minimum of -Ej*V_squid
    >>> nl = make_nonlinear_potential(V_squid, phi_min=phi_min)

    Notes
    -----
    Uses eigendecomposition of φ_op to evaluate *V* as a matrix function — exact for
    any analytic *V* but slightly slower than the Taylor-series cosine default for
    small φ_zpf.  For production runs with standard junctions, prefer
    :func:`cos_full_correction` (``use_full_cos=True``) or the default ``cos_trunc``
    path; use :func:`make_nonlinear_potential` when the junction potential is not
    a simple cosine centred at zero.

    See Also
    --------
    cos_full_correction : Exact nonlinear potential for the standard Josephson junction.
    epr_numerical_diagonalization : Top-level solver; accepts ``non_linear_potential``.

    References
    ----------
    * Z. K. Minev *et al.*, arXiv:2010.00620 — EPR quantisation of Josephson circuits.
    * arXiv:2411.15039 — EPR analysis for very anharmonic superconducting circuits.
    """
    V0 = V(phi_min)
    Vpp = (V(phi_min + _eps) - 2.0 * V0 + V(phi_min - _eps)) / _eps ** 2
    abs_Vpp = abs(Vpp)

    if abs_Vpp < 1e-12:
        raise ValueError(
            f"V''(phi_min={phi_min}) ≈ 0: phi_min is not a quadratic extremum of V. "
            "Provide the correct potential minimum via phi_min."
        )

    def nl(phi_op):
        """Nonlinear EPR correction nl(δφ) = [V(φ_min+δφ) - V(φ_min) - V''(φ_min)/2·δφ²] / |V''(φ_min)|."""
        import qutip
        V_op = MatrixOps.apply_scalar_function(phi_op, lambda x: V(phi_min + x))
        I = qutip.qeye(phi_op.dims[0])
        return (V_op - V0 * I - (Vpp / 2.0) * phi_op ** 2) / abs_Vpp

    return nl


# ==============================================================================
# ANALYSIS FUNCTIONS
# ==============================================================================


def epr_numerical_diagonalization(
    freqs,
    Ljs,
    ϕzpf,
    cos_trunc=8,
    fock_trunc=9,
    use_1st_order=False,
    return_H=False,
    non_linear_potential=None,
    use_full_cos=False,
):
    """Numerical diagonalization of the EPR Hamiltonian.

    Parameters
    ----------
    freqs : array-like
        Linearised normal-mode frequencies in **GHz** (not radians), length M.
    Ljs : array-like
        Junction linearised inductances in **Henries**, length J.
    ϕzpf : array-like
        Reduced zero-point flux fluctuations, shape M × J.
    cos_trunc : int, optional
        Order of the Taylor-series truncation of cos(φ).  Only used when
        ``use_full_cos=False`` (default 8).  Ignored when ``use_full_cos=True``.
    fock_trunc : int, optional
        Fock-space truncation (number of levels per mode).  Default 9.
    use_1st_order : bool, optional
        Use first-order perturbation theory instead of full diagonalization.
        Default ``False``.
    return_H : bool, optional
        If ``True``, also return the Hamiltonian Qobj.  Default ``False``.
    non_linear_potential : callable, optional
        Custom nonlinear potential function.  Must accept a qutip Qobj φ (the
        reduced phase-fluctuation operator) and return a qutip Qobj equal to the
        **correction beyond the quadratic** of the junction potential:

        .. math::

            \\mathrm{nl}(\\delta\\varphi) =
              \\frac{V(\\varphi_0 + \\delta\\varphi) - V(\\varphi_0)
                    - \\tfrac{1}{2}V''(\\varphi_0)\\,\\delta\\varphi^2}
                   {|V''(\\varphi_0)|}

        Use :func:`make_nonlinear_potential` to build this callable from any
        scalar Python function.  Overrides both ``cos_trunc`` and ``use_full_cos``.
    use_full_cos : bool, optional
        If ``True``, use the **exact** matrix-exponential cosine
        ``cos(φ) = (e^{iφ} + e^{-iφ}) / 2`` instead of the truncated Taylor
        series.  Recommended for strongly anharmonic circuits such as
        **fluxonium**, where large zero-point phase fluctuations make the
        truncated expansion inaccurate.  Default ``False`` (truncated series).

    Returns
    -------
    f_ND : numpy.ndarray
        Dressed mode frequencies in **GHz**.
    χ_ND : numpy.ndarray
        Cross-Kerr / anharmonicity matrix in **MHz** (positive = red shift).
    Hs : qutip.Qobj
        Full Hamiltonian Qobj — only returned when ``return_H=True``.

    Notes
    -----
    **Choosing the right potential path:**

    * Small φ_zpf (transmon, φ_zpf ≲ 0.3): default ``cos_trunc=8`` is fast and
      accurate.
    * Large φ_zpf (fluxonium, φ_zpf ≳ 1): set ``use_full_cos=True`` to avoid
      truncation errors.  See arXiv:2411.15039.
    * Non-cosine or flux-biased junction: use ``non_linear_potential=make_nonlinear_potential(V)``.

    Examples
    --------
    Fluxonium with exact cosine:

    >>> f_ND, chi_ND = epr_numerical_diagonalization(
    ...     freqs, Ljs, phi_zpf, fock_trunc=25, use_full_cos=True
    ... )

    Flux-biased junction at φ_ext = π/4:

    >>> import numpy as np
    >>> phi_ext = np.pi / 4
    >>> nl = make_nonlinear_potential(lambda phi: np.cos(phi - phi_ext),
    ...                               phi_min=phi_ext)
    >>> f_ND, chi_ND = epr_numerical_diagonalization(
    ...     freqs, Ljs, phi_zpf, fock_trunc=15, non_linear_potential=nl
    ... )
    """

    freqs, Ljs, ϕzpf = map(np.array, (freqs, Ljs, ϕzpf))
    assert all(freqs < 1e6), "Please input the frequencies in GHz. \N{nauseated face}"
    assert all(
        Ljs < 1e-3
    ), "Please input the inductances in Henries. \N{nauseated face}"

    if non_linear_potential is None and use_full_cos:
        non_linear_potential = cos_full_correction

    Hs = black_box_hamiltonian(
        freqs * 1e9,
        Ljs.astype(float),
        fluxQ * ϕzpf,
        cos_trunc,
        fock_trunc,
        individual=use_1st_order,
        non_linear_potential=non_linear_potential,
    )
    f_ND, χ_ND, _, _ = make_dispersive(
        Hs, fock_trunc, ϕzpf, freqs, use_1st_order=use_1st_order
    )
    χ_ND = (
        -1 * χ_ND * 1e-6
    )  # convert to MHz, and flip sign so that down shift is positive

    return (f_ND, χ_ND, Hs) if return_H else (f_ND, χ_ND)


def black_box_hamiltonian(
    fs,
    ljs,
    fzpfs,
    cos_trunc=5,
    fock_trunc=8,
    individual=False,
    non_linear_potential=None,
):
    r"""Build the EPR Hamiltonian matrix in Fock space.

    Constructs

    .. math::

        H = H_{\rm lin} + H_{\rm nl}
          = \sum_m f_m \hat n_m
            - \sum_j \frac{E_{J,j}}{h}\,{\rm nl}(\hat\varphi_j)

    where :math:`\hat\varphi_j = \sum_m \phi_{\rm zpf,jm} (\hat a_m + \hat a_m^\dag)`
    (dimensionless reduced phase), and the nonlinear correction ``nl`` defaults to
    :func:`cos_approx` (:math:`\cos\varphi - 1 + \varphi^2/2`, the Taylor series
    starting at :math:`\varphi^4`).

    Parameters
    ----------
    fs : array-like
        Linearised normal-mode frequencies in **Hz**, length M.
    ljs : array-like
        Junction linearised inductances in **Henries**, length J.
    fzpfs : array-like
        Generalised (not reduced) zero-point flux fluctuations in Wb, shape M×J.
        Internally divided by ``fluxQ = Φ_0/(2π)`` to obtain reduced phases.
    cos_trunc : int
        Cosine Taylor-series truncation order (default 5). Ignored when
        *non_linear_potential* is given.
    fock_trunc : int
        Fock-space truncation per mode (default 8).
    individual : bool
        If ``True`` return ``[H_lin, H_nl]`` separately instead of their sum.
    non_linear_potential : callable or None
        Custom nonlinear potential: ``nl(φ_op: Qobj) -> Qobj``. Must return the
        correction **beyond the quadratic**, i.e.
        :math:`V(\varphi) - V(0) - \tfrac12 V''(0)\varphi^2` normalised by
        :math:`|V''(0)|`. Use :func:`make_nonlinear_potential` to construct this
        from any scalar Python function.

    Returns
    -------
    qutip.Qobj or [qutip.Qobj, qutip.Qobj]
        Full Hamiltonian H/h in Hz, or [H_lin/h, H_nl/h] when *individual=True*.
    """
    n_modes = len(fs)
    njuncs = len(ljs)
    fs, ljs, fzpfs = map(np.array, (fs, ljs, fzpfs))
    ejs = fluxQ**2 / ljs
    fjs = ejs / h

    fzpfs = np.transpose(fzpfs)  # Take from MxJ  to JxM

    assert np.isnan(fzpfs).any() == False, (
        "Phi ZPF has NAN, this is NOT allowed! Fix me. \n%s" % fzpfs
    )
    assert np.isnan(ljs).any() == False, "Ljs has NAN, this is NOT allowed! Fix me."
    assert np.isnan(fs).any() == False, "freqs has NAN, this is NOT allowed! Fix me."
    assert fzpfs.shape == (
        njuncs,
        n_modes,
    ), "incorrect shape for zpf array, {} not {}".format(fzpfs.shape, (njuncs, n_modes))
    assert fs.shape == (n_modes,), "incorrect number of mode frequencies"
    assert ejs.shape == (njuncs,), "incorrect number of qubit frequencies"

    def tensor_out(op, loc):
        "Make operator <op> tensored with identities at locations other than <loc>"
        op_list = [qutip.qeye(fock_trunc) for i in range(n_modes)]
        op_list[loc] = op
        return reduce(qutip.tensor, op_list)

    a = qutip.destroy(fock_trunc)
    ad = a.dag()
    n = qutip.num(fock_trunc)
    mode_fields = [tensor_out(a + ad, i) for i in range(n_modes)]
    mode_ns = [tensor_out(n, i) for i in range(n_modes)]

    def cos(x):
        return cos_approx(x, cos_trunc=cos_trunc)

    if non_linear_potential is None:
        non_linear_potential = cos

    linear_part = dot(fs, mode_ns)
    cos_interiors = [dot(fzpf_row / fluxQ, mode_fields) for fzpf_row in fzpfs]
    nonlinear_part = dot(-fjs, map(non_linear_potential, cos_interiors))
    if individual:
        return linear_part, nonlinear_part
    else:
        return linear_part + nonlinear_part


bbq_hmt = black_box_hamiltonian


def make_dispersive(
    H, fock_trunc, fzpfs=None, f0s=None, chi_prime=False, use_1st_order=False
):
    r"""
    Input: Hamiltonian Matrix.
        Optional: phi_zpfs and normal mode frequencies, f0s.
        use_1st_order : deprecated
    Output:
        Return dressed mode frequencies, chis, chi prime, phi_zpf flux (not reduced), and linear frequencies
    Description:
        Takes the Hamiltonian matrix `H` from bbq_hmt. It them finds the eigenvalues/eigenvectors and  assigns quantum numbers to them --- i.e., mode excitations,  such as, for instance, for three mode, :math:`|0,0,0\rangle` or :math:`|0,0,1\rangle`, which correspond to no excitations in any of the modes or one excitation in the 3rd mode, resp.    The assignment is performed based on the maximum overlap between the eigenvectors of H_full and H_lin.
        Based on the assignment of the excitations, the function returns the dressed mode frequencies :math:`\omega_m^\prime`, and the cross-Kerr matrix (including anharmonicities) extracted from the numerical diagonalization, as well as from 1st order perturbation theory.
        Note, the diagonal of the CHI matrix is directly the anharmonicity term.
    """
    if isinstance(H, (list, tuple)):  # [H_lin, H_nl] from individual=True
        [H_lin, H_nl] = H
        H = H_lin + H_nl
    else:  # make sure its a quantum object
        from qutip import Qobj

        if not isinstance(H, Qobj):  # Validate that the input is a Qobj instance.
            raise TypeError(
                "Please pass in either a list of Qobjs or a Qobj for the Hamiltonian"
            )
        # assert type(
        #    H) == qutip.qobj.Qobj, "Please pass in either a list of Qobjs or Qobj for the Hamiltonian"

    from .. import logger as _logger
    _logger.info("Starting diagonalization (%d×%d matrix)", H.shape[0], H.shape[0])
    evals, evecs = H.eigenstates()
    _logger.info("Diagonalization complete")
    evals -= evals[0]

    N = int(np.log(H.shape[0]) / np.log(fock_trunc))  # number of modes
    assert H.shape[0] == fock_trunc**N

    def fock_state_on(d):
        """d={mode number: # of photons}"""
        return qutip.tensor(
            *[qutip.basis(fock_trunc, d.get(i, 0)) for i in range(N)]
        )  # give me the value d[i]  or 0 if d[i] does not exist

    if use_1st_order:
        num_modes = N
        _logger.debug("Using 1st-order perturbation theory")

        def multi_index_2_vector(d, num_modes, fock_trunc):
            return tensor([basis(fock_trunc, d.get(i, 0)) for i in range(num_modes)])
            """this function creates a vector representation a given fock state given the data for excitations per
                        mode of the form d={mode number: # of photons}"""

        def find_multi_indices(fock_trunc):
            multi_indices = [
                {ind: item for ind, item in enumerate([i, j, k])}
                for i in range(fock_trunc)
                for j in range(fock_trunc)
                for k in range(fock_trunc)
            ]
            return multi_indices
            """this function generates all possible multi-indices for three modes for a given fock_trunc"""

        def get_expect_number(left, middle, right):
            result = left.dag() * middle * right
            # qutip 4: returns 1x1 Qobj; qutip 5: returns complex scalar
            return result.full()[0, 0] if hasattr(result, "full") else complex(result)
            """this function calculates the expectation value of an operator called "middle" """

        def get_basis0(fock_trunc, num_modes):
            multi_indices = find_multi_indices(fock_trunc)
            basis0 = [
                multi_index_2_vector(multi_indices[i], num_modes, fock_trunc)
                for i in range(len(multi_indices))
            ]
            evalues0 = [get_expect_number(v0, H_lin, v0).real for v0 in basis0]
            return multi_indices, basis0, evalues0
            """this function creates a basis of fock states and their corresponding eigenvalues"""

        def closest_state_to(vector0):
            def PT_on_vector(original_vector, original_basis, pertub, energy0, evalue):
                new_vector = 0 * original_vector
                for i in range(len(original_basis)):
                    if (energy0[i] - evalue) > 1e-3:
                        mel = original_basis[i].dag() * H_nl * original_vector
                        # qutip 4: 1x1 Qobj; qutip 5: complex scalar
                        mel = mel.full()[0, 0] if hasattr(mel, "full") else complex(mel)
                        new_vector += mel * original_basis[i] / (evalue - energy0[i])
                    else:
                        pass
                return (new_vector + original_vector) / (
                    new_vector + original_vector
                ).norm()
                """this function calculates the normalized vector with the first order correction term
                   from the non-linear hamiltonian """

            [multi_indices, basis0, evalues0] = get_basis0(fock_trunc, num_modes)
            evalue0 = get_expect_number(vector0, H_lin, vector0)
            vector1 = PT_on_vector(vector0, basis0, H_nl, evalues0, evalue0)

            # qutip 4: dag()*ket returns 1×1 Qobj with .norm(); qutip 5: returns complex scalar
            index = np.argmax([
                abs(r.full()[0, 0]) if hasattr(r, "full") else abs(r)
                for r in (vector1.dag() * evec for evec in evecs)
            ])
            return evals[index], evecs[index]

    else:

        def closest_state_to(s):
            def distance(s2):
                return np.abs((s.dag() * s2[1]))

            return max(zip(evals, evecs), key=distance)

    f1s = [closest_state_to(fock_state_on({i: 1}))[0] for i in range(N)]
    chis = [[0] * N for _ in range(N)]
    chips = [[0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i, N):
            d = {k: 0 for k in range(N)}  # put 0 photons in each mode (k)
            d[i] += 1
            d[j] += 1
            # load ith mode and jth mode with 1 photon
            fs = fock_state_on(d)
            ev, evec = closest_state_to(fs)
            chi = ev - (f1s[i] + f1s[j])
            chis[i][j] = chi
            chis[j][i] = chi

            if chi_prime:
                d[j] += 1
                fs = fock_state_on(d)
                ev, evec = closest_state_to(fs)
                chip = ev - (f1s[i] + 2 * f1s[j]) - 2 * chis[i][j]
                chips[i][j] = chip
                chips[j][i] = chip

    if chi_prime:
        return (
            np.array(f1s),
            np.array(chis),
            np.array(chips),
            np.array(fzpfs),
            np.array(f0s),
        )
    else:
        return np.array(f1s), np.array(chis), np.array(fzpfs), np.array(f0s)


def black_box_hamiltonian_nq(
    freqs, zmat, ljs, cos_trunc=6, fock_trunc=8, show_fit=False
):
    """
    N-Qubit version of bbq, based on the full Z-matrix
    Currently reproduces 1-qubit data, untested on n-qubit data
    Assume: Solve the model without loss in HFSS.
    """
    nf = len(freqs)
    nj = len(ljs)
    assert zmat.shape == (nf, nj, nj)

    imY = (1 / zmat[:, 0, 0]).imag
    # zeros where the sign changes from negative to positive

    (zeros,) = np.where((imY[:-1] <= 0) & (imY[1:] > 0))
    nz = len(zeros)

    imYs = np.array([1 / zmat[:, i, i] for i in range(nj)]).imag
    f0s = np.zeros(nz)
    slopes = np.zeros((nj, nz))
    import matplotlib.pyplot as plt

    # Fit a second order polynomial in the region around the zero
    # Extract the exact location of the zero and the associated slope
    # If you need better than second order fit, you're not sampling finely enough
    for i, z in enumerate(zeros):
        f0_guess = (freqs[z + 1] + freqs[z]) / 2
        zero_polys = np.polyfit(
            freqs[z - 1 : z + 3] - f0_guess, imYs[:, z - 1 : z + 3].transpose(), 2
        )
        zero_polys = zero_polys.transpose()
        f0s[i] = f0 = min(np.roots(zero_polys[0]), key=lambda r: abs(r)) + f0_guess
        for j, p in enumerate(zero_polys):
            slopes[j, i] = np.polyval(np.polyder(p), f0 - f0_guess)
        if show_fit:
            plt.plot(
                freqs[z - 1 : z + 3] - f0_guess,
                imYs[:, z - 1 : z + 3].transpose(),
                lw=1,
                ls="--",
                marker="o",
                label=str(f0),
            )
            p = np.poly1d(zero_polys[0, :])
            p2 = np.poly1d(zero_polys[1, :])
            plt.plot(
                freqs[z - 1 : z + 3] - f0_guess, p(freqs[z - 1 : z + 3] - f0_guess)
            )
            plt.plot(
                freqs[z - 1 : z + 3] - f0_guess, p2(freqs[z - 1 : z + 3] - f0_guess)
            )
            plt.legend(loc=0)

    zeffs = 2 / (slopes * f0s[np.newaxis, :])
    # Take signs with respect to first port
    zsigns = np.sign(zmat[zeros, 0, :])
    fzpfs = zsigns.transpose() * np.sqrt(hbar * abs(zeffs) / 2)

    H = black_box_hamiltonian(f0s, ljs, fzpfs, cos_trunc, fock_trunc)
    return make_dispersive(H, fock_trunc, fzpfs, f0s)


black_box_hamiltonian_nq = black_box_hamiltonian_nq
