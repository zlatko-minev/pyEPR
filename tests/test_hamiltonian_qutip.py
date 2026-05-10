"""
Tests for qutip-dependent Hamiltonian construction and diagonalization.

Covers the code paths in:
  pyEPR/calcs/hamiltonian.py   — MatrixOps, HamOps
  pyEPR/calcs/back_box_numeric.py — black_box_hamiltonian, make_dispersive,
                                     epr_numerical_diagonalization

All tests use purely synthetic inputs and run without any HFSS connection.
Designed to catch qutip 4→5 API regressions (e.g. ket.dag()*ket now returns
a complex scalar in qutip 5, not a 1x1 Qobj).
"""
import numpy as np
import pytest

qutip = pytest.importorskip("qutip", reason="qutip not installed")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

FOCK = 6  # small truncation for fast tests


def simple_1mode_hamiltonian(fock_trunc=FOCK, freq_hz=5e9, lj_h=10e-9, zpf=0.15):
    """Return H (Qobj) for a 1-mode, 1-junction system."""
    from pyEPR.calcs.back_box_numeric import black_box_hamiltonian
    from pyEPR.calcs.constants import fluxQ

    fs = np.array([freq_hz])
    ljs = np.array([lj_h])
    fzpfs = np.array([[zpf * fluxQ]])  # generalized ZPF
    return black_box_hamiltonian(fs, ljs, fzpfs, cos_trunc=6, fock_trunc=fock_trunc)


def simple_2mode_hamiltonian(fock_trunc=FOCK):
    """Return H (Qobj) for a 2-mode, 1-junction system."""
    from pyEPR.calcs.back_box_numeric import black_box_hamiltonian
    from pyEPR.calcs.constants import fluxQ

    fs = np.array([5e9, 6.5e9])
    ljs = np.array([10e-9])
    fzpfs = np.array([[0.15 * fluxQ], [0.01 * fluxQ]])
    return black_box_hamiltonian(fs, ljs, fzpfs, cos_trunc=6, fock_trunc=fock_trunc)


# ---------------------------------------------------------------------------
# MatrixOps
# ---------------------------------------------------------------------------

class TestMatrixOps:
    def test_cos_approx_identity_at_zero(self):
        """cos_approx(0) ≈ identity (missing constant term from Taylor, but
        Taylor of cos(x) around 0 has leading term 1, which cos_approx omits
        — it returns only the non-trivial part starting at x^2)."""
        from pyEPR.calcs.hamiltonian import MatrixOps
        # Passing a zero operator: result should be a zero Qobj
        zero_op = qutip.qzero(FOCK)
        result = MatrixOps.cos_approx(zero_op, cos_trunc=5)
        # All terms are proportional to zero_op^(2i), so result is zero
        assert isinstance(result, qutip.Qobj)
        assert np.allclose(result.full(), 0)

    def test_cos_approx_returns_hermitian(self):
        """cos_approx of a Hermitian operator should be Hermitian."""
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(FOCK)
        x = a + a.dag()  # position-like, Hermitian
        result = MatrixOps.cos_approx(x, cos_trunc=6)
        assert isinstance(result, qutip.Qobj)
        diff = (result - result.dag()).norm()
        assert diff < 1e-10, f"cos_approx result not Hermitian, |H - H†| = {diff}"

    def test_cos_exact_hermitian(self):
        """MatrixOps.cos of a Hermitian operator should be Hermitian."""
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(FOCK)
        x = 0.1 * (a + a.dag())
        result = MatrixOps.cos(x)
        assert isinstance(result, qutip.Qobj)
        diff = (result - result.dag()).norm()
        assert diff < 1e-10

    def test_cos_exact_vs_approx_close_for_small_arg(self):
        """For small argument, exact cos and Taylor approx should agree.

        cos_approx starts at i=2 (x^4/4! term), intentionally omitting the
        constant (1) and quadratic (-x^2/2) terms which are handled elsewhere
        in the Hamiltonian.  The full cosine is 1 - x^2/2 + cos_approx(x).
        """
        from pyEPR.calcs.hamiltonian import MatrixOps
        from pyEPR.toolbox.pythonic import fact
        a = qutip.destroy(FOCK)
        x = 0.05 * (a + a.dag())  # small argument
        exact = MatrixOps.cos(x)
        # Reconstruct full Taylor series: 1 - x^2/2 + higher terms (cos_approx)
        full_approx = qutip.qeye(FOCK) - x**2 / float(fact(2)) + MatrixOps.cos_approx(x, cos_trunc=8)
        diff = (exact - full_approx).norm() / exact.norm()
        assert diff < 1e-6, f"cos exact vs approx differ by {diff:.2e}"

    def test_dot_product(self):
        """MatrixOps.dot should sum pairwise products correctly."""
        from pyEPR.calcs.hamiltonian import MatrixOps
        n = qutip.num(FOCK)
        I = qutip.qeye(FOCK)
        result = MatrixOps.dot([2.0, 3.0], [n, I])
        expected = 2.0 * n + 3.0 * I
        assert (result - expected).norm() < 1e-12


# ---------------------------------------------------------------------------
# HamOps — qutip 5 compatibility critical here
# ---------------------------------------------------------------------------

class TestHamOps:
    def test_fock_state_on_single_mode(self):
        """fock_state_on({0:n}) should give |n>."""
        from pyEPR.calcs.hamiltonian import HamOps
        ket0 = HamOps.fock_state_on({0: 0}, FOCK, N_modes=1)
        ket1 = HamOps.fock_state_on({0: 1}, FOCK, N_modes=1)
        # |0> and |1> should be orthogonal
        inner = ket0.dag() * ket1
        val = inner.norm() if hasattr(inner, "norm") else abs(inner)
        assert val < 1e-12

    def test_fock_state_on_two_mode(self):
        """|1,0> and |0,1> should be orthogonal; each should have unit norm."""
        from pyEPR.calcs.hamiltonian import HamOps
        s10 = HamOps.fock_state_on({0: 1, 1: 0}, FOCK, N_modes=2)
        s01 = HamOps.fock_state_on({0: 0, 1: 1}, FOCK, N_modes=2)
        inner = s10.dag() * s01
        val = inner.norm() if hasattr(inner, "norm") else abs(inner)
        assert val < 1e-12
        assert abs(s10.norm() - 1.0) < 1e-12
        assert abs(s01.norm() - 1.0) < 1e-12

    def test_closest_state_to_exact_eigenstate(self):
        """Given an exact Fock state as eigenvector, closest_state_to finds it."""
        from pyEPR.calcs.hamiltonian import HamOps
        # Simple 1-mode number operator — eigenstates are Fock states
        n_op = qutip.num(FOCK)
        evals, evecs = n_op.eigenstates()

        target = qutip.basis(FOCK, 2)  # |2>
        energy, vec = HamOps.closest_state_to(target, evals, evecs)
        assert abs(energy - 2.0) < 1e-10, f"Expected energy 2, got {energy}"

    def test_closest_state_to_idx_exact(self):
        """closest_state_to_idx returns the index of the best-matching eigenstate."""
        from pyEPR.calcs.hamiltonian import HamOps
        n_op = qutip.num(FOCK)
        evals, evecs = n_op.eigenstates()

        target = qutip.basis(FOCK, 3)  # |3>
        idx, _ = HamOps.closest_state_to_idx(target, evecs)
        assert idx == 3, f"Expected index 3, got {idx}"

    def test_identify_fock_levels_linear_hamiltonian(self):
        """For a linear 2-mode Hamiltonian, Fock level assignment should be trivial."""
        from pyEPR.calcs.hamiltonian import HamOps
        fock_trunc = 5
        I = qutip.qeye(fock_trunc)
        n = qutip.num(fock_trunc)
        # H = n ⊗ I + I ⊗ n  (equal mode frequencies)
        H = qutip.tensor(n, I) + qutip.tensor(I, n)
        _, evecs = H.eigenstates()

        fock_map = HamOps.identify_Fock_levels(fock_trunc, evecs, N_modes=2, Fock_max=3)
        # Ground state (index of |0,0>) should map to {0:0, 1:0}
        assert {0: 0, 1: 0} in fock_map.values()


# ---------------------------------------------------------------------------
# black_box_hamiltonian
# ---------------------------------------------------------------------------

class TestBlackBoxHamiltonian:
    def test_returns_qobj(self):
        """black_box_hamiltonian should return a qutip.Qobj."""
        H = simple_1mode_hamiltonian()
        assert isinstance(H, qutip.Qobj)

    def test_dimension_1mode(self):
        """1-mode fock_trunc=6 → H should be 6×6."""
        H = simple_1mode_hamiltonian(fock_trunc=6)
        assert H.shape == (6, 6)

    def test_dimension_2mode(self):
        """2-mode fock_trunc=5 → H should be 25×25."""
        H = simple_2mode_hamiltonian(fock_trunc=5)
        assert H.shape == (25, 25)

    def test_hermitian(self):
        """Hamiltonian must be Hermitian (H = H†)."""
        H = simple_1mode_hamiltonian()
        diff = (H - H.dag()).norm()
        assert diff < 1e-10, f"|H - H†| = {diff:.2e}"

    def test_hermitian_2mode(self):
        """2-mode Hamiltonian must be Hermitian."""
        H = simple_2mode_hamiltonian()
        diff = (H - H.dag()).norm()
        assert diff < 1e-10

    def test_eigenvalues_real(self):
        """All eigenvalues of a Hermitian operator should be real."""
        H = simple_1mode_hamiltonian()
        evals = H.eigenenergies()
        assert np.allclose(evals.imag, 0, atol=1e-10), "Eigenvalues have imaginary part"

    def test_ground_state_minimal(self):
        """Ground state energy should be the minimum eigenvalue."""
        H = simple_1mode_hamiltonian()
        evals = H.eigenenergies()
        assert evals[0] == pytest.approx(min(evals), rel=1e-10)

    def test_individual_mode_returns_tuple(self):
        """individual=True should return (H_lin, H_nl) tuple."""
        from pyEPR.calcs.back_box_numeric import black_box_hamiltonian
        from pyEPR.calcs.constants import fluxQ
        fs = np.array([5e9])
        ljs = np.array([10e-9])
        fzpfs = np.array([[0.15 * fluxQ]])
        result = black_box_hamiltonian(fs, ljs, fzpfs, individual=True)
        assert isinstance(result, tuple) and len(result) == 2
        H_lin, H_nl = result
        assert isinstance(H_lin, qutip.Qobj)
        assert isinstance(H_nl, qutip.Qobj)

    def test_nan_in_fzpfs_raises(self):
        """NaN in fzpfs should raise AssertionError."""
        from pyEPR.calcs.back_box_numeric import black_box_hamiltonian
        from pyEPR.calcs.constants import fluxQ
        fzpfs = np.array([[np.nan * fluxQ]])
        with pytest.raises(AssertionError):
            black_box_hamiltonian(np.array([5e9]), np.array([10e-9]), fzpfs)


# ---------------------------------------------------------------------------
# make_dispersive
# ---------------------------------------------------------------------------

class TestMakeDispersive:
    def test_output_shapes_1mode(self):
        """make_dispersive: 1-mode system → f1s shape (1,), chis shape (1,1)."""
        from pyEPR.calcs.back_box_numeric import black_box_hamiltonian, make_dispersive
        from pyEPR.calcs.constants import fluxQ
        H = simple_1mode_hamiltonian()
        phi_zpf = np.array([[0.15]])
        f0s = np.array([5.0])
        f1s, chis, _, _ = make_dispersive(H, fock_trunc=FOCK, fzpfs=phi_zpf, f0s=f0s)
        assert f1s.shape == (1,)
        assert np.array(chis).shape == (1, 1)

    def test_output_shapes_2mode(self):
        """make_dispersive: 2-mode system → f1s shape (2,), chis shape (2,2)."""
        from pyEPR.calcs.back_box_numeric import make_dispersive
        H = simple_2mode_hamiltonian(fock_trunc=FOCK)
        phi_zpf = np.array([[0.15], [0.01]])
        f0s = np.array([5.0, 6.5])
        f1s, chis, _, _ = make_dispersive(H, fock_trunc=FOCK, fzpfs=phi_zpf, f0s=f0s)
        assert f1s.shape == (2,)
        assert np.array(chis).shape == (2, 2)

    def test_chi_matrix_symmetric(self):
        """Cross-Kerr matrix should be symmetric."""
        from pyEPR.calcs.back_box_numeric import make_dispersive
        H = simple_2mode_hamiltonian(fock_trunc=FOCK)
        phi_zpf = np.array([[0.15], [0.01]])
        f0s = np.array([5.0, 6.5])
        _, chis, _, _ = make_dispersive(H, fock_trunc=FOCK, fzpfs=phi_zpf, f0s=f0s)
        chis = np.array(chis)
        assert np.isclose(chis[0, 1].real, chis[1, 0].real, rtol=1e-6)

    def test_accepts_list_input(self):
        """make_dispersive should accept [H_lin, H_nl] list (individual=True path)."""
        from pyEPR.calcs.back_box_numeric import black_box_hamiltonian, make_dispersive
        from pyEPR.calcs.constants import fluxQ
        fs = np.array([5e9])
        ljs = np.array([10e-9])
        fzpfs = np.array([[0.15 * fluxQ]])
        # Must use the same fock_trunc in both calls so dimensions match
        H_list = black_box_hamiltonian(fs, ljs, fzpfs, fock_trunc=FOCK, individual=True)
        phi_zpf = np.array([[0.15]])
        f0s = np.array([5.0])
        # Should not raise
        f1s, chis, _, _ = make_dispersive(H_list, fock_trunc=FOCK, fzpfs=phi_zpf, f0s=f0s)
        assert f1s.shape == (1,)

    def test_rejects_non_qobj(self):
        """Passing a numpy array should raise TypeError."""
        from pyEPR.calcs.back_box_numeric import make_dispersive
        with pytest.raises(TypeError):
            make_dispersive(np.eye(FOCK), fock_trunc=FOCK)


# ---------------------------------------------------------------------------
# Full pipeline: epr_numerical_diagonalization
# ---------------------------------------------------------------------------

class TestEprNumericalDiagonalizationFull:
    def test_return_H_flag(self):
        """return_H=True should return (f_ND, chi_ND, H) triple."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        result = epr_numerical_diagonalization(
            np.array([5.0]), np.array([10e-9]), np.array([[0.15]]),
            cos_trunc=6, fock_trunc=FOCK, return_H=True
        )
        assert len(result) == 3
        f_ND, chi_ND, H = result
        assert isinstance(H, qutip.Qobj)

    def test_custom_nonlinear_potential(self):
        """Passing a custom non_linear_potential should run without error."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        from pyEPR.calcs.hamiltonian import MatrixOps

        def quartic(x):
            return x**4 / 24.0  # trivial test potential

        f_ND, chi_ND = epr_numerical_diagonalization(
            np.array([5.0]), np.array([10e-9]), np.array([[0.15]]),
            cos_trunc=6, fock_trunc=FOCK,
            non_linear_potential=quartic
        )
        assert f_ND.shape == (1,)

    def test_physical_transmon_anharmonicity_sign(self):
        """
        Transmon anharmonicity should be positive in pyEPR sign convention
        (down-shift = positive value in chi_ND diagonal).
        Verified for a range of physically reasonable ZPF values.
        """
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        for zpf in [0.10, 0.15, 0.20]:
            _, chi_ND = epr_numerical_diagonalization(
                np.array([5.0]), np.array([10e-9]), np.array([[zpf]]),
                cos_trunc=8, fock_trunc=9
            )
            assert chi_ND[0, 0].real > 0, (
                f"Anharmonicity should be positive for zpf={zpf}, got {chi_ND[0,0].real}"
            )

    def test_larger_zpf_gives_larger_anharmonicity(self):
        """Larger ZPF → stronger coupling → larger anharmonicity."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        _, chi_small = epr_numerical_diagonalization(
            np.array([5.0]), np.array([10e-9]), np.array([[0.05]]),
            cos_trunc=8, fock_trunc=9
        )
        _, chi_large = epr_numerical_diagonalization(
            np.array([5.0]), np.array([10e-9]), np.array([[0.20]]),
            cos_trunc=8, fock_trunc=9
        )
        assert chi_large[0, 0].real > chi_small[0, 0].real
