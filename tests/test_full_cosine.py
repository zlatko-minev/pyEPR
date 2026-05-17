"""
Tests for the full-cosine (exact matrix-exponential) diagonalization path,
relevant for strongly anharmonic circuits like fluxonium.

Physics note
------------
In the EPR Hamiltonian, the linear eigenfrequencies already account for the
harmonic (φ²/2) Josephson energy.  The nonlinear potential that is subtracted
is therefore the CORRECTION beyond quadratic:

    H_nl = -Ej * [cos(φ) - 1 + φ²/2]
         = -Ej * Σ_{n≥2} (-1)^n φ^{2n} / (2n)!

``cos_approx`` truncates this series.  ``cos_full_correction`` = cos(φ) - I + φ²/2
computes it exactly via matrix exponential — no truncation error.

Using the full cosine cos(φ) itself as the nonlinear potential would be WRONG
because it would double-count the harmonic term already in H_lin.

All tests run without Ansys.

Reference: arXiv:2411.15039 — EPR analysis for very anharmonic circuits.
"""
import numpy as np
import pytest


# ── helpers ──────────────────────────────────────────────────────────────────

def _transmon_params():
    """Transmon-like params: phi_zpf ≈ 0.15 (weakly anharmonic)."""
    freqs = np.array([5.0])       # GHz
    Ljs   = np.array([10e-9])     # H
    phi_zpf = np.array([[0.15]])  # shape (n_modes, n_junctions)
    return freqs, Ljs, phi_zpf


def _fluxonium_params():
    """Fluxonium-like params: phi_zpf ≈ 2.0 (strongly anharmonic)."""
    freqs = np.array([0.5])       # GHz
    Ljs   = np.array([500e-9])    # H
    phi_zpf = np.array([[2.0]])
    return freqs, Ljs, phi_zpf


# ── smoke tests ───────────────────────────────────────────────────────────────

class TestFullCosineSmoke:
    def test_truncated_cos_runs(self):
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _transmon_params()
        f_ND, chi_ND = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, cos_trunc=8, fock_trunc=9
        )
        assert f_ND is not None and chi_ND is not None

    def test_full_cos_flag_runs(self):
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _transmon_params()
        f_ND, chi_ND = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=9, use_full_cos=True
        )
        assert f_ND is not None and chi_ND is not None

    def test_full_cos_returns_finite_values(self):
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _transmon_params()
        f_ND, chi_ND = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=9, use_full_cos=True
        )
        assert np.all(np.isfinite(np.real(f_ND)))
        assert np.all(np.isfinite(np.real(chi_ND)))

    def test_full_cos_fluxonium_runs(self):
        """Full cosine must not crash for large phi_zpf (fluxonium regime)."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _fluxonium_params()
        f_ND, chi_ND = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=20, use_full_cos=True
        )
        assert np.all(np.isfinite(np.real(f_ND)))
        assert np.all(np.isfinite(np.real(chi_ND)))


# ── convergence: small phi_zpf → high-order truncated ≈ full ─────────────────

class TestFullCosConvergence:
    """For small phi_zpf, high-order truncated series converges to exact result."""

    def test_transmon_freqs_agree_high_trunc(self):
        """cos_trunc=16 and use_full_cos should agree to <0.1% for transmon."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _transmon_params()
        f_trunc, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, cos_trunc=16, fock_trunc=15
        )
        f_full, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=15, use_full_cos=True
        )
        np.testing.assert_allclose(
            np.real(f_full), f_trunc, rtol=1e-3,
            err_msg="High-order truncated and full cosine should agree for transmon",
        )

    def test_transmon_chi_agree_high_trunc(self):
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _transmon_params()
        _, chi_trunc = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, cos_trunc=16, fock_trunc=15
        )
        _, chi_full = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=15, use_full_cos=True
        )
        np.testing.assert_allclose(
            np.real(chi_full), chi_trunc, rtol=1e-2,
            err_msg="High-order truncated and full cosine chi should agree for transmon",
        )

    def test_low_trunc_differs_from_full(self):
        """Low truncation order gives inaccurate results even for transmon."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _transmon_params()
        _, chi_low = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, cos_trunc=4, fock_trunc=15
        )
        _, chi_full = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=15, use_full_cos=True
        )
        # At 4th order there should be a measurable difference (even for small phi_zpf)
        diff = abs(np.real(chi_full[0, 0]) - chi_low[0, 0])
        assert diff > 0, "cos_trunc=4 should differ from full cosine"


# ── large phi_zpf: truncated series diverges, full cosine converges ───────────

class TestFullCosFluxonium:
    def test_fluxonium_freqs_differ_from_low_trunc(self):
        """Low truncation gives wrong frequencies for fluxonium (phi_zpf~2)."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _fluxonium_params()
        f_low, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, cos_trunc=4, fock_trunc=25
        )
        f_full, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=25, use_full_cos=True
        )
        max_rel_diff = np.max(np.abs(np.real(f_full) - f_low) / np.abs(np.real(f_full)))
        assert max_rel_diff > 0.01, (
            f"Expected large disagreement for fluxonium at cos_trunc=4; "
            f"got only {max_rel_diff:.4f} relative diff"
        )

    def test_fluxonium_anharmonicity_sign(self):
        """Fluxonium anharmonicity (chi[0,0]) should be positive (red shift convention)."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs, Ljs, phi_zpf = _fluxonium_params()
        _, chi_full = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=25, use_full_cos=True
        )
        assert np.real(chi_full[0, 0]) > 0, "Fluxonium anharmonicity should be positive"


# ── MatrixOps.cos and cos_full_correction unit tests ─────────────────────────

class TestMatrixOpsCos:
    def test_cos_zero_is_identity(self):
        """cos(0) = I."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        n = 5
        zero_op = qutip.qzero(n)
        result = MatrixOps.cos(zero_op)
        np.testing.assert_allclose(result.full(), np.eye(n), atol=1e-12)

    def test_cos_is_hermitian(self):
        """cos(H) is Hermitian when H is Hermitian."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(8)
        H = 0.3 * (a + a.dag())
        result = MatrixOps.cos(H)
        assert (result - result.dag()).norm() < 1e-10

    def test_cos_agrees_with_numpy_on_diagonal(self):
        """For a diagonal operator, cos should match numpy.cos element-wise."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        angles = np.linspace(0, np.pi, 5)
        diag_op = qutip.Qobj(np.diag(angles))
        result = MatrixOps.cos(diag_op).full().real
        np.testing.assert_allclose(result, np.diag(np.cos(angles)), atol=1e-10)

    def test_cos_full_correction_zero_is_zero(self):
        """cos(0) - I + 0²/2 = I - I + 0 = 0."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        zero_op = qutip.qzero(6)
        result = MatrixOps.cos_full_correction(zero_op)
        np.testing.assert_allclose(result.full(), np.zeros((6, 6)), atol=1e-12)

    def test_cos_full_correction_is_hermitian(self):
        """cos_full_correction(H) is Hermitian when H is Hermitian."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(8)
        H = 0.5 * (a + a.dag())
        result = MatrixOps.cos_full_correction(H)
        assert (result - result.dag()).norm() < 1e-10

    def test_cos_full_correction_matches_series_small_arg(self):
        """cos(x) - I + x²/2 should match high-order cos_approx for small x."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(10)
        x = 0.1 * (a + a.dag())
        full = MatrixOps.cos_full_correction(x)
        approx = MatrixOps.cos_approx(x, cos_trunc=12)
        np.testing.assert_allclose(full.full(), approx.full(), atol=1e-8)


# ── flag behaviour ────────────────────────────────────────────────────────────

class TestFullCosFlag:
    def test_use_full_cos_produces_different_result_than_low_trunc(self):
        """use_full_cos=True differs from cos_trunc=4 at moderate phi_zpf."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        freqs = np.array([3.0])
        Ljs   = np.array([20e-9])
        phi_zpf = np.array([[0.8]])  # moderate — series converges slowly
        _, chi_low = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, cos_trunc=4, fock_trunc=12
        )
        _, chi_full = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=12, use_full_cos=True
        )
        assert abs(np.real(chi_full[0, 0]) - chi_low[0, 0]) > 0

    def test_explicit_non_linear_potential_overrides_use_full_cos(self):
        """Explicit non_linear_potential takes precedence over use_full_cos flag."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
        from pyEPR.calcs.hamiltonian import MatrixOps

        freqs, Ljs, phi_zpf = _transmon_params()
        f_full, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=9, use_full_cos=True
        )
        # Pass the same function explicitly — should give identical result
        f_explicit, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=9,
            non_linear_potential=MatrixOps.cos_full_correction,
            use_full_cos=True,  # flag ignored because non_linear_potential is set
        )
        np.testing.assert_allclose(np.real(f_full), np.real(f_explicit), rtol=1e-10)
