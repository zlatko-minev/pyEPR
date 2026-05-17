"""
Tests for make_nonlinear_potential — generic scalar potential API.

Physics background
------------------
In the EPR Hamiltonian H = H_lin + H_nl, H_lin already encodes the harmonic
(quadratic) Josephson energy via the linearised inductance Lj.  The nonlinear
correction that must be added is

    nl(δφ) = [V(φ_min + δφ) - V(φ_min) - V''(φ_min)/2 · δφ²] / |V''(φ_min)|

For the standard Josephson junction V(φ) = cos(φ), φ_min = 0:
    V(0) = 1, V''(0) = -1  →  nl(δφ) = cos(δφ) - 1 + δφ²/2

which is exactly cos_full_correction.  make_nonlinear_potential recovers this
automatically for any valid V.

All tests run without Ansys.
"""
import numpy as np
import pytest


# ── helpers ──────────────────────────────────────────────────────────────────

def _transmon_params():
    return np.array([5.0]), np.array([10e-9]), np.array([[0.15]])


def _fluxonium_params():
    return np.array([0.5]), np.array([500e-9]), np.array([[2.0]])


# ── MatrixOps.apply_scalar_function ──────────────────────────────────────────

class TestApplyScalarFunction:
    def test_cos_matches_matrix_exponential(self):
        """apply_scalar_function(phi, np.cos) should equal MatrixOps.cos(phi)."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(8)
        phi = 0.3 * (a + a.dag())
        via_expm = MatrixOps.cos(phi)
        via_eig  = MatrixOps.apply_scalar_function(phi, np.cos)
        np.testing.assert_allclose(via_eig.full(), via_expm.full(), atol=1e-10)

    def test_identity_function_returns_operator(self):
        """apply_scalar_function(H, lambda x: x) should be a no-op."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(6)
        phi = 0.4 * (a + a.dag())
        result = MatrixOps.apply_scalar_function(phi, lambda x: x)
        np.testing.assert_allclose(result.full(), phi.full(), atol=1e-12)

    def test_constant_function_gives_scaled_identity(self):
        """apply_scalar_function(H, lambda x: 3.0) should be 3·I."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(5)
        phi = 0.2 * (a + a.dag())
        result = MatrixOps.apply_scalar_function(phi, lambda x: 3.0)
        np.testing.assert_allclose(result.full(), 3.0 * np.eye(5), atol=1e-12)

    def test_hermitian_preserved(self):
        """f(H) is Hermitian when H is Hermitian and f is real."""
        import qutip
        from pyEPR.calcs.hamiltonian import MatrixOps
        a = qutip.destroy(8)
        phi = 0.5 * (a + a.dag())
        result = MatrixOps.apply_scalar_function(phi, np.cos)
        assert (result - result.dag()).norm() < 1e-10


# ── make_nonlinear_potential: basic physics ───────────────────────────────────

class TestMakeNonlinearPotentialPhysics:
    def test_cosine_matches_cos_full_correction(self):
        """make_nonlinear_potential(cos) must agree with cos_full_correction.

        Both compute the same mathematical object via different numerical paths
        (eigendecomposition vs. matrix exponential), so we allow ~1 ppm tolerance.
        """
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential, cos_full_correction
        a = qutip.destroy(10)
        phi = 0.4 * (a + a.dag())
        nl_generic = make_nonlinear_potential(np.cos)(phi)
        nl_exact   = cos_full_correction(phi)
        np.testing.assert_allclose(nl_generic.full(), nl_exact.full(), atol=1e-6)

    def test_zero_operator_gives_zero(self):
        """nl(0) = V(phi_min+0) - V(phi_min) - 0 = 0 for any V."""
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        nl = make_nonlinear_potential(np.cos)
        zero_op = qutip.qzero(6)
        result = nl(zero_op)
        np.testing.assert_allclose(result.full(), np.zeros((6, 6)), atol=1e-12)

    def test_result_is_hermitian(self):
        """nl(H) must be Hermitian when H is Hermitian."""
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        a = qutip.destroy(8)
        phi = 0.5 * (a + a.dag())
        nl = make_nonlinear_potential(np.cos)
        result = nl(phi)
        assert (result - result.dag()).norm() < 1e-10

    def test_leading_term_is_quartic(self):
        """For small φ, nl ≈ φ⁴/24; the quadratic term is subtracted by construction.

        We verify this by comparing nl with the leading φ⁴ approximation on the
        diagonal of the phase operator (where comparison is exact).
        """
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        # Diagonal phase operator: eigenvalues are just the diagonal entries
        angles = np.array([0.0, 0.05, 0.1, 0.15, 0.2])
        diag_op = qutip.Qobj(np.diag(angles))
        nl = make_nonlinear_potential(np.cos)
        result = nl(diag_op).full().real
        # Expected: cos(x) - 1 + x^2/2 evaluated at each diagonal
        expected = np.diag(np.cos(angles) - 1 + angles**2 / 2)
        # Tolerance accounts for finite-difference error in V'' estimation (~eps^2~1e-9)
        np.testing.assert_allclose(result, expected, atol=1e-7)
        # Verify leading term really is phi^4/24 (no phi^2 contribution)
        x = angles[2]  # x = 0.1
        leading = x**4 / 24
        actual = np.cos(x) - 1 + x**2 / 2
        np.testing.assert_allclose(actual, leading, rtol=0.01)  # 1% — cos6 term is small

    def test_normalisation_independent_of_overall_scale(self):
        """Scaling V by a constant must not change nl (Ej_eff absorbs it)."""
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        a = qutip.destroy(8)
        phi = 0.3 * (a + a.dag())
        nl1 = make_nonlinear_potential(np.cos)(phi)
        nl2 = make_nonlinear_potential(lambda x: 3.0 * np.cos(x))(phi)
        np.testing.assert_allclose(nl1.full(), nl2.full(), atol=1e-8)

    def test_invalid_phi_min_raises(self):
        """V''(phi_min) ≈ 0 should raise ValueError."""
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        # V(phi) = phi has V''=0 everywhere
        with pytest.raises(ValueError, match="quadratic extremum"):
            make_nonlinear_potential(lambda phi: phi, phi_min=0.0)


# ── flux-biased junction ──────────────────────────────────────────────────────

class TestFluxBiasedJunction:
    """
    For V(phi) = cos(phi - phi_ext) biased at phi_min = phi_ext,
    nl should match the unbiased cosine nl because the expansion at the
    minimum is always cos(delta_phi).
    """

    def test_biased_matches_unbiased_at_minimum(self):
        """cos(φ - φ_ext) expanded at φ_ext == cos(δφ) expanded at 0."""
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        phi_ext = 0.4
        a = qutip.destroy(10)
        # The phase fluctuation operator (small ZPF)
        phi_op = 0.15 * (a + a.dag())
        nl_unbiased = make_nonlinear_potential(np.cos, phi_min=0.0)(phi_op)
        nl_biased   = make_nonlinear_potential(
            lambda phi: np.cos(phi - phi_ext), phi_min=phi_ext
        )(phi_op)
        np.testing.assert_allclose(nl_biased.full(), nl_unbiased.full(), atol=1e-8)

    def test_biased_half_flux_is_sin_correction(self):
        """At phi_ext = pi/2, V(phi)=cos(phi - pi/2)=sin(phi); nl is well-defined."""
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        phi_ext = np.pi / 2
        a = qutip.destroy(8)
        phi_op = 0.15 * (a + a.dag())
        nl = make_nonlinear_potential(
            lambda phi: np.cos(phi - phi_ext), phi_min=phi_ext
        )
        result = nl(phi_op)
        assert np.all(np.isfinite(result.full()))
        assert (result - result.dag()).norm() < 1e-10


# ── asymmetric SQUID ──────────────────────────────────────────────────────────

class TestAsymmetricSQUID:
    """
    Asymmetric SQUID at external flux phi_ext, asymmetry d = (Ej1-Ej2)/(Ej1+Ej2):
        V(phi) = cos(phi)*cos(phi_ext) + d*sin(phi)*sin(phi_ext)

    The effective Josephson energy is Ej_eff = Ej_total * |V''(phi_min)|.
    At d=0 this reduces to a symmetric SQUID; at phi_ext=0 it reduces to a JJ.
    """

    def test_d0_matches_standard_cosine(self):
        """Symmetric SQUID (d=0) at phi_ext=0 → standard cos correction."""
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential, cos_full_correction
        phi_ext, d = 0.0, 0.0
        def V_squid(phi):
            return np.cos(phi) * np.cos(phi_ext) + d * np.sin(phi) * np.sin(phi_ext)
        phi_min = 0.0  # minimum of -Ej*V_squid at phi_ext=0, d=0

        a = qutip.destroy(10)
        phi_op = 0.2 * (a + a.dag())
        nl_squid = make_nonlinear_potential(V_squid, phi_min=phi_min)(phi_op)
        nl_std   = cos_full_correction(phi_op)
        np.testing.assert_allclose(nl_squid.full(), nl_std.full(), atol=1e-7)

    def test_asymmetric_squid_returns_finite_hermitian(self):
        """Asymmetric SQUID at moderate bias gives finite, Hermitian result."""
        import qutip
        from pyEPR.calcs.back_box_numeric import make_nonlinear_potential
        phi_ext, d = 0.3, 0.1
        def V_squid(phi):
            return np.cos(phi) * np.cos(phi_ext) + d * np.sin(phi) * np.sin(phi_ext)
        # phi_min: maximum of V_squid (i.e., minimum of -Ej*V_squid)
        phi_min = np.arctan(-d * np.tan(phi_ext)) % np.pi

        a = qutip.destroy(10)
        phi_op = 0.2 * (a + a.dag())
        nl = make_nonlinear_potential(V_squid, phi_min=phi_min)(phi_op)
        assert np.all(np.isfinite(nl.full()))
        assert (nl - nl.dag()).norm() < 1e-10


# ── end-to-end via epr_numerical_diagonalization ─────────────────────────────

class TestEndToEndGenericPotential:
    def test_generic_cos_matches_builtin_use_full_cos(self):
        """make_nonlinear_potential(cos) passed as non_linear_potential must
        give the same result as use_full_cos=True."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization, make_nonlinear_potential
        freqs, Ljs, phi_zpf = _transmon_params()
        nl = make_nonlinear_potential(np.cos)
        f_generic, chi_generic = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=15, non_linear_potential=nl
        )
        f_builtin, chi_builtin = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=15, use_full_cos=True
        )
        # Both methods compute the same quantity via different floating-point paths;
        # eigendecomposition vs. matrix exponential differ at ~1e-7 level.
        np.testing.assert_allclose(np.real(f_generic), np.real(f_builtin), rtol=1e-6)
        np.testing.assert_allclose(np.real(chi_generic), np.real(chi_builtin), rtol=1e-4)

    def test_fluxonium_generic_returns_finite(self):
        """make_nonlinear_potential works for large phi_zpf (fluxonium regime)."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization, make_nonlinear_potential
        freqs, Ljs, phi_zpf = _fluxonium_params()
        nl = make_nonlinear_potential(np.cos)
        f_ND, chi_ND = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=20, non_linear_potential=nl
        )
        assert np.all(np.isfinite(np.real(f_ND)))
        assert np.all(np.isfinite(np.real(chi_ND)))

    def test_fluxonium_generic_matches_use_full_cos(self):
        """Generic cos potential and use_full_cos=True must agree for fluxonium."""
        from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization, make_nonlinear_potential
        freqs, Ljs, phi_zpf = _fluxonium_params()
        nl = make_nonlinear_potential(np.cos)
        f_gen, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=20, non_linear_potential=nl
        )
        f_full, _ = epr_numerical_diagonalization(
            freqs, Ljs, phi_zpf, fock_trunc=20, use_full_cos=True
        )
        # Both use eigendecomposition — should be machine-precision identical
        np.testing.assert_allclose(np.real(f_gen), np.real(f_full), rtol=1e-8)
