"""
Coverage tests for pure-logic calc modules that don't require an HFSS session.

Modules covered:
  - pyEPR.calcs.quantum   (was 0%)
  - pyEPR.calcs.basic     (was 53%, missing epr_to_zpf / epr_cap_to_nzpf)
  - pyEPR.calcs.transmon  (was 36%, missing dispersiveH / transmon_get_all_params /
                            charge_dispersion_approx)
  - pyEPR.calcs.convert   (was 80%, missing Ic_from_Lj, Lj_from_Ic, ZPF_from_LC,
                            ZPF_from_EPR)
  - core_distributed_analysis pure-logic paths: _parse_listvariations,
    get_nominal_variation_index
"""
import types

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# calcs.quantum — operator matrices
# ---------------------------------------------------------------------------

from pyEPR.calcs.quantum import create, destroy, number, basis


class TestQuantumOperators:

    def test_destroy_shape(self):
        a = destroy(4)
        assert a.shape == (4, 4)

    def test_create_is_destroy_dag(self):
        n = 5
        assert np.allclose(create(n), destroy(n).T)

    def test_destroy_matrix_elements(self):
        # a |1⟩ = |0⟩, a |2⟩ = sqrt(2)|1⟩, a |0⟩ = 0
        a = destroy(4)
        assert np.isclose(a[0, 1], 1.0)
        assert np.isclose(a[1, 2], np.sqrt(2))
        assert np.isclose(a[2, 3], np.sqrt(3))
        assert np.isclose(a[0, 0], 0.0)

    def test_create_matrix_elements(self):
        adag = create(4)
        assert np.isclose(adag[1, 0], 1.0)
        assert np.isclose(adag[2, 1], np.sqrt(2))

    def test_number_operator(self):
        N = number(5)
        assert N.shape == (5, 5)
        assert np.allclose(np.diag(N), [0, 1, 2, 3, 4])
        assert np.allclose(N - np.diag(np.diag(N)), 0)  # off-diagonal zero

    def test_number_from_ladder(self):
        n = 6
        N_constructed = create(n) @ destroy(n)
        assert np.allclose(N_constructed, number(n), atol=1e-12)

    def test_basis_shape(self):
        v = basis(2, 5)
        assert v.shape == (5, 1)

    def test_basis_single_nonzero(self):
        v = basis(3, 6)
        assert np.isclose(v[3, 0], 1.0)
        assert np.isclose(np.sum(np.abs(v)), 1.0)

    def test_basis_orthogonality(self):
        N = 4
        vecs = [basis(i, N) for i in range(N)]
        for i in range(N):
            for j in range(N):
                inner = (vecs[i].T @ vecs[j]).item()
                assert np.isclose(inner, 1.0 if i == j else 0.0)

    def test_commutator_a_adag(self):
        """[a, a†] = I  (bosonic commutation relation)"""
        n = 8
        a = destroy(n)
        adag = create(n)
        comm = a @ adag - adag @ a
        # Only exact for infinite dim; for finite n the last diagonal element is wrong
        assert np.allclose(comm[:-1, :-1], np.eye(n - 1))


# ---------------------------------------------------------------------------
# calcs.basic — epr_to_zpf
# ---------------------------------------------------------------------------

from pyEPR.calcs.basic import CalcsBasic


class TestEprToZpf:

    def _simple_single_mode(self):
        """Single mode, single junction — known analytic result."""
        p   = np.array([[1.0]])   # full participation
        s   = np.array([[1.0]])   # positive sign
        Om  = np.array([[5.0]])   # GHz (diagonal matrix)
        Ej  = np.array([[20.0]])  # GHz (diagonal matrix)
        return p, s, Om, Ej

    def test_shape(self):
        p, s, Om, Ej = self._simple_single_mode()
        result = CalcsBasic.epr_to_zpf(p, s, Om, Ej)
        assert result.shape == (1, 1)

    def test_analytic_value(self):
        """phi_zpf = sign * sqrt(0.5 * Om * p / Ej)"""
        p, s, Om, Ej = self._simple_single_mode()
        expected = np.sqrt(0.5 * 5.0 * 1.0 / 20.0)
        result = CalcsBasic.epr_to_zpf(p, s, Om, Ej)
        assert np.isclose(result[0, 0], expected)

    def test_sign_flip(self):
        p, s, Om, Ej = self._simple_single_mode()
        s_neg = -s
        pos = CalcsBasic.epr_to_zpf(p, s, Om, Ej)
        neg = CalcsBasic.epr_to_zpf(p, s_neg, Om, Ej)
        assert np.isclose(pos, -neg)

    def test_two_mode_shape(self):
        p  = np.array([[0.9, 0.05], [0.05, 0.8]])
        s  = np.ones((2, 2))
        Om = np.diag([5.0, 7.0])
        Ej = np.diag([20.0, 18.0])
        result = CalcsBasic.epr_to_zpf(p, s, Om, Ej)
        assert result.shape == (2, 2)

    def test_epr_cap_to_nzpf_shape(self):
        p  = np.array([[0.5]])
        s  = np.array([[1.0]])
        Om = np.array([[5.0]])
        Ec = np.array([[0.3]])
        result = CalcsBasic.epr_cap_to_nzpf(p, s, Om, Ec)
        assert result.shape == (1, 1)


# ---------------------------------------------------------------------------
# calcs.transmon — dispersive Hamiltonian and parameter helpers
# ---------------------------------------------------------------------------

from pyEPR.calcs.transmon import CalcsTransmon
from pyEPR.calcs.convert import Convert


class TestCalcsTransmon:

    def test_dispersive_pt_single_mode(self):
        """Single mode, single junction — chi should be negative (self-Kerr)."""
        Pmj = np.array([[0.9]])
        Ωm  = np.diag([5.0])   # GHz
        Ej  = np.diag([20.0])  # GHz
        f_O1, chi_O1 = CalcsTransmon.dispersiveH_params_PT_O1(Pmj, Ωm, Ej)
        assert f_O1.shape == (1,)
        assert chi_O1.shape == (1, 1)
        # chi_O1[0,0] is the self-Kerr (anharmonicity) — positive in this convention
        assert chi_O1[0, 0] != 0

    def test_dispersive_pt_two_modes(self):
        Pmj = np.array([[0.9, 0.05], [0.05, 0.8]])
        Ωm  = np.diag([5.0, 7.0])
        Ej  = np.diag([20.0, 18.0])
        f_O1, chi_O1 = CalcsTransmon.dispersiveH_params_PT_O1(Pmj, Ωm, Ej)
        assert f_O1.shape == (2,)
        assert chi_O1.shape == (2, 2)

    def test_dispersive_pt_dimension_mismatch(self):
        """Shape mismatch should raise AssertionError."""
        Pmj = np.array([[0.9, 0.05]])   # 1x2 — inconsistent with 1x1 Ej
        Ωm  = np.diag([5.0])
        Ej  = np.diag([20.0])
        with pytest.raises(AssertionError):
            CalcsTransmon.dispersiveH_params_PT_O1(Pmj, Ωm, Ej)

    def test_transmon_get_all_params_keys(self):
        params = CalcsTransmon.transmon_get_all_params(Ej_MHz=20_000, Ec_MHz=300)
        expected_keys = {"Ej_MHz", "Ec_MHz", "Lj_H", "Cs_F", "Lj_nH", "Cs_fF",
                         "Phi_ZPF", "Q_ZPF", "phi_ZPF", "n_ZPF", "Omega_MHz",
                         "f_MHz", "Z_Ohms"}
        assert expected_keys == set(params.keys())

    def test_transmon_get_all_params_values_positive(self):
        params = CalcsTransmon.transmon_get_all_params(Ej_MHz=20_000, Ec_MHz=300)
        for key, val in params.items():
            assert val > 0, f"{key} should be positive, got {val}"

    def test_transmon_get_all_params_freq_reasonable(self):
        """Plasma frequency sqrt(8 Ej Ec) / h for Ej=20GHz, Ec=300MHz ≈ 6.9 GHz"""
        params = CalcsTransmon.transmon_get_all_params(Ej_MHz=20_000, Ec_MHz=300)
        # Note: f_MHz key returns in GHz despite the label name
        f_GHz = params["f_MHz"]
        assert 5.0 < f_GHz < 10.0

    def test_charge_dispersion_approx_m0(self):
        """m=0 should give a positive value for a transmon (large Ej/Ec)."""
        eps = CalcsTransmon.charge_dispersion_approx(m=0, Ec=300, Ej=20_000)
        assert eps > 0

    def test_charge_dispersion_approx_sign_alternates(self):
        """Koch eq 2.5: sign is (-1)^m."""
        Ec, Ej = 300, 20_000
        eps0 = CalcsTransmon.charge_dispersion_approx(0, Ec, Ej)
        eps1 = CalcsTransmon.charge_dispersion_approx(1, Ec, Ej)
        assert eps0 > 0
        assert eps1 < 0

    def test_charge_dispersion_decreases_with_ej(self):
        """Larger Ej/Ec → smaller charge dispersion (deeper in transmon regime)."""
        Ec = 300
        eps_small = abs(CalcsTransmon.charge_dispersion_approx(0, Ec, Ej=5_000))
        eps_large = abs(CalcsTransmon.charge_dispersion_approx(0, Ec, Ej=50_000))
        assert eps_small > eps_large


# ---------------------------------------------------------------------------
# calcs.convert — missing paths
# ---------------------------------------------------------------------------


class TestConvertMissingPaths:

    def test_ic_from_lj_roundtrip(self):
        """Lj → Ic → Lj should round-trip."""
        Lj_nH = 12.0
        Ic_nA = Convert.Ic_from_Lj(Lj_nH, units_in="nH", units_out="nA")
        Lj_back = Convert.Lj_from_Ic(Ic_nA, units_in="nA", units_out="nH")
        assert np.isclose(Lj_back, Lj_nH, rtol=1e-6)

    def test_ic_from_lj_positive(self):
        assert Convert.Ic_from_Lj(10.0, "nH", "nA") > 0

    def test_lj_from_ic_positive(self):
        assert Convert.Lj_from_Ic(100.0, "nA", "nH") > 0

    def test_zpf_from_lc_shape(self):
        Phi_ZPF, Q_ZPF = Convert.ZPF_from_LC(10e-9, 100e-15)
        assert np.isscalar(Phi_ZPF) or Phi_ZPF.ndim == 0
        assert np.isscalar(Q_ZPF) or Q_ZPF.ndim == 0

    def test_zpf_from_lc_positive(self):
        Phi, Q = Convert.ZPF_from_LC(10e-9, 100e-15)
        assert Phi > 0
        assert Q > 0

    def test_zpf_from_lc_heisenberg(self):
        """Phi_ZPF * Q_ZPF = hbar/2 (uncertainty product at minimum)."""
        import scipy.constants as const
        Phi, Q = Convert.ZPF_from_LC(10e-9, 100e-15)
        assert np.isclose(Phi * Q, const.hbar / 2, rtol=1e-6)

    def test_zpf_from_epr_shape(self):
        freqs  = np.array([5.0, 7.0])       # GHz
        epr    = np.array([[0.9, 0.02],
                           [0.02, 0.85]])
        signs  = np.ones((2, 2))
        Ljs    = np.array([10e-9, 12e-9])   # H
        zpfs, mats = Convert.ZPF_from_EPR(freqs, epr, signs, Ljs, Lj_units_in="H")
        assert zpfs.shape == (2, 2)
        Ωd, Ej, epr_out, signs_out = mats
        assert Ωd.shape == (2, 2)
        assert Ej.shape == (2, 2)

    def test_zpf_from_epr_values_reasonable(self):
        """phi_zpf for a transmon (p≈1, f=5 GHz, Lj=10nH) should be ~0.4."""
        freqs = np.array([5.0])
        epr   = np.array([[1.0]])
        signs = np.array([[1.0]])
        Ljs   = np.array([10e-9])
        zpfs, _ = Convert.ZPF_from_EPR(freqs, epr, signs, Ljs, Lj_units_in="H")
        assert 0.1 < zpfs[0, 0] < 1.0


# ---------------------------------------------------------------------------
# core_distributed_analysis — pure-logic paths (no HFSS)
# ---------------------------------------------------------------------------

from pyEPR.core_distributed_analysis import DistributedAnalysis


def _make_da_stub(list_variations=None, nominal=None, variations=None):
    list_variations = list_variations or ("Cj='2fF' Lj='12nH'", "Cj='2fF' Lj='13nH'")
    obj = types.SimpleNamespace()
    obj._list_variations = list_variations
    obj._nominal_variation = nominal or list_variations[0]
    obj.variations = variations or [str(i) for i in range(len(list_variations))]
    # bind unbound methods
    obj._parse_listvariations = DistributedAnalysis._parse_listvariations.__get__(obj)
    obj._variation_index      = DistributedAnalysis._variation_index.__get__(obj)
    obj.get_nominal_variation_index = DistributedAnalysis.get_nominal_variation_index.__get__(obj)
    return obj


class TestParseListvariations:

    def setup_method(self):
        self.stub = _make_da_stub()

    def test_simple(self):
        result = self.stub._parse_listvariations("Cj='2fF' Lj='13.5nH'")
        assert result == ["Cj:=", "2fF", "Lj:=", "13.5nH"]

    def test_empty_string(self):
        result = self.stub._parse_listvariations("")
        # empty string → [""] after split
        assert result == [""]

    def test_single_variable(self):
        result = self.stub._parse_listvariations("Lj='12nH'")
        assert "Lj:=" in result
        assert "12nH" in result

    def test_three_variables(self):
        result = self.stub._parse_listvariations("A='1mm' B='2nH' C='3fF'")
        assert len(result) == 6
        assert result[0] == "A:="
        assert result[1] == "1mm"


class TestGetNominalVariationIndex:

    def test_nominal_is_first(self):
        stub = _make_da_stub(
            list_variations=("Lj='12nH'", "Lj='13nH'", "Lj='14nH'"),
            nominal="Lj='12nH'",
        )
        assert stub.get_nominal_variation_index() == "0"

    def test_nominal_is_second(self):
        stub = _make_da_stub(
            list_variations=("Lj='12nH'", "Lj='13nH'", "Lj='14nH'"),
            nominal="Lj='13nH'",
        )
        assert stub.get_nominal_variation_index() == "1"

    def test_nominal_not_in_list_returns_zero(self):
        """If the nominal variation isn't solved, fall back to '0'."""
        stub = _make_da_stub(
            list_variations=("Lj='12nH'", "Lj='13nH'"),
            nominal="Lj='99nH'",  # not in solved list
        )
        assert stub.get_nominal_variation_index() == "0"
