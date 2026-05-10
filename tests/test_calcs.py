"""
Tests for pure-computation functions in pyEPR.calcs.
No HFSS or Ansys connection required — all inputs are synthetic.
"""
import numpy as np
import pytest

# ---------------------------------------------------------------------------
# epr_numerical_diagonalization
# ---------------------------------------------------------------------------

def test_epr_diag_single_mode_shapes():
    """Output shapes are correct for a 1-mode 1-junction system."""
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
    freqs   = np.array([5.0])
    Ljs     = np.array([10e-9])
    phi_zpf = np.array([[0.15]])
    f_ND, chi_ND = epr_numerical_diagonalization(freqs, Ljs, phi_zpf,
                                                 cos_trunc=8, fock_trunc=9)
    assert f_ND.shape == (1,), f"Expected (1,), got {f_ND.shape}"
    assert chi_ND.shape == (1, 1), f"Expected (1,1), got {chi_ND.shape}"


def test_epr_diag_single_mode_physical_values():
    """
    1-mode transmon: dressed frequency is close to linearized frequency,
    and self-Kerr (anharmonicity) is positive (sign convention: down shift = positive).
    Reference values captured from a known-good run.
    """
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
    freqs   = np.array([5.0])
    Ljs     = np.array([10e-9])
    phi_zpf = np.array([[0.15]])
    f_ND, chi_ND = epr_numerical_diagonalization(freqs, Ljs, phi_zpf,
                                                 cos_trunc=8, fock_trunc=9)
    # Dressed frequency should be within 0.5% of the bare 5 GHz (in Hz)
    assert abs(f_ND[0].real - 5e9) / 5e9 < 0.005
    # Self-Kerr (anharmonicity) should be positive and physically reasonable (~1-100 MHz)
    kerr = chi_ND[0, 0].real
    assert kerr > 0, "Self-Kerr should be positive (sign convention)"
    assert 0.1 < kerr < 500, f"Self-Kerr {kerr:.3f} MHz outside expected range [0.1, 500]"


def test_epr_diag_two_mode_shapes_and_symmetry():
    """
    2-mode system (qubit + resonator): output shapes and dispersive matrix symmetry.
    chi_ND should be symmetric (chi_01 == chi_10).
    """
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
    freqs   = np.array([5.0, 6.5])
    Ljs     = np.array([10e-9])
    phi_zpf = np.array([[0.15], [0.01]])
    f_ND, chi_ND = epr_numerical_diagonalization(freqs, Ljs, phi_zpf,
                                                 cos_trunc=8, fock_trunc=9)
    assert f_ND.shape == (2,)
    assert chi_ND.shape == (2, 2)
    # Dispersive coupling is symmetric
    assert np.isclose(chi_ND[0, 1].real, chi_ND[1, 0].real, rtol=1e-6)
    # Qubit anharmonicity (chi_00) >> resonator anharmonicity (chi_11)
    assert chi_ND[0, 0].real > chi_ND[1, 1].real


def test_epr_diag_two_mode_regression():
    """Regression: known-good values for a 2-mode qubit+resonator system."""
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
    freqs   = np.array([5.0, 6.5])
    Ljs     = np.array([10e-9])
    phi_zpf = np.array([[0.15], [0.01]])
    f_ND, chi_ND = epr_numerical_diagonalization(freqs, Ljs, phi_zpf,
                                                 cos_trunc=8, fock_trunc=9)
    # Frequencies: qubit ~5 GHz, resonator ~6.5 GHz (in Hz)
    assert np.isclose(f_ND[0].real, 4.99586e9, rtol=1e-4)
    assert np.isclose(f_ND[1].real, 6.49998e9, rtol=1e-4)
    # Kerr anharmonicity: qubit ~4.1 MHz
    assert np.isclose(chi_ND[0, 0].real, 4.106, rtol=1e-2)
    # Cross-Kerr qubit-resonator ~0.036 MHz
    assert np.isclose(chi_ND[0, 1].real, 0.0362, rtol=5e-2)


def test_epr_diag_bad_units_raises():
    """Passing frequencies in Hz instead of GHz should raise AssertionError."""
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
    with pytest.raises(AssertionError):
        epr_numerical_diagonalization(np.array([5e9]), np.array([10e-9]),
                                      np.array([[0.15]]))


def test_epr_diag_bad_lj_units_raises():
    """Passing Lj in nH instead of Henries should raise AssertionError."""
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
    with pytest.raises(AssertionError):
        epr_numerical_diagonalization(np.array([5.0]), np.array([10.0]),
                                      np.array([[0.15]]))


# ---------------------------------------------------------------------------
# Convert: unit-conversion roundtrips and known physical values
# ---------------------------------------------------------------------------

def test_convert_lj_ej_roundtrip():
    """Lj -> Ej -> Lj roundtrip should be lossless."""
    from pyEPR.calcs.convert import Convert
    Lj = 10e-9  # 10 nH
    Ej = Convert.Ej_from_Lj(Lj, units_in='H', units_out='GHz')
    Lj_back = Convert.Lj_from_Ej(Ej, units_in='GHz', units_out='H')
    assert np.isclose(Lj, Lj_back, rtol=1e-9)


def test_convert_ec_cs_roundtrip():
    """Ec -> Cs -> Ec roundtrip should be lossless."""
    from pyEPR.calcs.convert import Convert
    Ec = 0.2  # GHz
    Cs = Convert.Cs_from_Ec(Ec, units_in='GHz', units_out='fF')
    Ec_back = Convert.Ec_from_Cs(Cs, units_in='fF', units_out='GHz')
    assert np.isclose(Ec, Ec_back, rtol=1e-9)


def test_convert_ej_physical_range():
    """10 nH junction should have Ej in the physically expected range (~16 GHz)."""
    from pyEPR.calcs.convert import Convert
    Ej_GHz = Convert.Ej_from_Lj(10e-9, units_in='H', units_out='GHz')
    assert 10 < Ej_GHz < 25, f"Ej {Ej_GHz:.1f} GHz outside expected range for 10 nH"


def test_convert_omega_from_lc():
    """LC resonator frequency formula: 1/sqrt(LC) should match expected GHz."""
    from pyEPR.calcs.convert import Convert
    import numpy as np
    # L=1 nH, C=1 fF -> f = 1/(2*pi*sqrt(LC)) ~ 159 GHz
    L, C = 1e-9, 1e-15
    omega = Convert.Omega_from_LC(L, C)
    f_GHz = omega / (2 * np.pi) / 1e9
    assert np.isclose(f_GHz, 159.15, rtol=1e-3), f"Got {f_GHz:.2f} GHz, expected ~159.15"


# ---------------------------------------------------------------------------
# pyEPR package-level import and version
# ---------------------------------------------------------------------------

def test_pyepr_imports():
    """Core pyEPR modules should import without error."""
    import pyEPR as epr
    assert hasattr(epr, '__version__')
    assert epr.__version__  # non-empty string


def test_pyepr_version_format():
    """Version string should follow semver-like format X.Y.Z."""
    import pyEPR as epr
    parts = epr.__version__.split('.')
    assert len(parts) >= 2
    assert all(p.isdigit() for p in parts)
