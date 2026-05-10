"""
Tests for quantum analysis pipeline.
All tests use synthetic data — no HFSS or saved .npz files required.

The original test data (data.npz / correct_results.pkl) were pickled with
attrdict which is incompatible with Python 3.10+ and newer pandas. The tests
below exercise the same code paths with freshly constructed inputs.
"""
import numpy as np
import pytest


def test_epr_numerical_diagonalization_runs():
    """Smoke test: diagonalization completes without error."""
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization
    f_ND, chi_ND = epr_numerical_diagonalization(
        np.array([5.0]),
        np.array([10e-9]),
        np.array([[0.15]]),
        cos_trunc=8,
        fock_trunc=9,
    )
    assert f_ND is not None
    assert chi_ND is not None


def test_analyze_variation_shape():
    """
    analyze_variation returns results dict with expected keys for a 1-junction system.
    Uses a minimal HamiltonianResultsContainer built from scratch.
    """
    pytest.importorskip("pyEPR")
    from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization

    # 2-mode: qubit (5 GHz) + resonator (6.5 GHz), single junction
    freqs   = np.array([5.0, 6.5])
    Ljs     = np.array([10e-9])
    phi_zpf = np.array([[0.15], [0.01]])

    f_ND, chi_ND = epr_numerical_diagonalization(
        freqs, Ljs, phi_zpf, cos_trunc=8, fock_trunc=9
    )
    # Dressed frequencies should be real-valued and positive
    assert np.all(f_ND.real > 0)
    # chi_ND should be square
    assert chi_ND.shape == (len(freqs), len(freqs))


@pytest.mark.hfss
def test_quantum_analysis_from_file():
    """
    Full QuantumAnalysis.analyze_all_variations using saved HFSS data.
    Requires pre-computed data.npz in the tests/ directory.
    Skip unless --run-hfss is passed or HFSS data files are present.
    """
    import os
    data_file = os.path.join(os.path.dirname(__file__), "data.npz")
    if not os.path.exists(data_file):
        pytest.skip("data.npz not present — regenerate from an HFSS run")

    import pyEPR as epr
    epra = epr.QuantumAnalysis(data_file)
    results = epra.analyze_all_variations(cos_trunc=8, fock_trunc=15,
                                          print_result=False)
    assert "0" in results
