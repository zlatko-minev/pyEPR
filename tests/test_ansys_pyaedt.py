"""
Tests for pyEPR.ansys_pyaedt — the PyAEDT (gRPC) HFSS backend.

The backend talks to HFSS through Ansys's PyAEDT API over gRPC, but its
*physics* (the junction participation formula), its unit parsing, and its
session-targeting logic are pure Python and are tested here without Ansys,
COM, PyAEDT, or Windows.

Only the single ``@pytest.mark.hfss`` test needs a live AEDT session; it is
skipped in CI and is also skipped locally unless the env var
``PYEPR_PYAEDT_TEST_PROJECT`` points at a solved demo ``.aedt`` project.

Designed to run on Linux/macOS/Windows with no optional dependencies.
"""
import os
from types import SimpleNamespace

import pytest

from pyEPR.ansys_pyaedt import (
    PyaedtDistributedAnalysis,
    compute_p_mj,
    _parse_henries,
    _owning_session_from_lock,
    _require_pyaedt,
)


# ── compute_p_mj — the junction participation formula ───────────────────────

# Golden values from the demo transmon, validated digit-for-digit against the
# COM path (qubit eigenmode 4.815509 GHz, Lj 12.9201 nH).
GOLDEN_V = -2.200983e-05
GOLDEN_FREQ_HZ = 4.815509e9
GOLDEN_LJ_H = 12.9201e-9
GOLDEN_UE_RAW = 4.198685e-23


def test_compute_p_mj_matches_golden_transmon():
    p, sign = compute_p_mj(
        V_peak=GOLDEN_V,
        freq_hz=GOLDEN_FREQ_HZ,
        Lj_henries=GOLDEN_LJ_H,
        U_E_raw=GOLDEN_UE_RAW,
    )
    assert p == pytest.approx(0.9755, abs=2e-3)
    assert sign == -1


def test_compute_p_mj_sign_follows_voltage():
    common = dict(freq_hz=GOLDEN_FREQ_HZ, Lj_henries=GOLDEN_LJ_H, U_E_raw=GOLDEN_UE_RAW)
    p_pos, s_pos = compute_p_mj(V_peak=abs(GOLDEN_V), **common)
    p_neg, s_neg = compute_p_mj(V_peak=-abs(GOLDEN_V), **common)
    assert s_pos == +1 and s_neg == -1
    # participation is sign-independent (V enters squared)
    assert p_pos == pytest.approx(p_neg)


def test_compute_p_mj_capacitance_lowers_participation():
    base, _ = compute_p_mj(
        V_peak=GOLDEN_V, freq_hz=GOLDEN_FREQ_HZ,
        Lj_henries=GOLDEN_LJ_H, U_E_raw=GOLDEN_UE_RAW,
    )
    with_cj, _ = compute_p_mj(
        V_peak=GOLDEN_V, freq_hz=GOLDEN_FREQ_HZ,
        Lj_henries=GOLDEN_LJ_H, U_E_raw=GOLDEN_UE_RAW,
        Cj_farads=1e-15,
    )
    # adding junction capacitance grows the normalization → smaller p_mj
    assert with_cj < base


def test_compute_p_mj_rejects_bad_inputs():
    with pytest.raises(ValueError):
        compute_p_mj(V_peak=1.0, freq_hz=0.0, Lj_henries=GOLDEN_LJ_H, U_E_raw=1e-23)
    with pytest.raises(ValueError):
        compute_p_mj(V_peak=1.0, freq_hz=GOLDEN_FREQ_HZ, Lj_henries=0.0, U_E_raw=1e-23)
    with pytest.raises(ValueError):
        # non-positive electric energy → undefined normalization
        compute_p_mj(V_peak=1.0, freq_hz=GOLDEN_FREQ_HZ, Lj_henries=GOLDEN_LJ_H, U_E_raw=0.0)


# ── _parse_henries — inductance string → SI Henries ─────────────────────────

@pytest.mark.parametrize("expr,expected", [
    ("12.9201nH", 12.9201e-9),
    ("5uH", 5e-6),
    ("2mH", 2e-3),
    ("1.5H", 1.5),
    ("12.9201NH", 12.9201e-9),   # case-insensitive
    ("  3nH ", 3e-9),            # surrounding whitespace
    ("4e-9", 4e-9),             # unit-less / already SI
    ("4.2", 4.2),
])
def test_parse_henries(expr, expected):
    assert _parse_henries(expr) == pytest.approx(expected)


# ── _owning_session_from_lock — .aedt.lock targeting ────────────────────────

def test_owning_session_reads_pid_from_lock(tmp_path):
    project = tmp_path / "Demo.aedt"
    lock = tmp_path / "Demo.aedt.lock"
    lock.write_text("SomeHeader=1\nDesktopProcessID=4242\nOther=x\n")
    grpc = {4242: 50051, 9999: 50060}
    assert _owning_session_from_lock(str(project), grpc) == (4242, 50051)


def test_owning_session_none_when_pid_not_running(tmp_path):
    project = tmp_path / "Demo.aedt"
    (tmp_path / "Demo.aedt.lock").write_text("DesktopProcessID=4242\n")
    assert _owning_session_from_lock(str(project), {9999: 50051}) is None


def test_owning_session_none_when_no_lock(tmp_path):
    project = tmp_path / "Demo.aedt"
    assert _owning_session_from_lock(str(project), {4242: 50051}) is None


# ── lazy PyAEDT import — backend must import without pyaedt ──────────────────

def test_backend_imports_and_constructs_without_pyaedt():
    """Importing the backend and constructing the object must not need PyAEDT."""
    pinfo = SimpleNamespace(
        junctions={}, project_path=None, project_name=None,
        design_name=None, setup_name=None,
    )
    obj = PyaedtDistributedAnalysis(pinfo, aedt_version="2026.1")
    assert obj.junction_names == []
    assert obj.PJ is None  # not analyzed yet


def test_require_pyaedt_raises_clear_error_when_missing(monkeypatch):
    """When pyaedt is absent, the error must name PyAEDT and how to fix it."""
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name.startswith("ansys.aedt"):
            raise ImportError("simulated missing PyAEDT")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="PyAEDT is required"):
        _require_pyaedt()


# ── live extraction (needs AEDT + the solved demo project) ──────────────────

@pytest.mark.hfss
def test_live_pmj_matches_golden():
    """End-to-end p_mj from a live AEDT session equals the golden 0.9755.

    Skipped unless ``PYEPR_PYAEDT_TEST_PROJECT`` points at the solved demo
    ``.aedt``; requires AEDT open with that project. Mirrors the demo notebook.
    """
    project = os.environ.get("PYEPR_PYAEDT_TEST_PROJECT")
    if not project:
        pytest.skip("set PYEPR_PYAEDT_TEST_PROJECT to the solved demo .aedt to run")

    import pyEPR as epr

    pinfo = epr.ProjectInfo(
        project_path=os.path.dirname(project),
        project_name=os.path.splitext(os.path.basename(project))[0],
        design_name="EPR_Sample_Demo",
        setup_name="EPR_Scan",
        do_connect=False,
    )
    pinfo.junctions["j1"] = {"Lj_variable": "Lj_Transmon", "line": "Junction_line"}

    eprd = PyaedtDistributedAnalysis(pinfo, aedt_version="2026.1")
    try:
        eprd.do_EPR_analysis()
    finally:
        eprd.disconnect()

    # qubit mode is the one with participation near unity
    assert eprd.PJ.max() == pytest.approx(0.9755, abs=2e-3)
