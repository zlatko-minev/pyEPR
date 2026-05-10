"""
Tests for pyEPR.solution_types — alias sets, predicates, and normalize().

These tests act as contracts: they pin every known GetSolutionType() alias
to its canonical family and assert that the sets are disjoint, so a typo or
duplicate alias is caught immediately.

Designed to run without Ansys, COM, Qt, or Windows-specific dependencies.
"""
import pytest

from pyEPR.solution_types import (
    DRIVEN_MODAL,
    DRIVEN_MODAL_NAMES,
    DRIVEN_TERMINAL,
    DRIVEN_TERMINAL_NAMES,
    EIGENMODE,
    EIGENMODE_NAMES,
    Q3D,
    Q3D_NAMES,
    canonical_kind,
    is_drivenmodal,
    is_driventerminal,
    is_eigenmode,
    is_q3d,
    normalize,
)


# ── Canonical constants ────────────────────────────────────────────────────────

def test_canonical_constants():
    assert EIGENMODE == "Eigenmode"
    assert DRIVEN_MODAL == "DrivenModal"
    assert DRIVEN_TERMINAL == "DrivenTerminal"
    assert Q3D == "Q3D"


# ── Alias membership ───────────────────────────────────────────────────────────

class TestEigenmodeMembership:
    def test_eigenmode(self):
        assert "Eigenmode" in EIGENMODE_NAMES

    def test_canonical_is_member(self):
        assert EIGENMODE in EIGENMODE_NAMES


class TestDrivenModalMembership:
    def test_legacy(self):
        assert "DrivenModal" in DRIVEN_MODAL_NAMES

    def test_aedt_2021_2(self):
        assert "HFSS Modal Network" in DRIVEN_MODAL_NAMES

    def test_hybrid(self):
        assert "HFSS Hybrid Modal Network" in DRIVEN_MODAL_NAMES

    def test_canonical_is_member(self):
        assert DRIVEN_MODAL in DRIVEN_MODAL_NAMES


class TestDrivenTerminalMembership:
    def test_legacy(self):
        assert "DrivenTerminal" in DRIVEN_TERMINAL_NAMES

    def test_aedt_2021_2(self):
        assert "HFSS Terminal Network" in DRIVEN_TERMINAL_NAMES

    def test_hybrid(self):
        assert "HFSS Hybrid Terminal Network" in DRIVEN_TERMINAL_NAMES

    def test_canonical_is_member(self):
        assert DRIVEN_TERMINAL in DRIVEN_TERMINAL_NAMES


class TestQ3DMembership:
    def test_q3d(self):
        assert "Q3D" in Q3D_NAMES

    def test_canonical_is_member(self):
        assert Q3D in Q3D_NAMES


# ── Disjointness ───────────────────────────────────────────────────────────────

def test_sets_are_disjoint():
    """No alias may belong to two families at once."""
    all_sets = [EIGENMODE_NAMES, DRIVEN_MODAL_NAMES, DRIVEN_TERMINAL_NAMES, Q3D_NAMES]
    for i, a in enumerate(all_sets):
        for j, b in enumerate(all_sets):
            if i != j:
                assert a.isdisjoint(b), (
                    f"Sets {i} and {j} share: {a & b}"
                )


# ── Predicates ─────────────────────────────────────────────────────────────────

class TestPredicates:
    @pytest.mark.parametrize("s", list(EIGENMODE_NAMES))
    def test_is_eigenmode(self, s):
        assert is_eigenmode(s)
        assert not is_drivenmodal(s)
        assert not is_driventerminal(s)
        assert not is_q3d(s)

    @pytest.mark.parametrize("s", list(DRIVEN_MODAL_NAMES))
    def test_is_drivenmodal(self, s):
        assert is_drivenmodal(s)
        assert not is_eigenmode(s)
        assert not is_driventerminal(s)
        assert not is_q3d(s)

    @pytest.mark.parametrize("s", list(DRIVEN_TERMINAL_NAMES))
    def test_is_driventerminal(self, s):
        assert is_driventerminal(s)
        assert not is_eigenmode(s)
        assert not is_drivenmodal(s)
        assert not is_q3d(s)

    @pytest.mark.parametrize("s", list(Q3D_NAMES))
    def test_is_q3d(self, s):
        assert is_q3d(s)
        assert not is_eigenmode(s)
        assert not is_drivenmodal(s)
        assert not is_driventerminal(s)


# ── canonical_kind ─────────────────────────────────────────────────────────────

class TestCanonicalKind:
    @pytest.mark.parametrize("s", list(EIGENMODE_NAMES))
    def test_eigenmode_kind(self, s):
        assert canonical_kind(s) == "eigenmode"

    @pytest.mark.parametrize("s", list(DRIVEN_MODAL_NAMES))
    def test_drivenmodal_kind(self, s):
        assert canonical_kind(s) == "drivenmodal"

    @pytest.mark.parametrize("s", list(DRIVEN_TERMINAL_NAMES))
    def test_driventerminal_kind(self, s):
        assert canonical_kind(s) == "driventerminal"

    @pytest.mark.parametrize("s", list(Q3D_NAMES))
    def test_q3d_kind(self, s):
        assert canonical_kind(s) == "q3d"

    def test_unknown_returns_none(self):
        assert canonical_kind("SBR+") is None
        assert canonical_kind("Transient") is None
        assert canonical_kind("") is None

    def test_case_sensitive(self):
        assert canonical_kind("drivenmodal") is None
        assert canonical_kind("DRIVENMODAL") is None
        assert canonical_kind("eigenmode") is None


# ── normalize ──────────────────────────────────────────────────────────────────

class TestNormalize:
    def test_legacy_drivenmodal(self):
        assert normalize("DrivenModal") == "DrivenModal"

    def test_aedt_modal_network(self):
        assert normalize("HFSS Modal Network") == "DrivenModal"

    def test_hybrid_modal(self):
        assert normalize("HFSS Hybrid Modal Network") == "DrivenModal"

    def test_legacy_driventerminal(self):
        assert normalize("DrivenTerminal") == "DrivenTerminal"

    def test_aedt_terminal_network(self):
        assert normalize("HFSS Terminal Network") == "DrivenTerminal"

    def test_hybrid_terminal(self):
        assert normalize("HFSS Hybrid Terminal Network") == "DrivenTerminal"

    def test_eigenmode(self):
        assert normalize("Eigenmode") == "Eigenmode"

    def test_q3d_passthrough(self):
        assert normalize("Q3D") == "Q3D"

    def test_unknown_passthrough(self):
        assert normalize("SBR+") == "SBR+"
        assert normalize("Transient") == "Transient"

    def test_normalize_output_is_canonical(self):
        """Every known alias normalises to its canonical constant."""
        for s in DRIVEN_MODAL_NAMES:
            assert normalize(s) == DRIVEN_MODAL, f"normalize({s!r}) != {DRIVEN_MODAL!r}"
        for s in DRIVEN_TERMINAL_NAMES:
            assert normalize(s) == DRIVEN_TERMINAL, f"normalize({s!r}) != {DRIVEN_TERMINAL!r}"
        for s in EIGENMODE_NAMES:
            assert normalize(s) == EIGENMODE, f"normalize({s!r}) != {EIGENMODE!r}"
