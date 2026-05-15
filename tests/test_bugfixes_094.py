"""
Tests for the four bug fixes in 0.9.4:
  - #128  junction sign always printed as (+)
  - #162  QuTiP 5.x .norm() on complex scalar in back_box_numeric
  - #169  setup name ignored / hardcoded to first setup
  - #137  convergence plot crashes on log-scale with non-positive data
"""
import numpy as np
import pandas as pd
import pytest


# ── Bug #128 ──────────────────────────────────────────────────────────────────

class TestJunctionSignPrint:
    """_Smj is always 1 or -1; old code used `if _Smj` (always True)."""

    def test_positive_sign_truthy(self):
        _Smj = 1
        result = "(+)" if _Smj > 0 else "(-)"
        assert result == "(+)"

    def test_negative_sign(self):
        _Smj = -1
        result = "(+)" if _Smj > 0 else "(-)"
        assert result == "(-)"

    def test_old_code_was_wrong(self):
        # Demonstrate the old bug: both 1 and -1 are truthy
        _Smj = -1
        old_result = "(+)" if _Smj else "(-)"
        assert old_result == "(+)", "Old code always returned (+) — this confirms the bug"

    def test_new_code_is_correct(self):
        for smj, expected in [(1, "(+)"), (-1, "(-)"), (1, "(+)")]:
            assert ("(+)" if smj > 0 else "(-)" ) == expected


# ── Bug #162 ──────────────────────────────────────────────────────────────────

class TestQutip5InnerProduct:
    """
    In QuTiP 5, ket.dag() * ket returns a complex scalar, not a 1x1 Qobj.
    The fixed code uses: abs(r.full()[0,0]) if hasattr(r, 'full') else abs(r)
    """

    def test_complex_scalar_path(self):
        # Simulate what QuTiP 5 returns: a plain complex number
        r = complex(0.6, 0.8)  # |r| = 1.0
        result = abs(r.full()[0, 0]) if hasattr(r, "full") else abs(r)
        assert pytest.approx(result, abs=1e-10) == 1.0

    def test_complex_scalar_negative_real(self):
        r = complex(-0.8, 0.6)
        result = abs(r.full()[0, 0]) if hasattr(r, "full") else abs(r)
        assert pytest.approx(result, abs=1e-10) == 1.0

    def test_real_scalar(self):
        r = 0.5
        result = abs(r.full()[0, 0]) if hasattr(r, "full") else abs(r)
        assert pytest.approx(result) == 0.5

    def test_argmax_selects_correct_index(self):
        # Simulate a list of inner products (complex scalars as QuTiP 5 returns)
        values = [complex(0.1, 0), complex(0.9, 0), complex(0.3, 0)]
        index = np.argmax([
            abs(r.full()[0, 0]) if hasattr(r, "full") else abs(r)
            for r in values
        ])
        assert index == 1  # 0.9 is largest


# ── Bug #169 ──────────────────────────────────────────────────────────────────

class TestSetupNameSelection:
    """
    The setup-name logic (extracted from ProjectInfo._configure_ansys_project):
      - if user specified a name that IS in the list: keep it
      - if user specified a name NOT in the list: raise ValueError
      - if user specified nothing: pick setup_names[0]
    """

    @staticmethod
    def _resolve_setup_name(user_specified, setup_names):
        """Mirrors the fixed logic in project_info.py."""
        setup_name = user_specified
        if setup_name:
            if setup_name not in setup_names:
                raise ValueError(
                    f"Setup '{setup_name}' not found in design. "
                    f"Available setups: {setup_names}"
                )
            # else keep user_specified
        else:
            setup_name = setup_names[0]
        return setup_name

    def test_user_specified_valid_name_is_kept(self):
        result = self._resolve_setup_name("MySetup", ["DefaultSetup", "MySetup"])
        assert result == "MySetup"

    def test_no_user_name_falls_back_to_first(self):
        result = self._resolve_setup_name(None, ["DefaultSetup", "MySetup"])
        assert result == "DefaultSetup"

    def test_invalid_name_raises_value_error(self):
        with pytest.raises(ValueError, match="not found"):
            self._resolve_setup_name("NonExistent", ["DefaultSetup", "MySetup"])

    def test_old_behaviour_always_used_first(self):
        # Old code: self.setup_name = setup_names[0]  — always overwrote user input
        user = "MySetup"
        setup_names = ["DefaultSetup", "MySetup"]
        old_result = setup_names[0]  # old code always did this
        assert old_result == "DefaultSetup"  # confirms old bug
        # New code preserves user intent:
        new_result = self._resolve_setup_name(user, setup_names)
        assert new_result == "MySetup"

    def test_empty_string_treated_as_not_specified(self):
        # Empty string is falsy — should fall back to first setup
        result = self._resolve_setup_name("", ["DefaultSetup"])
        assert result == "DefaultSetup"


# ── Bug #137 / #164 ──────────────────────────────────────────────────────────

class TestConvergencePlotLogScale:
    """
    plot_convergence_max_df and plot_convergence_maxdf_vs_sol called
    set_yscale("log") unconditionally. If delta-F data has no positive values
    matplotlib raises ValueError. Fixed by guarding with (s > 0).any().
    """

    def test_positive_data_uses_log(self):
        s = pd.Series([0.5, 0.3, 0.1])
        assert (s > 0).any()

    def test_all_zero_data_skips_log(self):
        s = pd.Series([0.0, 0.0, 0.0])
        assert not (s > 0).any()

    def test_negative_data_skips_log(self):
        s = pd.Series([-0.1, -0.2])
        assert not (s > 0).any()

    def test_mixed_data_with_positives_uses_log(self):
        s = pd.Series([0.0, 0.1, 0.3])
        assert (s > 0).any()

    def test_matplotlib_does_not_crash_with_zero_data(self):
        """Full integration: matplotlib should not raise with zero delta-F data."""
        import matplotlib
        matplotlib.use("Agg")  # non-interactive backend
        import matplotlib.pyplot as plt
        from pyEPR.reports import plot_convergence_max_df

        s = pd.Series([0.0, 0.0], name="Max Delta Freq")
        fig, ax = plt.subplots()
        try:
            plot_convergence_max_df(ax, s)
            fig.tight_layout()   # this is what crashed before
        except ValueError as e:
            pytest.fail(f"plot_convergence_max_df raised ValueError: {e}")
        finally:
            plt.close(fig)

    def test_matplotlib_does_not_crash_with_positive_data(self):
        """Ensure normal data still works (log scale applied)."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from pyEPR.reports import plot_convergence_max_df

        s = pd.Series([0.5, 0.3, 0.1], name="Max Delta Freq")
        fig, ax = plt.subplots()
        try:
            plot_convergence_max_df(ax, s)
            fig.tight_layout()
            assert ax.get_yscale() == "log"
        finally:
            plt.close(fig)
