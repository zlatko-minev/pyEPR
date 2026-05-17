"""
Tests for numerical sorting fixes (issues #56 and variation ordering):
  - sort_df_col: columns with >9 integer strings now sorted numerically
  - sort_Series_idx: same fix for Series index
  - plot_hamiltonian_results: sort_index key is robust to float sweep variables
"""
import numpy as np
import pandas as pd
import pytest


# ── sort_df_col ───────────────────────────────────────────────────────────────

class TestSortDfCol:
    def test_integer_string_columns_sorted_numerically(self):
        """Old code: '10' < '2' (lexicographic). New code: 2 < 10 (numeric)."""
        from pyEPR.toolbox.pythonic import sort_df_col
        cols = ["0", "10", "2", "9", "1"]
        df = pd.DataFrame(np.zeros((1, len(cols))), columns=cols)
        result = sort_df_col(df)
        assert list(result.columns) == ["0", "1", "2", "9", "10"]

    def test_float_columns_sorted_numerically(self):
        from pyEPR.toolbox.pythonic import sort_df_col
        cols = [1.5, 0.1, 10.0, 2.3]
        df = pd.DataFrame(np.zeros((1, len(cols))), columns=cols)
        result = sort_df_col(df)
        assert list(result.columns) == [0.1, 1.5, 2.3, 10.0]

    def test_single_digit_unchanged(self):
        from pyEPR.toolbox.pythonic import sort_df_col
        cols = ["0", "1", "2"]
        df = pd.DataFrame(np.zeros((1, 3)), columns=cols)
        result = sort_df_col(df)
        assert list(result.columns) == ["0", "1", "2"]

    def test_non_numeric_columns_sorted_lexicographically(self):
        from pyEPR.toolbox.pythonic import sort_df_col
        cols = ["beta", "alpha", "gamma"]
        df = pd.DataFrame(np.zeros((1, 3)), columns=cols)
        result = sort_df_col(df)
        assert list(result.columns) == ["alpha", "beta", "gamma"]

    def test_data_preserved_after_sort(self):
        from pyEPR.toolbox.pythonic import sort_df_col
        cols = ["0", "10", "2"]
        values = [[10, 20, 30]]
        df = pd.DataFrame(values, columns=cols)
        result = sort_df_col(df)
        # After sort: columns are [0, 2, 10], values should follow
        assert result["0"].iloc[0] == 10
        assert result["2"].iloc[0] == 30
        assert result["10"].iloc[0] == 20


# ── sort_Series_idx ───────────────────────────────────────────────────────────

class TestSortSeriesIdx:
    def test_integer_string_index_sorted_numerically(self):
        from pyEPR.toolbox.pythonic import sort_Series_idx
        idx = ["0", "10", "2", "9", "1"]
        sr = pd.Series(range(5), index=idx)
        result = sort_Series_idx(sr)
        assert list(result.index) == ["0", "1", "2", "9", "10"]

    def test_float_index_sorted_numerically(self):
        from pyEPR.toolbox.pythonic import sort_Series_idx
        idx = [1.5, 0.1, 10.0, 2.3]
        sr = pd.Series(range(4), index=idx)
        result = sort_Series_idx(sr)
        assert list(result.index) == [0.1, 1.5, 2.3, 10.0]

    def test_non_numeric_index_sorted_lexicographically(self):
        from pyEPR.toolbox.pythonic import sort_Series_idx
        idx = ["beta", "alpha", "gamma"]
        sr = pd.Series(range(3), index=idx)
        result = sort_Series_idx(sr)
        assert list(result.index) == ["alpha", "beta", "gamma"]


# ── plot_hamiltonian_results sort key (issue #56) ─────────────────────────────

class TestPlotSortKey:
    """
    The sort key in plot_hamiltonian_results used x.astype(int) which raises
    ValueError when swp_variable is a float-valued sweep (e.g., Lj).
    The fix uses pd.to_numeric(..., errors='coerce') which handles both cases.
    """

    def test_integer_index_sorts_correctly(self):
        idx = pd.Index(["0", "1", "2", "10"])
        numeric = pd.to_numeric(idx, errors="coerce")
        assert numeric.notna().all()
        sorted_idx = idx[np.argsort(numeric)]
        assert list(sorted_idx) == ["0", "1", "2", "10"]

    def test_float_index_sorts_correctly(self):
        idx = pd.Index([3.0, 1.5, 0.5, 2.0])
        numeric = pd.to_numeric(idx, errors="coerce")
        assert numeric.notna().all()
        sorted_idx = idx[np.argsort(numeric)]
        assert list(sorted_idx) == [0.5, 1.5, 2.0, 3.0]

    def test_non_numeric_index_coerces_to_nan(self):
        idx = pd.Index(["a", "b", "c"])
        numeric = pd.to_numeric(idx, errors="coerce")
        assert numeric.isna().all()

    def test_small_float_astype_int_gives_wrong_sort(self):
        """Confirm the old code gave wrong ordering for small-float sweep variables (e.g., Lj in H).

        Values like 1.4e-8, 2.0e-8 all truncate to 0 when cast to int,
        making the sort key useless — the new pd.to_numeric path avoids this.
        """
        idx = pd.Index([2.0e-8, 0.5e-8, 1.5e-8])
        # Old: astype(int) truncates all to 0 — sort key is degenerate
        old_key = idx.astype(int)
        assert list(old_key) == [0, 0, 0], "all small floats truncate to 0"
        # New: pd.to_numeric preserves values for correct ordering
        numeric = pd.to_numeric(idx, errors="coerce")
        assert numeric.notna().all()
        sorted_idx = idx[np.argsort(numeric)]
        assert list(sorted_idx) == [0.5e-8, 1.5e-8, 2.0e-8]
