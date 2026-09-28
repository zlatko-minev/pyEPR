"""pyEPR must not use names NumPy 2.0 removed from its main namespace.

NumPy 2 raises ``AttributeError`` on first use of such a name, so a
removed alias in a rarely run helper fails only when a user reaches it:
``print_matrix`` (``np.mat``) runs at the end of every
``analyze_variation(print_result=True)``, the default.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pyEPR.toolbox.pythonic import df_find_index, df_interpolate_value, print_matrix

# Removed in NumPy 2.0 (NEP 52); each has a replacement that also works in 1.x.
REMOVED = (
    "NaN Inf Infinity infty NINF PINF NZERO PZERO float_ complex_ cfloat "
    "longfloat clongfloat longcomplex singlecomplex string_ unicode_ mat "
    "product cumproduct sometrue alltrue round_ asfarray issubclass_ "
    "find_common_type msort"
).split()


def test_no_removed_numpy_names_in_source():
    pattern = re.compile(r"\bnp\.(" + "|".join(REMOVED) + r")\b(?!\w)")
    root = Path(__file__).resolve().parent.parent / "pyEPR"
    hits = [
        f"{path.relative_to(root)}:{n}: {line.strip()}"
        for path in sorted(root.rglob("*.py"))
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if pattern.search(line) and not line.lstrip().startswith("#")
    ]
    assert not hits, "\n".join(hits)


def test_print_matrix(capsys):
    print_matrix(np.array([[1.0, 2.0], [3.0, 4.0]]), frmt="{:5.1f}")
    print_matrix([1.0, 2.0], frmt="{:5.1f}")
    assert capsys.readouterr().out == "   1.0  2.0\n   3.0  4.0\n   1.0  2.0\n"


def test_df_interpolate_value():
    s = pd.Series([1.0, 3.0], index=[10.0, 20.0])
    value, _ = df_interpolate_value(s, 15.0)
    assert value == 2.0


def test_df_find_index():
    # Frequencies (values) 5 and 6 at Lj (index) 10 and 20: the Lj for 5.5 is
    # interpolated. The range check used to compare 5.5 with the Lj range, so
    # an in-range target went to the extrapolation branch instead.
    s = pd.Series([5.0, 6.0], index=[10.0, 20.0])
    value, _ = df_find_index(s, 5.5)
    assert value == 15.0


def test_df_find_index_extrapolates_outside_the_values():
    s = pd.Series([5.0, 6.0], index=[10.0, 20.0])
    value, _ = df_find_index(s, 7.0, degree=1)
    assert value == pytest.approx(30.0)
