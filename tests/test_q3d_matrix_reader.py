"""Reading a Q3D matrix export (``AnsysQ3DSetup.load_q3d_matrix``).

The sample is the one in ``_readin_Q3D_matrix``'s docstring. The reader
used ``pd.read_csv(delim_whitespace=True)``, deprecated in pandas 2.2 and
removed in 3.0; FutureWarning is an error here so a deprecated argument
fails on pandas 2.x already.

No COM or Ansys is involved: ``_readin_Q3D_matrix`` parses the exported
text file, so this test needs no ``hfss`` mark.
"""
import warnings

import numpy as np
import pytest

from pyEPR.ansys import AnsysQ3DSetup

NAMES = [
    "ground_plane",
    "Q1_bus_Q0_connector_pad",
    "Q1_bus_Q2_connector_pad",
    "Q1_pad_bot",
    "Q1_pad_top1",
    "Q1_readout_connector_pad",
]
CAP = [
    [2.8829e-13, -3.254e-14, -3.1978e-14, -4.0063e-14, -4.3842e-14, -3.0053e-14],
    [-3.254e-14, 4.7257e-14, -2.2765e-16, -1.269e-14, -1.3351e-15, -1.451e-16],
    [-3.1978e-14, -2.2765e-16, 4.5327e-14, -1.218e-15, -1.1552e-14, -5.0414e-17],
    [-4.0063e-14, -1.269e-14, -1.218e-15, 9.5831e-14, -3.2415e-14, -8.3665e-15],
    [-4.3842e-14, -1.3351e-15, -1.1552e-14, -3.2415e-14, 9.132e-14, -1.0199e-15],
    [-3.0053e-14, -1.451e-16, -5.0414e-17, -8.3665e-15, -1.0199e-15, 3.9884e-14],
]


def _table(rows):
    header = "\t" + "\t".join(NAMES)
    body = ["\t".join([name] + [f"{v:.5G}" for v in row]) for name, row in zip(NAMES, rows)]
    return "\n".join([header] + body)


def _export(tmp_path):
    text = "\n".join(
        [
            "DesignVariation:$BBoxL='650um' $boxH='750um' Lj_1='13nH'",
            "Setup1:LastAdaptive",
            "Problem Type:C",
            "C Units:farad, G Units:mSie",
            "Reduce Matrix:Original",
            "Frequency: 5.5E+09 Hz",
            "",
            "Capacitance Matrix",
            _table(CAP),
            "",
            "Conductance Matrix",
            _table([[0] * 6] * 6),
            "",
        ]
    )
    path = tmp_path / "q3d_export.txt"
    path.write_text(text)
    return path


def test_readin_q3d_matrix(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        df_cmat, units, variation, df_cond, units_cond = AnsysQ3DSetup._readin_Q3D_matrix(
            _export(tmp_path)
        )
    assert list(df_cmat.index) == NAMES
    assert list(df_cmat.columns) == NAMES
    np.testing.assert_allclose(df_cmat.values, CAP, rtol=1e-12)
    assert units == "farad"
    assert variation == "$BBoxL='650um' $boxH='750um' Lj_1='13nH'"
    assert df_cond.shape == (6, 6)
    assert units_cond == "mSie"


def test_load_q3d_matrix_converts_units(tmp_path):
    df_cmat, units, _, _ = AnsysQ3DSetup.load_q3d_matrix(_export(tmp_path), user_units="fF")
    assert units == "fF"
    assert df_cmat.loc["ground_plane", "ground_plane"] == pytest.approx(288.29)
