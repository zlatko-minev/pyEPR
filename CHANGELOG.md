# Changelog

All notable changes to pyEPR are documented here.
Versions follow [Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### New features

- **PyAEDT (gRPC) HFSS backend** (`pyEPR.ansys_pyaedt.PyaedtDistributedAnalysis`).
  Runs the same Energy-Participation-Ratio field extraction as
  `DistributedAnalysis`, but through Ansys's official PyAEDT library
  (`pyaedt`) entirely over gRPC — no COM / `pywin32`. PyAEDT can attach
  to an already-running AEDT session that owns the project (via its `.aedt.lock`),
  avoiding stale-session and project-locked errors common with COM. The extracted
  participations feed pyEPR's own physics (`epr_to_zpf`,
  `epr_numerical_diagonalization`, `QuantumAnalysis`) unchanged; validated
  digit-for-digit against the COM path (`p_mj = 0.9755` on a demo transmon).
  PyAEDT is imported lazily behind the new `[pyaedt]` optional-dependency extra,
  so `import pyEPR` never requires it. See `_tutorial_notebooks/` and the new
  docs page.

## [0.9.5] — 2026-05-17

### New features

- **Exact cosine diagonalization** (`use_full_cos=True` in `epr_numerical_diagonalization` and `analyze_variation`).
  For strongly anharmonic circuits such as fluxonium (φ_zpf ≳ 1) the truncated Taylor series diverges;
  the new path evaluates cos(φ) − 1 + φ²/2 exactly via matrix exponential.
  See `MatrixOps.cos_full_correction` and arXiv:2411.15039.

- **Generic junction potential** (`make_nonlinear_potential`).
  Converts any scalar Python function V(φ) — including flux-biased junctions,
  asymmetric SQUIDs, or exotic elements — into the EPR-compatible operator-valued
  correction, automatically expanding around the bias minimum and normalising by |V′′(φ₀)|.
  Also adds `MatrixOps.apply_scalar_function` for evaluating arbitrary scalar functions
  on Hermitian operators via eigendecomposition.

- **Public Q3D convergence plot API.**
  `pyEPR.reports.plot_q3d_convergence_main` and `plot_q3d_convergence_chi_f`
  are now public, documented functions. The old underscore names remain as
  backward-compatible aliases.

- **Tutorial 5** — *Generic junction potential and fluxonium EPR* (`_tutorial_notebooks/`):
  transmon vs. exact-cosine convergence, fluxonium regime, flux-biased and asymmetric-SQUID
  examples, potential visualisation. Fully self-contained (no HFSS required).

### Improvements

- **Numerical sort for >9 variations** (`sort_df_col`, `sort_Series_idx`):
  fixed long-standing bug where variation "10" sorted before "2" (lexicographic order).
  Now uses `pd.to_numeric` for correct numeric ordering.

- **Float sweep variable sort** (`plot_hamiltonian_results`):
  the old `x.astype(int)` sort key truncated all small floats (e.g., Lj = 1.4e-8) to 0,
  making the sort degenerate. Fixed with `pd.to_numeric(..., errors='coerce')`.

- **DeprecationWarnings on legacy aliases.**
  `pyEPR.Project_Info`, `pyEPR.pyEPR_HFSSAnalysis`, `pyEPR.pyEPR_Analysis` now emit
  `DeprecationWarning` on first access (PEP 562 module-level `__getattr__`).
  The aliases themselves continue to work for backward compatibility.

- **Logger replaces print.**
  28 `print()` calls in `DistributedAnalysis` and 12 in `QuantumAnalysis` replaced
  with structured `logger.info/debug/warning/error` calls. Downstream code can now
  control pyEPR verbosity via the standard Python logging hierarchy.

- **NumPy-style docstrings** added to `QuantumAnalysis`, `DistributedAnalysis`,
  `ProjectInfo`, `get_frequencies`, `get_chis`, `get_quality_factors`,
  `get_participations`, `analyze_variation`, `analyze_all_variations`,
  `plot_hamiltonian_results`, `epr_numerical_diagonalization`, `black_box_hamiltonian`,
  `make_nonlinear_potential`, and all new `MatrixOps` methods.

### Bug fixes

- Rebase conflict in 0.9.5 docstring pass restored missing `QuantumAnalysis`-level
  docstrings that were dropped when cherry-picking onto master.

### Upgrading from 0.9.4

No breaking changes. All existing call signatures are unchanged.

- `epr_numerical_diagonalization` and `analyze_variation` accept a new
  `use_full_cos=False` keyword; old callers are unaffected.
- `make_nonlinear_potential` and `cos_full_correction` are new additions to
  `pyEPR.calcs.back_box_numeric`.
- `Project_Info`, `pyEPR_HFSSAnalysis`, `pyEPR_Analysis` now emit
  `DeprecationWarning`. They still work — update to `ProjectInfo`,
  `DistributedAnalysis`, `QuantumAnalysis` at your convenience.
- `_plot_q3d_convergence_main` and `_plot_q3d_convergence_chi_f` remain as
  aliases for the new public names; no caller changes required.

---

## [0.9.4] — 2024

### New features

- `HfssDesign` context manager (`with design:`) for automatic resource cleanup.
- `ProjectInfo.junctions` keyword argument for cleaner junction specification.
- `HfssDesign.get_variable_value()` helper.
- `new_dt_design()` (DrivenTerminal) added alongside existing `new_dm_design()`.

### Bug fixes

- Fix `HfssDesign.new_dm_design()` version check and docstrings (#128, #162, #169, #137).
- Fix assigning Y coordinate to `YAxisZvec` in `create_relative_coordinate_system_both`.
- Fix missing return in `get_excitations`.
- Fix pandas `FutureWarning` and per-curve line-width loop in reports (#140, #141).
- Fix `InsertDesign` / `SetSolutionType` for AEDT 2024.1 hybrid-modal default (#182).

---

## [0.9.3] — 2023

- Bump for PyPI compatibility fixes (`pyproject.toml` license field).
- CI: opt into Node.js 24 for GitHub Actions.
- Remove stale `.gitmodules` submodule entry.

---

## [0.9.2] and earlier

See the [GitHub commit history](https://github.com/zlatko-minev/pyEPR/commits/master)
for changes prior to 0.9.3.
