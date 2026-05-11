# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install (editable)
pip install -e ".[test]"

# Run all non-HFSS tests
pytest

# Run a single test file
pytest tests/test_solution_types.py -v

# Run tests that require live Ansys HFSS (skipped in CI by default)
pytest -m hfss

# Lint (CI runs errors-only; locally use full output)
pylint pyEPR/                    # full report
pylint --errors-only pyEPR/      # CI mode

# Build distribution
pip install build
python -m build
```

Key pytest config (`pyproject.toml`): default `addopts = "-m 'not hfss'"`, so HFSS tests never run in CI unless explicitly requested.

## Architecture

### Package layout

```
pyEPR/
  __init__.py              # public API, version string, import checks
  ansys.py                 # COM wrappers for Ansys HFSS / Q3D
  solution_types.py        # canonical names, alias frozensets, normalize()
  core.py                  # ProjectInfo + re-exports
  project_info.py          # ProjectInfo dataclass
  core_distributed_analysis.py   # DistributedAnalysis (field/EPR extraction)
  core_quantum_analysis.py       # QuantumAnalysis (Hamiltonian diagonalization)
  calcs/                   # EPR math utilities
  toolbox/                 # logging, plotting, pandas helpers
```

### COM wrapper hierarchy (`ansys.py`)

```
HfssApp          – wraps oDesktop application object
  HfssDesktop    – wraps desktop; open/close projects
    HfssProject  – wraps project; holds _ansys_version cache
      HfssDesign – wraps design; normalises solution_type at __init__
        HfssSetup (abstract)
          HfssDMSetup   – driven-modal setup
          HfssDTSetup   – driven-terminal setup
          HfssEMSetup   – eigenmode setup
          AnsysQ3DSetup – Q3D setup
```

All COM objects are accessed via a thin `COMWrapper` base that delegates attribute access to `_wrapped_COM_object`.

`HfssProject._ansys_version` is cached at init from `self.parent.version` (format `"YYYY.N"`, e.g. `"2024.1"`). String comparison (`>=`) works for gating version-specific behaviour.

### EPR analysis pipeline

1. `ProjectInfo` (alias `Project_Info`) — user-facing config: HFSS project path, junction list, sweep settings.
2. `DistributedAnalysis` (`pyEPR_HFSSAnalysis`) — connects to HFSS, runs field extraction, computes EPR participation ratios per mode/junction.
3. `QuantumAnalysis` (`pyEPR_Analysis`) — loads EPR results, performs numerical diagonalization (via qutip), returns Hamiltonian parameters (χ, α, g).

### solution_types module

`pyEPR/solution_types.py` is the single source of truth for HFSS solution-type string handling. Import from here rather than duplicating alias sets:

```python
from pyEPR.solution_types import normalize, DRIVEN_MODAL_NAMES, is_drivenmodal
```

## HFSS version compatibility

Two layers of fixes are required for AEDT compatibility:

### READ side (normalisation at init)

`HfssDesign.__init__` calls `normalize(design.GetSolutionType())` on the raw COM string. The `normalize()` function uses substring matching (same approach as PyAEDT):

| Raw string (AEDT 2021.2+)       | Canonical output   |
|----------------------------------|--------------------|
| `"HFSS Modal Network"`           | `"DrivenModal"`    |
| `"HFSS Hybrid Modal Network"`    | `"DrivenModal"`    |
| `"HFSS Terminal Network"`        | `"DrivenTerminal"` |
| `"HFSS Hybrid Terminal Network"` | `"DrivenTerminal"` |

All downstream pyEPR code always receives the canonical pre-2021.2 string.

### CREATE side (AEDT 2024.1+ hybrid default)

Starting with AEDT 2024.1, `InsertDesign("HFSS", name, "DrivenModal", "")` silently creates an **HFSS Hybrid Modal Network** instead of a plain DrivenModal. Fix: call `SetSolutionType` immediately after creation.

```python
# Use these helpers instead of new_design() directly:
project.new_dm_design(name)   # DrivenModal  – calls SetSolutionType on AEDT >= 2024.1
project.new_dt_design(name)   # DrivenTerminal – same fix
project.new_em_design(name)   # Eigenmode – no hybrid variant exists, no workaround needed
```

When checking Ansys scripting API for new issues, look at:
- `InsertDesign` — signature changes between AEDT versions
- `SetSolutionType` — accepted string values differ by AEDT version
- `GetSolutionType` — return value changes (new aliases are added, not removed)

The Ansys AEDT scripting guide (IronPython/CPython) is the authoritative source; PyAEDT source code is a useful secondary reference for how they handle the same changes.

## Release workflow

1. Bump `__version__` in `pyEPR/__init__.py` — `pyproject.toml` reads version dynamically from there.
2. Commit and push to master (via PR, not direct push).
3. Create a GitHub Release (tag e.g. `v0.9.4`). The `publish-to-pypi.yml` workflow triggers on `release: created` and uses OIDC Trusted Publishing (no API token needed; configured on PyPI under Trusted Publishers).
4. Verify the published version on PyPI matches the intended `__version__`. If `__version__` was not bumped before tagging, the wheel will carry the old version.

CI workflow (`ci.yaml`) runs on every push: pylint (errors-only), pytest (no hfss marker), and docs build.

## Backwards compatibility

- All public API must remain stable across minor versions. New helpers like `new_dm_design` / `new_dt_design` are additive; they do not change existing `new_design()` behaviour.
- `Project_Info`, `pyEPR_HFSSAnalysis`, `pyEPR_Analysis` are deprecated aliases kept in `__init__.py` for backwards compatibility — do not remove them.
- The `solution_types` module is new (added in 0.9.x) — downstream packages (qiskit-metal) can import from it without importing the full COM stack.
- When adding new HFSS functionality, always check that the COM calls gracefully handle older AEDT versions where the API may not exist.

## Testing notes

- Mark any test requiring a live Ansys HFSS session with `@pytest.mark.hfss` and they will be skipped in CI automatically.
- `tests/correct_results.pkl` and `tests/data*.npz` are reference fixtures for numerical regression tests.
- qutip 5.x changed `ket.dag() * ket` to return a complex scalar instead of a 1×1 Qobj — guard with `hasattr(inner, "norm")` when writing quantum analysis tests.
