# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Agent context — read these before starting work

| File | When to read |
|------|-------------|
| `.claude/context/lessons-learned.md` | Before touching docs, CI, Sphinx config, or Ansys version code. Contains every hard-won fix from real build failures and regressions. |
| `.claude/context/ecosystem.md` | Before changing public API, release timing, imports, or anything that could affect downstream users. Explains who uses pyEPR, the quantum-metal relationship, and the no-HFSS adoption path. |

## Slash commands

| Command | What it does |
|---------|-------------|
| `/health-check` | Full 10-section maintenance audit: docs, tests, deps, API stability, README/PyPI, release, CI, tutorials, ecosystem, code hygiene. |
| `/release` | Step-by-step release workflow: verify readiness → version bump → PR → merge → tag → PyPI → post-release check. |
| `/docstring-audit` | Module-by-module NumPy docstring audit with priority ordering and RST pitfall reminders. |

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

# Build docs locally (must copy notebooks first — see Docs section)
rm -rf docs/source/_tutorial_notebooks
cp -r _tutorial_notebooks docs/source/_tutorial_notebooks
cd docs && make html
# open docs/build/html/index.html

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

## Testability boundary

This is the most important constraint for any agent working in this repo.

| Area | Testable in CI | Notes |
|------|---------------|-------|
| `calcs/` | Yes | Pure math, numpy/scipy/qutip only |
| `solution_types.py` | Yes | No external deps |
| `QuantumAnalysis` (post-HDF5) | Yes | Uses fixtures in `tests/` |
| `toolbox/` | Yes | Logging, plotting, pandas helpers |
| `project_info.py` | Yes | Config object, no COM |
| `ansys_pyaedt.py` (physics + unit parsing + lock targeting) | Yes | `compute_p_mj`, `_parse_henries`, `_owning_session_from_lock` are pure Python — tested in `tests/test_ansys_pyaedt.py` without PyAEDT installed |
| `ansys_pyaedt.py` (live gRPC extraction) | **No** | Requires live AEDT session; mark with `@pytest.mark.hfss` |
| `ansys.py` | **No** | Requires live Ansys HFSS COM session |
| `DistributedAnalysis` (field extraction) | **No** | Requires live HFSS solve |
| `HfssDesign`, `HfssSetup`, etc. | **No** | COM-only |

**Never write a test that calls into `ansys.py` without `@pytest.mark.hfss`.** Tests without that mark run in CI on every PR and will fail without a physical Ansys licence.

For the `ansys.py` code: you can read, reason about, and document it. You can flag logic errors or version-compatibility issues. But do not modify COM-facing behaviour without a human who can validate against a real AEDT session.

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

## Documentation standards

pyEPR uses the **PyData Sphinx Theme** with `sphinx-design` and `myst-nb`. These were chosen to match the broader scientific Python ecosystem (NumPy, SciPy, pandas, QuTiP) and to support rich landing pages and notebook galleries without requiring Ansys.

### Hard rules

- **Zero warnings on `make html`** — this is a CI gate. Every new docstring, RST file, or notebook added must not introduce warnings. Run `cd docs && make html 2>&1 | grep -c WARNING` before committing docs changes; it must return `0`.
- **No print() in analysis code** — use `logger = logging.getLogger(__name__)`. The toolbox provides logging helpers.
- **NumPy-style docstrings** — all public methods in `core_*.py`, `project_info.py`, and `calcs/` must use the NumPy docstring convention (Parameters / Returns / Raises sections).

### Notebook copy step (critical)

Sphinx does not follow symlinks outside its source root. `docs/source/_tutorial_notebooks` is a git symlink to `../../_tutorial_notebooks/`. **Always copy before building:**

```bash
rm -rf docs/source/_tutorial_notebooks
cp -r _tutorial_notebooks docs/source/_tutorial_notebooks
```

This step is required in:
- Local `make html` (do it manually)
- `.readthedocs.yml` `pre_build` job
- `.github/workflows/ci.yaml` `test_docs` job

If you add a new CI job that builds docs, add the copy step.

### RST / docstring pitfalls

These have each caused real build failures. Know them before editing docs or docstrings:

- **Pipe characters in docstrings** — `|x⟩` bra-ket notation triggers RST substitution reference parsing. Use `:math:`|x\\rangle`` instead.
- **`*args` / `**kwargs` in docstrings** — bare `*` and `**` inside paragraphs trigger RST emphasis markup. Wrap in double backticks or a code block.
- **Separator lines in docstrings** — a line of `----` at the start of a docstring section triggers "Transition must be child of section" warning. Use a proper RST section header instead.
- **Indented first lines** — a docstring whose first line is indented relative to the opening `"""` creates a block-quote RST node, which Sphinx renders oddly.
- **Duplicate object descriptions** — if a class appears in both `key_classes_reference.rst` (via `autoclass`) and `api/` (via `automodule`), Sphinx emits a duplicate-object warning. Use `:no-index:` on one, or replace `autoclass` with a plain cross-reference link (preferred).
- **Unknown directives** — `.. tabs::` (sphinx-tabs) and `.. tab::` are not installed; the project uses `.. tab-set::` / `.. tab-item::` (sphinx-design). Similarly, `nbsphinx` is not installed; notebooks are handled by `myst-nb`.
- **`.. code-block::` syntax** — must be `.. code-block:: python` with a space before the language. `..codeblock python` is silently ignored.

## Release workflow

1. Bump `__version__` in `pyEPR/__init__.py` — `pyproject.toml` reads version dynamically from there.
2. Commit and push to master via PR (never direct-push to master).
3. After the PR merges, create a GitHub Release with tag `vX.Y.Z` (e.g. `v0.9.6`). The `publish-to-pypi.yml` workflow triggers on `release: created` and uses OIDC Trusted Publishing — no API token needed.
4. Verify the published version on PyPI matches `__version__`. If the version was not bumped before tagging, the wheel will carry the wrong version.

CI (`ci.yaml`) runs on every push: pylint (errors-only), pytest (non-HFSS), and docs build.

**The tag must be created on the commit that has the bumped version number.** Create the GitHub Release pointing at master only after the version-bump PR has merged.

## Backwards compatibility

- All public API must remain stable across minor versions. New helpers are additive; they do not change existing behaviour.
- `Project_Info`, `pyEPR_HFSSAnalysis`, `pyEPR_Analysis` are deprecated aliases kept in `__init__.py` — do not remove them. Downstream packages (including qiskit-metal / quantum-metal) import these.
- The `solution_types` module was designed so downstream packages can import it without pulling in the COM stack (`import win32com` is Windows-only). Never add a top-level import of COM-related modules to `solution_types.py` or `calcs/`.
- When adding HFSS functionality, always check that COM calls gracefully handle older AEDT versions where the API may not exist. Guard with `self._ansys_version >= "YYYY.N"` checks.

## Ecosystem context

pyEPR is a dependency of **quantum-metal** (formerly qiskit-metal), IBM's open-source quantum chip design framework. That package declares `pyEPR-quantum >= 0.9.5` in its `pyproject.toml`. Changes to pyEPR's public API, import structure, or `solution_types` module may silently break quantum-metal until their next release — be conservative.

The primary adoption path for users without Ansys is:
1. Tutorial 6 (`_tutorial_notebooks/`) — numerical EPR without HFSS
2. `QuantumAnalysis` loaded from a pre-computed HDF5 file
3. `calcs/` subpackage for direct formula access

These paths must remain importable on Linux and macOS without `win32com` or any Windows-only dependency.

## Testing notes

- Mark any test requiring a live Ansys HFSS session with `@pytest.mark.hfss` — CI skips these automatically.
- `tests/correct_results.pkl` and `tests/data*.npz` are reference fixtures for numerical regression tests.
- qutip 5.x changed `ket.dag() * ket` to return a complex scalar instead of a 1×1 Qobj — guard with `hasattr(inner, "norm")` when writing quantum analysis tests.
- When adding a new feature, add a test that would have caught the most obvious misuse. The test suite is the primary protection against regressions from upstream dependency changes (qutip, numpy, pandas, scipy all release breaking changes).
