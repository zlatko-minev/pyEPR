# /health-check

Perform a systematic maintenance and health audit of the pyEPR repository. Work through each section below in order. For each item: check the current state, report what you find, fix anything clearly broken, and flag anything that needs human judgment (especially anything touching Ansys/HFSS COM code).

At the end, produce a prioritised punch list: what was fixed, what needs attention, and what was skipped and why.

---

## 1. Docs build — zero warnings required

```bash
rm -rf docs/source/_tutorial_notebooks
cp -r _tutorial_notebooks docs/source/_tutorial_notebooks
cd docs && make html 2>&1 | tee /tmp/sphinx_build.log
grep -c "WARNING\|ERROR" /tmp/sphinx_build.log
```

- If warning count > 0: read each warning and fix it before proceeding. Common causes are listed in CLAUDE.md under "RST / docstring pitfalls".
- Check that all six tutorial notebooks appear in `tutorials.rst` and render without toctree warnings.
- Check that the landing page (`index.rst`) hero image and xmon GIF render — both require absolute `_static/` paths, not symlinks.
- Verify `docs/source/conf.py` extensions list matches what's installed in `[docs]` extras in `pyproject.toml`. Any `Unknown directive` or `No module named` error means a mismatch.

## 2. Test suite — baseline must hold

```bash
pytest -q 2>&1 | tail -20
```

- Zero failures, zero errors on the non-HFSS suite. If anything is failing, it is a blocker — fix it before anything else.
- Check that `@pytest.mark.hfss` tests are properly skipped (not erroring). Run `pytest -q --co -m hfss` to list them; they should show as collected but never run in CI.
- Review `tests/` for tests that call into `ansys.py` without the `hfss` mark — this would cause silent failures for users without Ansys. Add the mark if missing.
- Check `tests/correct_results.pkl` and `tests/data*.npz` — if any are missing from the repo, the numerical regression tests will silently pass vacuously.

## 3. Dependency health

```bash
pip install -e ".[test]"
python -W error::DeprecationWarning -m pytest tests/ -q 2>&1 | grep -i deprecat | head -20
```

- Fix any DeprecationWarnings that come from pyEPR's own code (not from third-party libraries).
- Check `pyproject.toml` for version pins. The following are intentional and should not be changed without explicit review:
  - Any upper-bound pin on a scientific library (numpy, scipy, pandas) — these exist because upstream broke something.
- Verify that the package installs cleanly on Python 3.10 and 3.12:
  ```bash
  python --version  # check current
  python -c "import pyEPR; print(pyEPR.__version__)"
  ```
- Check for any import of `win32com`, `pythoncom`, or `pywintypes` at the top level of any module outside `ansys.py`. These are Windows-only and will break Linux/macOS installs.

## 4. API stability

```bash
python -c "
import pyEPR as epr
# All of these must import without error or DeprecationWarning
_ = epr.ProjectInfo
_ = epr.DistributedAnalysis
_ = epr.QuantumAnalysis
_ = epr.Project_Info          # deprecated alias
_ = epr.pyEPR_HFSSAnalysis    # deprecated alias
_ = epr.pyEPR_Analysis        # deprecated alias
from pyEPR.solution_types import normalize, DRIVEN_MODAL_NAMES, is_drivenmodal
print('All public API imports OK')
"
```

- If any alias raises `ImportError` or `AttributeError`, fix it — downstream packages depend on them.
- Check that `solution_types` is importable independently without triggering COM imports:
  ```bash
  python -c "import sys; import pyEPR.solution_types; print([m for m in sys.modules if 'win32' in m or 'ansys' in m.lower()])"
  # should print []
  ```

## 5. README and PyPI rendering

- Open `README.md` and check every image URL. Images must use absolute `https://raw.githubusercontent.com/zlatko-minev/pyEPR/master/...` URLs — relative paths do not render on PyPI.
- Check that the xmon GIF and any other media files referenced in README actually exist in the `imgs/` directory.
- Verify the PyPI badge URLs (shields.io) are well-formed and point to the correct package name (`pyEPR-quantum`, not `pyepr` or `pyEPR`).
- Check that `docs/source/about.rst` has no dead image URLs (frapsoft.com and rawgit.com are known-dead badge hosts — replace with sphinx-design badges if found).

## 6. Version and release consistency

```bash
python -c "import pyEPR; print(pyEPR.__version__)"
grep 'version' pyproject.toml
git tag --list | sort -V | tail -5
```

- `__version__` in `pyEPR/__init__.py` must match the most recent git tag.
- `pyproject.toml` should read version dynamically: `dynamic = ["version"]` with `[tool.setuptools.dynamic] version = {attr = "pyEPR.__version__"}`. If it's hardcoded, flag this — version drift between the two is a real release bug.
- Check that `CHANGELOG.md` (if it exists) has an entry for the current version.
- Check the GitHub Releases page (via `gh release list` or MCP tool) to confirm the latest release tag points to a commit that has the matching `__version__`. A tag on the wrong commit publishes the wrong version to PyPI.

## 7. CI workflow correctness

Read `.github/workflows/ci.yaml` and check:

- The `test_docs` job includes the notebook copy step before `make html`:
  ```yaml
  - name: Copy tutorial notebooks into docs source
    run: |
      rm -rf docs/source/_tutorial_notebooks
      cp -r _tutorial_notebooks docs/source/_tutorial_notebooks
  ```
  If this step is missing, the docs CI job will report false-pass with toctree warnings.
- The `test` job matrix covers at least Python 3.10 and 3.12 on Ubuntu and macOS.
- The `pylint` job uses `--errors-only` (not full report) to keep CI signal-to-noise ratio high.
- The `publish-to-pypi.yml` workflow triggers on `release: created` (not `push`). If it triggers on push, it will attempt to publish on every commit.

## 8. Tutorials and onboarding experience

```bash
ls _tutorial_notebooks/
```

- All six tutorial notebooks must be present: Tutorials 1–6.
- Tutorials 3, 5, and 6 are the no-HFSS path — these are the onboarding notebooks for users without Ansys. Check that they run correctly:
  ```bash
  jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=120 \
    "_tutorial_notebooks/Tutorial 6. EPR without HFSS — purely numerical workflow.ipynb"
  ```
- Check that Tutorial 6 has a working Binder badge in `docs/source/tutorials.rst`. The URL should point to the correct branch (`master`) and notebook path.
- Check that `docs/source/without_hfss.rst` exists and accurately describes the no-HFSS workflow. This is the primary landing page for users without Ansys licences and must be kept current.

## 9. Ecosystem compatibility

```bash
pip show pyEPR-quantum | grep Version
```

- Check the current version against what quantum-metal (formerly qiskit-metal) declares as its minimum: `pyEPR-quantum >= 0.9.5`. If pyEPR's current version is behind this, users of quantum-metal on older pyEPR will get import errors.
- Run a quick compatibility smoke test:
  ```python
  from pyEPR.solution_types import normalize, DRIVEN_MODAL_NAMES, is_drivenmodal
  from pyEPR import ProjectInfo, DistributedAnalysis, QuantumAnalysis
  # These are the imports quantum-metal uses
  ```
- If any of the above fail, that is a breaking change for downstream users — treat as P0.

## 10. Code hygiene

```bash
# Check for remaining print() calls in analysis code (not ansys.py, not toolbox)
grep -rn "^\s*print(" pyEPR/ --include="*.py" \
  | grep -v ansys.py | grep -v toolbox/ | grep -v __pycache__ | grep -v ".pyc"
```

- Any `print()` in `core_distributed_analysis.py`, `core_quantum_analysis.py`, `project_info.py`, or `calcs/` should be converted to `logger.info()` / `logger.warning()`.
- Check for bare `except:` clauses (should be `except Exception:` or more specific).
- Check for any `TODO` or `FIXME` comments that reference a specific version or issue that has since been resolved.

---

## Reporting format

After completing all checks, report in this format:

### Fixed
- [item] — what was wrong, what was changed

### Needs attention (flag for human)
- [item] — what was found, why human judgment is needed (especially anything touching HFSS/Ansys or public API)

### Clean
- [item] — checked, no action needed

### Skipped
- [item] — why (e.g. requires live HFSS session, requires PyPI access, etc.)
