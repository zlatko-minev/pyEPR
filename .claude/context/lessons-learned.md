# Lessons Learned — Hard-Won Maintenance Knowledge

Every item here caused a real build failure, a silent test pass, a broken
PyPI page, or a regression for downstream users. Read this before touching
docs, CI, or anything related to Ansys version compatibility.

---

## Sphinx / Documentation

### Sphinx won't follow symlinks outside its source root

`docs/source/_tutorial_notebooks` is a git symlink pointing to
`../../_tutorial_notebooks/`. Sphinx silently ignores it — no warning,
no error, just missing content and toctree "nonexisting document" errors.

**Fix:** Copy before every build. This step is required in three places:
- Local `make html` (do it manually)
- `.readthedocs.yml` under `jobs.pre_build`
- `.github/workflows/ci.yaml` in the `test_docs` job

If you add a new CI job that builds docs, add the copy step. Forgetting
it produces a false-green CI run (zero errors) because toctree entries
for missing files are silently dropped, not errored.

### sphinx-design `link-type: doc` silently strips spaces and periods

When using `:link: _tutorial_notebooks/Tutorial 3.  toolbox_circuits`
with `link-type: doc`, sphinx-design normalises the path by stripping
spaces and periods, producing `Tutorial3toolbox_circuits` — a path that
does not exist.

**Fix:** Use `link-type: url` with a relative `.html` path. Browsers
percent-encode spaces in hrefs automatically:
```rst
:link: _tutorial_notebooks/Tutorial 3.  toolbox_circuits.html
:link-type: url
```

### Duplicate object description warnings

If a class appears in both `key_classes_reference.rst` (via `.. autoclass::`)
and in `api/` (via `.. automodule::`), Sphinx emits a duplicate-object warning
that counts against the zero-warning gate.

**Fix:** Remove `.. autoclass::` from narrative RST files. Replace with a
plain cross-reference link: `:class:`pyEPR.core_quantum_analysis.QuantumAnalysis``.
The `api/` automodule page handles the full API documentation.

If removing autoclass is not possible, add `:no-index:` to the second
occurrence, but prefer the cross-reference approach.

### RST pipe characters in docstrings

`|x⟩` bra-ket notation with a pipe character triggers RST substitution
reference parsing. The warning is `Substitution_reference "x⟩" not found`
and the rendered output is broken.

**Fix:** Use `:math:`|x\\rangle`` instead of plain Unicode.

### `**kwargs` in docstring paragraphs

Bare `**` inside a paragraph triggers RST strong emphasis markup.
`**kwargs` renders as bold `kwargs` and may cause "Inline strong start-string
without end-string" warnings.

**Fix:** Wrap in double backticks: `` ``**kwargs`` ``.

### Separator lines in docstrings

A line of `----` at the start of a docstring section causes:
`WARNING: Transition must be child of a section or document`

**Fix:** Use a proper RST section header (underlined with `~` or `-`),
or remove the separator entirely. The NumPy docstring convention uses
section headers (`Parameters`, `Returns`, etc.), not horizontal rules.

### Indented first line of docstring

A docstring whose first line is indented relative to the opening `"""`
creates a block-quote RST node. Sphinx renders it oddly and may warn.

**Fix:** Start the first line immediately after `"""`, at the same
indentation level as the opening quotes.

### Wrong directive names

- `.. tabs::` / `.. tab::` — from `sphinx-tabs`, which is **not installed**.
  Use `.. tab-set::` / `.. tab-item::` from `sphinx-design` instead.
- `.. nbsphinx::` — not installed. Notebooks are handled by `myst-nb`.
- `..codeblock python` (no space) — silently ignored, not rendered.
  Must be `.. code-block:: python` with a space before the language.

### `suppress_warnings` in `conf.py`

Some warnings are unavoidable (myst-nb header level conventions, notebook
syntax highlighting fallbacks). These are suppressed in `conf.py`:
```python
suppress_warnings = ["myst.header", "misc.highlighting_failure", "ref.python"]
```
Do not add new entries here without first investigating whether the warning
can be fixed at the source. This list should stay short.

---

## Version management

### Two version fields in `__init__.py` — both must be updated

`pyEPR/__init__.py` contains two separate version references:

1. `__version__ = "X.Y.Z"` (~line 93) — the machine-readable string that
   `pyproject.toml` reads via `{attr = "pyEPR.__version__"}`. This is what pip,
   PyPI, and `import pyEPR; pyEPR.__version__` see.

2. `@version: X.Y.Z` (~line 62) — a hand-maintained field in the module-level
   docstring header block (`@author`, `@site`, `@license`, `@version`, …).

When bumping the version, **both must be updated**. The `@version` field is easy
to miss because it is inside a string literal and `grep __version__` won't find it.
Always run:

```bash
grep -n "@version\|__version__" pyEPR/__init__.py
```

and confirm both lines show the new version before committing.

The `@maintainer` field in the same docstring block should also be kept current
when the maintainer list changes.

---

## Python / Dependencies

### qutip 5.x — `ket.dag() * ket` returns a scalar, not a 1×1 Qobj

In qutip 4.x, `ket.dag() * ket` returned a 1×1 `Qobj` with a `.norm()` method.
In qutip 5.x, the same expression returns a complex scalar.

Code that calls `.norm()` or `.tr()` on the result will raise `AttributeError`
on qutip 5.x.

**Fix:**
```python
inner = ket.dag() * ket
if hasattr(inner, "norm"):   # qutip 4.x
    val = inner.norm()
else:                        # qutip 5.x — already a scalar
    val = abs(inner)
```

### pandas FutureWarnings

pandas regularly deprecates indexing patterns. The most common offenders:
- `df.loc[:, col]` vs `df[col]` — silent in old versions, warning in newer
- `DataFrame.append()` — removed in pandas 2.0; use `pd.concat()`
- Chained assignment (`df[col][mask] = val`) — use `.loc` instead

Run `python -W error::FutureWarning -m pytest tests/` to surface these.

### numpy 2.x — `np.bool`, `np.int`, `np.float` removed

These aliases were deprecated in numpy 1.20 and removed in 2.0.
Any code using `np.bool`, `np.int`, or `np.float` directly will raise
`AttributeError` on numpy 2.x.

The current `pyproject.toml` pins numpy to a range that may not include 2.x.
Do not lift a numpy upper-bound pin without auditing all such usages first.

### Windows-only imports must stay confined to `ansys.py`

`win32com`, `pythoncom`, and `pywintypes` are Windows-only. If any of these
appear at the top level of any module outside `ansys.py`, the package will
fail to import on Linux and macOS — silently, with an `ImportError`.

`solution_types.py` and `calcs/` must never import from `ansys.py` or
from any Windows-only library. This is a hard constraint: downstream
packages (quantum-metal) run on Linux and import these modules.

---

## CI / Release

### CI docs job needs the notebook copy step

The `test_docs` CI job initially did not have the notebook copy step, so
it ran `make html` against a source tree where `_tutorial_notebooks/` was
an unresolved symlink. The job showed zero errors (Sphinx silently dropped
missing toctree entries) but the built docs were missing all six tutorials.

**Fix:** The copy step must be in the CI YAML, not only in `.readthedocs.yml`.

### Tag must point to the commit with the bumped version

The PyPI publish workflow reads `__version__` from `pyEPR/__init__.py`
at the time the wheel is built. If the GitHub Release tag points to a
commit where `__version__` has not yet been bumped (e.g. the merge commit
of the feature work, not the version-bump commit), the wheel will carry
the old version number.

**Correct order:**
1. Merge all feature PRs
2. Open and merge a version-bump PR (change only `__version__`)
3. Create the GitHub Release tag pointing to the post-merge master commit

### Dead badge image URLs

- `frapsoft.com` (Open Source Love badges) — dead, returns 404
- `rawgit.com/sindresorhus/awesome` — dead CDN, returns 404

Both appeared in `docs/source/about.rst`. Replace with sphinx-design
inline badges (`:bdg-success:`, `:bdg-info:`, etc.) which have no
external runtime dependency.

---

## gRPC field calculator — `CalculatorWrite` vs `ClcEval`

When driving HFSS via PyAEDT over gRPC (the `ansys_pyaedt` backend),
most field-calculator operations work fine: `ClcMaterial`, `EnterVol`,
`EnterLine`, `Integrate`, and `Solutions.EditSources` all run over gRPC.

The one operation that does **not** survive gRPC is the stateful read-back
round-trip `ClcEval` + `GetTopEntryValue`. This is what COM uses to
pull a scalar result off the calculator stack.

**Fix:** Use `CalculatorWrite` instead — write the result to a `.fld` file
in the working directory, then read the last line. This is what
`_GrpcFieldCalc._evaluate()` in `pyEPR/ansys_pyaedt.py` does, and it was
validated digit-for-digit against the COM path (`p_mj = 0.9755` on a demo
transmon). Any future work that needs a scalar from the HFSS field calculator
over gRPC must use `CalculatorWrite`, not `ClcEval`.

---

## Ansys AEDT version compatibility

### AEDT 2024.1+ silently creates Hybrid Modal Network

`InsertDesign("HFSS", name, "DrivenModal", "")` creates an
`HFSS Hybrid Modal Network` design on AEDT 2024.1+, not a plain
`DrivenModal`. The string `"DrivenModal"` is accepted without error
but the design type is wrong.

**Fix:** Call `SetSolutionType` immediately after `InsertDesign` on
AEDT >= 2024.1. The helpers `project.new_dm_design()` and
`project.new_dt_design()` already do this correctly.

### `GetSolutionType` returns new strings on AEDT 2021.2+

AEDT 2021.2 introduced renamed solution-type strings
(`"HFSS Modal Network"` instead of `"DrivenModal"`). New aliases are
added with each AEDT release; the old ones are not removed.

**Fix:** Always pass `GetSolutionType()` output through `normalize()`
from `pyEPR.solution_types` before any string comparison. This is
already done in `HfssDesign.__init__`. Never add raw string comparisons
elsewhere.

---

## PyPI / README rendering

### Relative image paths do not render on PyPI

`README.md` is rendered by PyPI's pipeline. Relative paths like
`imgs/xmon-example.gif` resolve against the local filesystem during
development but produce broken images on PyPI.

**Fix:** Use absolute `https://raw.githubusercontent.com/zlatko-minev/pyEPR/master/imgs/...`
URLs for every image in `README.md`. GitHub's raw CDN is reliable and
renders correctly on both GitHub and PyPI.
