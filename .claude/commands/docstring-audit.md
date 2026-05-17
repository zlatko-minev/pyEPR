# /docstring-audit

Systematically audit and improve docstrings across the pyEPR codebase. Work
module by module. For each module, check coverage, style, and accuracy —
then fix what you can safely fix (pure text changes, no logic changes).

After completing each module, run the docs build to confirm zero new warnings
were introduced before moving to the next module.

---

## Style standard

All public methods in `core_*.py`, `project_info.py`, and `calcs/` must use
**NumPy-style docstrings**. The structure is:

```python
def example(self, x, y=None):
    """One-line summary, ending with a period.

    Optional extended description. Explain *why* non-obvious constraints exist,
    not *what* the code does (the code already shows that).

    Parameters
    ----------
    x : float
        Description of x, including units if physical (e.g. "Junction
        inductance in henries.").
    y : int, optional
        Description of y. Default is None, which means ...

    Returns
    -------
    result : np.ndarray, shape (n_modes, n_junctions)
        Description of the return value, including shape/dtype when relevant.

    Raises
    ------
    ValueError
        If x is negative.

    Notes
    -----
    Any mathematical context or derivation references go here. Use
    :math:`...` for inline LaTeX and `.. math::` for display equations.

    Examples
    --------
    >>> example(1.0, y=2)
    array([...])
    """
```

---

## RST pitfalls to avoid while writing docstrings

Read `.claude/context/lessons-learned.md` first. The most common failures:

- `|x⟩` bra-ket notation → use `:math:`|x\\rangle``
- `**kwargs` in paragraphs → wrap in `` ``**kwargs`` ``
- `----` separator lines → use NumPy section headers instead
- First line indented relative to `"""` → start on the same line as `"""`
- Duplicate class in both narrative RST and `api/` → use cross-reference in narrative

---

## Module priority order

Work through modules in this order. Higher priority = more user-facing.

### Priority 1 — `calcs/` subpackage

These are the most user-facing for the no-HFSS audience. Pure math, no COM.

```bash
python -m pydocstyle pyEPR/calcs/ --convention=numpy 2>&1 | head -40
```

Check each module:
- `calcs/basic.py` — general EM and qubit parameter formulas
- `calcs/transmon.py` — transmon charge-basis diagonalization
- `calcs/convert.py` — unit conversion helpers
- `calcs/hamiltonian.py` — Hamiltonian matrix construction
- `calcs/back_box_numeric.py` — numerical black-box diagonalization
- `calcs/constants.py` — physical constants (usually just needs a module docstring)

For each function/class: does it have a docstring? Is the Parameters section
complete? Does it describe units? Is there a Returns section?

After editing each file:
```bash
rm -rf docs/source/_tutorial_notebooks
cp -r _tutorial_notebooks docs/source/_tutorial_notebooks
cd docs && make html 2>&1 | grep -c WARNING
# Must still be 0
```

### Priority 2 — `project_info.py`

This is the user's first touch point — they configure `ProjectInfo` before
anything else. Its docstrings are the most important for onboarding.

Check:
- Class docstring explains what `ProjectInfo` is and its role in the pipeline.
- `__init__` documents every parameter including `project_path`, `project_name`,
  `design_name`, and `setup_name`.
- `junctions` attribute: document the dict schema
  (`{name: {'Lj_variable': str, 'rect': str, 'line': str, 'Cj_variable': str}}`).
  Misconfigured junctions are the most common user error.
- `dissipative` attribute: document the dict keys
  (`dielectrics_bulk`, `dielectric_surfaces`, `resistive_surfaces`, `seams`).

### Priority 3 — `core_quantum_analysis.py`

This is the post-processing half that runs without HFSS. Used heavily by
the no-HFSS audience.

Check key methods:
- `analyze_all_variations()` — most important method. Document `cos_trunc`,
  `fock_trunc` parameters with their typical ranges and what each controls.
- `report_results()` — document `swp_variable` and `numeric` parameters.
- `plot_hamiltonian_results()` — document what is plotted.
- `get_Pmj()` — document return shape.

### Priority 4 — `core_distributed_analysis.py`

The HFSS-dependent half. Lower priority for the no-HFSS audience, but
important for HFSS users.

Check key methods:
- `do_EPR_analysis()` — the main entry point. Document what it does, what
  it saves, and where.
- `calc_p_junction()` — document the physical meaning of the return value.

Do not add docstrings to private methods (leading `_`) unless they are
complex enough to warrant it.

### Priority 5 — `toolbox/`

Lower priority. Check for missing module docstrings and class docstrings.

---

## Checking coverage

```bash
# List all public functions/classes missing docstrings
python -c "
import ast, sys
from pathlib import Path

def check_file(path):
    src = path.read_text()
    tree = ast.parse(src)
    missing = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if not node.name.startswith('_'):
                if not (node.body and isinstance(node.body[0], ast.Expr)
                        and isinstance(node.body[0].value, ast.Constant)):
                    missing.append(f'{path}:{node.lineno} {node.name}')
    return missing

for p in sorted(Path('pyEPR').rglob('*.py')):
    if '__pycache__' in str(p):
        continue
    for m in check_file(p):
        print(m)
" 2>/dev/null | grep -v ansys.py | head -50
```

Work through the list starting from `calcs/` and `project_info.py`.

---

## What not to do

- Do not add docstrings to `ansys.py` COM methods without human review.
  The parameter names and types in COM calls are version-sensitive and
  a wrong docstring is worse than no docstring.
- Do not change method signatures while adding docstrings.
- Do not add `Examples` sections that require Ansys HFSS — they will not
  run in doctests and will mislead users.
- Do not add `.. automethod::` directives to narrative RST pages — they
  duplicate API pages and cause duplicate-object warnings.

---

## Reporting

After completing each module, report:
- Module name
- Number of functions/classes updated
- Any warnings that appeared in the docs build (should be zero)
- Any docstrings you skipped and why
