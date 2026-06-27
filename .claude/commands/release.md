# /release

Walk through the full pyEPR release process for the next version. Work through
each step in order and do not skip steps, even if they seem already done.

---

## Step 1 — Confirm readiness

Before touching any version numbers, verify the codebase is clean:

```bash
# All non-HFSS tests pass
pytest -q

# Docs build with zero warnings
rm -rf docs/source/_tutorial_notebooks
cp -r _tutorial_notebooks docs/source/_tutorial_notebooks
cd docs && make html 2>&1 | grep -c WARNING
# Must print 0

# No uncommitted changes
git status
```

If anything fails, stop and fix it first. A release with a broken test suite
or broken docs is worse than no release.

## Step 2 — Determine the new version number

```bash
python -c "import pyEPR; print('Current version:', pyEPR.__version__)"
git tag --list | sort -V | tail -5
```

Decide the new version following these rules:
- **Patch** (`0.9.x → 0.9.x+1`): bug fixes, documentation updates, dependency
  compatibility fixes, CI improvements. No new features, no deprecations.
- **Minor** (`0.9.x → 0.10.0`): new public API (additive), new tutorials,
  significant documentation overhaul. Should not break any existing user code.
- **Major** (`0.9.x → 1.0.0`): reserved for architectural changes that may
  break existing workflows.

When in doubt, use a patch release. It is always safe to release patch versions
frequently.

## Step 3 — Summarise changes since last release

```bash
# Find the commit SHA of the last release tag
git log --oneline $(git tag --list | sort -V | tail -1)..HEAD
```

Group the commits into categories for the release notes:
- New features
- Bug fixes
- Documentation
- CI / packaging
- Dependency compatibility

Draft release notes. Format (markdown, for GitHub Releases body):

```markdown
## What's new in X.Y.Z

### Bug fixes
- ...

### Documentation
- ...

### CI / packaging
- ...

### Dependencies
- ...
```

## Step 4 — Bump the version

There are **two** version fields in `pyEPR/__init__.py` — both must be updated:

```python
# 1. The machine-readable version (line ~93) — what pyproject.toml and pip read:
__version__ = "X.Y.Z"

# 2. The hand-maintained module docstring header (line ~62):
@version: X.Y.Z
```

`pyproject.toml` reads dynamically via `version = {attr = "pyEPR.__version__"}` —
do not edit `pyproject.toml` directly.

Also update `@maintainer` in the docstring if the maintainer list has changed.

Verify the dynamic read works:
```bash
python -c "import pyEPR; print(pyEPR.__version__)"
# Must print the new version
```

Check both fields are in sync:
```bash
grep -n "@version\|__version__" pyEPR/__init__.py
# Both lines must show X.Y.Z
```

## Step 5 — Commit and open a PR

```bash
git checkout -b release/vX.Y.Z
git add pyEPR/__init__.py
git commit -m "bump version to X.Y.Z"
git push -u origin release/vX.Y.Z
```

Open a PR titled `bump version to X.Y.Z`. The PR body should include the
release notes drafted in Step 3 (reviewers need to know what they are approving).

Wait for CI to pass (pylint, pytest, docs build). Do not merge if any CI job fails.

## Step 6 — Merge the PR

After CI passes and any review is complete, merge the PR to master.

Confirm the merge landed:
```bash
git fetch origin master
git log --oneline origin/master -3
python -c "
import subprocess, sys
result = subprocess.run([sys.executable, '-c',
    'import importlib.util; spec = importlib.util.spec_from_file_location(\"pyEPR\", \"pyEPR/__init__.py\"); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); print(m.__version__)'],
    capture_output=True, text=True)
print(result.stdout.strip())
"
```

## Step 7 — Create the GitHub Release

Use the GitHub MCP tool or the GitHub web UI. Key requirements:
- **Tag:** `vX.Y.Z` (with the `v` prefix, e.g. `v0.9.6`)
- **Target:** must point to the master commit that has the bumped `__version__`
- **Release title:** `pyEPR vX.Y.Z`
- **Body:** the release notes from Step 3
- **Not a pre-release** (unless intentionally releasing a beta)

The `publish-to-pypi.yml` GitHub Actions workflow triggers automatically on
`release: created` and publishes to PyPI via OIDC Trusted Publishing. No API
token is needed.

## Step 8 — Verify the PyPI publish

Wait 3–5 minutes, then:

```bash
pip index versions pyEPR-quantum 2>/dev/null | head -3
# or
curl -s https://pypi.org/pypi/pyEPR-quantum/json | python -c "import sys,json; d=json.load(sys.stdin); print(list(d['releases'].keys())[-5:])"
```

The new version must appear. If it does not appear after 10 minutes:
- Check the `publish-to-pypi.yml` Actions run for errors.
- Common issue: the tag was created on the wrong commit (before the version bump).
  If this happened, the wheel was built with the old `__version__`. You must
  yank the release on PyPI and redo from Step 6.

## Step 9 — Post-release check

```bash
# Verify clean install from PyPI
pip install pyEPR-quantum==X.Y.Z
python -c "import pyEPR; print(pyEPR.__version__)"
# Must print X.Y.Z

# Verify the no-HFSS import path still works
python -c "
from pyEPR.solution_types import normalize, DRIVEN_MODAL_NAMES
from pyEPR import ProjectInfo, QuantumAnalysis
from pyEPR.calcs.transmon import transmon_get_spectrum_charge_basis
print('All imports OK')
"
```

If anything fails at this stage, file a hotfix immediately — users are now
installing the broken version.

## Step 10 — Update downstream if needed

If this release contains a breaking change (rare) or a significant new feature
that downstream packages should adopt:
- Check quantum-metal's `pyproject.toml` for their `pyEPR-quantum >= X.Y.Z` pin.
- If their minimum is now outdated, consider opening an issue or PR on quantum-metal
  to update their minimum version. Do not do this silently.
