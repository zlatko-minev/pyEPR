# Project State — pyEPR

This file is a living snapshot of the project's current state, key decisions,
and open items. Update it when making significant changes so future agents can
orient quickly without reading the full commit history.

---

## Current version

**1.0.0** — released 2026-06-27.

First major version. Milestone: cross-platform HFSS access via the PyAEDT gRPC
backend. All 0.9.x COM workflows are unchanged and fully backwards compatible.

---

## Maintainers

| Name | GitHub | Role |
|------|--------|------|
| Zlatko Minev | [@zlatko-minev](https://github.com/zlatko-minev) | Author, lead maintainer |
| Joey Yaker | [@joeyyaker](https://github.com/joeyyaker) | Maintainer — PyAEDT gRPC backend |

---

## Architecture snapshot

Two HFSS transport backends coexist:

| Backend | Module | Class | Install | Platform |
|---------|--------|-------|---------|----------|
| COM (classic) | `pyEPR.ansys` | `DistributedAnalysis` | `pip install pyEPR-quantum` | Windows only |
| gRPC (new) | `pyEPR.ansys_pyaedt` | `PyaedtDistributedAnalysis` | `pip install "pyEPR-quantum[pyaedt]"` | Linux · macOS · Windows |

The gRPC backend was contributed by Joey Yaker (PR #207, merged 2026-06-27).
Physics and results are identical — validated digit-for-digit (`p_mj = 0.9755`
on a demo transmon with `pyaedt 1.1.0` + AEDT 2026.1).

Key technical detail: `CalculatorWrite` (write to `.fld` file, read last line)
must replace `ClcEval`/`GetTopEntryValue` for results to survive gRPC. See
`lessons-learned.md` for the full explanation.

---

## Key downstream dependency

**quantum-metal** (PyPI: `quantum-metal`, formerly `qiskit-metal`) declares
`pyepr-quantum>=0.9.5` with **no upper bound** in its `ansys` and `full` extras.
As of 2026-06-27, quantum-metal is at v0.7.4. The `>=0.9.5` pin is satisfied by
1.0.0 — no quantum-metal update required.

---

## Recent significant changes (since 0.9.5)

| PR | What | Merged |
|----|------|--------|
| #199–203 | PyData Sphinx theme migration, myst-nb, tutorials gallery, zero-warning build | 2026-05-xx |
| #204 | Version bump 0.9.5 → 0.9.6 | 2026-05-17 |
| #207 | PyAEDT gRPC backend (`pyEPR.ansys_pyaedt`) — Joey Yaker | 2026-06-27 |
| #208 | Companion docs: dual-backend diagram, README feature section, Joey's credit | 2026-06-27 |
| #209 | Version bump 0.9.6 → 1.0.0 | 2026-06-27 |

---

## CI jobs

| Job | What it tests |
|-----|---------------|
| `test` (4 matrix: py3.10+3.12 × ubuntu+macos) | Full non-HFSS test suite |
| `test_pyaedt` | Installs `[pyaedt]` extra, runs offline PyAEDT backend tests |
| `test_docs` | Sphinx build, zero-warning gate |
| `pylint` | `--errors-only`, includes `pyEPR.ansys_pyaedt` |
| `greeting` | First-time contributor welcome bot |
| `publish-to-pypi` | Triggers on `release: created` via OIDC Trusted Publishing |

HFSS live tests (`@pytest.mark.hfss`) are never run in CI — they require a
physical Ansys licence. The test_pyaedt CI job covers only the offline-testable
functions: `compute_p_mj`, `_parse_henries`, `_owning_session_from_lock`.

---

## Open items (as of 2026-06-27)

- [ ] **GitHub Release v1.0.0** — create after PR #209 merges to trigger PyPI publish.
      Tag: `v1.0.0`, target: master post-merge commit.
- [ ] **Fork CI auto-run** — confirm `github.com/zlatko-minev/pyEPR/settings/actions`
      is set to "Run workflows automatically" for outside collaborators.
- [ ] **Issue #206** — close after v1.0.0 PyPI publish is confirmed.
- [ ] **quantum-metal** — no action needed for 1.0.0 (no breaking changes, no
      upper-bound pin). Monitor if a 2.0.0 ever becomes necessary.

---

## Release checklist (condensed)

1. Bump `__version__` in `pyEPR/__init__.py` (the `__version__ = "X.Y.Z"` line).
2. Also update `@version: X.Y.Z` in the module docstring (~line 62). **Both must match.**
3. Update `CHANGELOG.md`: rename `[Unreleased]` to `[X.Y.Z] — YYYY-MM-DD`.
4. Open PR `release/vX.Y.Z` → master, wait for CI green, merge.
5. Create GitHub Release with tag `vX.Y.Z` pointing at the post-merge master commit.
6. Verify PyPI publish (~5 min): `pip index versions pyEPR-quantum | head -1`.

Full procedure: see `.claude/commands/release.md`.
