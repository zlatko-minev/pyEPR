# Ecosystem, Adoption, and Developer Relations Context

Understanding *why* pyEPR exists, who uses it, and what its place in the
broader quantum computing stack is — this is necessary context for making
good decisions about API design, documentation, and release timing.

---

## What pyEPR is and is not

**pyEPR is** a post-processing and quantization library. It takes eigenmode
field simulation results (from Ansys HFSS or any compatible solver) and
produces quantum Hamiltonian parameters: mode frequencies, anharmonicities,
dispersive shifts (χ matrix), and zero-point fluctuations.

**pyEPR is not:**
- A geometry or layout tool — that is quantum-metal (formerly qiskit-metal) or KLayout.
- A general AEDT scripting library — that is PyAEDT (though pyEPR optionally uses PyAEDT as a transport layer via `ansys_pyaedt`).
- A circuit simulator — it post-processes distributed-field solutions.
- An alternative to QuTiP — it uses QuTiP internally for diagonalization.

### Two HFSS transport backends

pyEPR has two ways to talk to HFSS:

| Backend | Module | Transport | Platform | Install |
|---------|--------|-----------|----------|---------|
| Classic | `pyEPR.ansys` / `DistributedAnalysis` | COM / pywin32 | Windows only | `pip install pyEPR-quantum` |
| PyAEDT | `pyEPR.ansys_pyaedt` / `PyaedtDistributedAnalysis` | gRPC | Linux, macOS, Windows | `pip install "pyEPR-quantum[pyaedt]"` |

The PyAEDT backend was added in 0.9.6+. It is fully additive — the COM backend is
unchanged. The physics (participation formula, diagonalization) is identical between
the two; only the transport layer differs. Use the PyAEDT backend when:
- You are on Linux or macOS
- You are hitting COM stale-session or project-locked errors
- You want to use Ansys's officially maintained API rather than the COM interface

The key technical detail: reading scalar results from the HFSS field calculator over
gRPC requires `CalculatorWrite` (write to a `.fld` file, read last line) rather than
`ClcEval`/`GetTopEntryValue`, which is the stateful round-trip that does not survive gRPC.
See `.claude/context/lessons-learned.md` for the full explanation.

The distinction matters because feature requests often conflate these roles.
When someone asks pyEPR to "draw a qubit" or "run a SPICE simulation", the
answer is "that is not this tool's scope."

---

## The user spectrum

**Group 1: Experimental hardware groups (primary)**
These are the core users. They have Ansys HFSS licences, they run eigenmode
simulations of 3D cavities or planar chips, and they use pyEPR to extract
Hamiltonian parameters. They are physicists first, programmers second.
Their pain points: COM connection stability, AEDT version compatibility,
clear error messages when junction parameters are misconfigured.

**Group 2: Theorists and students (high adoption potential)**
No Ansys licence. They want to run the numerical EPR workflow (Tutorial 6)
to understand the method, or they load a pre-computed HDF5 file from a
collaborator. This group is the largest potential audience and the most
underserved. They are served by: Tutorial 6, `QuantumAnalysis` loaded from
HDF5, `calcs/` subpackage, Binder links, and headless docs.

**Group 3: Chip design automation users (strategic)**
Using quantum-metal (formerly qiskit-metal) for layout design and relying
on pyEPR for the EPR analysis step. These users rarely interact with pyEPR
directly — quantum-metal calls it internally. They care about: import
stability, no Windows-only deps on the analysis path, correct pyEPR version
pinning in quantum-metal's requirements.

---

## The quantum-metal relationship

quantum-metal (PyPI: `quantum-metal`, formerly `qiskit-metal`) is the main
downstream consumer of pyEPR. It:
- Declares `pyEPR-quantum >= 0.9.5` as a dependency
- Imports `pyEPR.solution_types` (must be importable without COM stack)
- Imports `pyEPR.ProjectInfo`, `pyEPR.DistributedAnalysis`, `pyEPR.QuantumAnalysis`
- Uses the deprecated aliases `pyEPR_HFSSAnalysis` and `pyEPR_Analysis`

**Constraints this imposes on pyEPR:**
- Do not remove deprecated aliases. quantum-metal imports them.
- Do not add Windows-only imports to `solution_types.py` or `calcs/`. quantum-metal
  runs on Linux/macOS and these modules must import cleanly.
- Treat any change to the public API of `ProjectInfo`, `DistributedAnalysis`, or
  `QuantumAnalysis` as potentially breaking for quantum-metal users until
  quantum-metal ships a new release.
- When pyEPR releases a new version, quantum-metal's `pyEPR-quantum >= X.Y.Z`
  pin means users on older pyEPR will still run. Breaking changes in a new
  pyEPR minor release can be masked for a long time if quantum-metal's minimum
  is not bumped. Be conservative.

---

## The no-HFSS path — the most important adoption lever

Most researchers who hear about pyEPR are blocked from trying it because they
assume they need Ansys HFSS. They do not, for a large class of use cases.

The no-HFSS path:
1. **Tutorial 6** — supply frequencies, junction inductances, and φ_zpf manually;
   get the full χ matrix back. No EM solver at all.
2. **`QuantumAnalysis` from HDF5** — if a collaborator runs the HFSS extraction
   and saves the HDF5 file, the quantum analysis runs anywhere.
3. **`calcs/` subpackage** — direct formula access. Transmon spectrum, EPR formulas,
   unit conversions. Pure Python/numpy/scipy.

These paths must be the most polished, best-documented, and most reliably
tested parts of the codebase. They are the onramp. If Tutorial 6 has a broken
import or a unit bug, the user's first impression of pyEPR is broken.

**Binder:** Tutorial 6 must have a working Binder badge. Binder allows a new
user to run the notebook in the browser with zero install. This is the lowest
possible friction onramp. Keep the badge URL correct and test it periodically;
Binder builds can fail if the `requirements.txt` or `environment.yml` goes stale.

---

## The documentation philosophy

The docs use the **PyData Sphinx Theme** — the same theme as NumPy, SciPy,
pandas, and QuTiP. This is intentional: it signals to researchers that
pyEPR is part of the scientific Python ecosystem they already know. It also
provides a modern, navigable layout without custom CSS work.

**sphinx-design** provides the visual components: feature cards on the
landing page, install tabs, HFSS/no-HFSS badges on tutorials. These make
the docs look maintained and signal capability quickly to a new visitor.

**myst-nb** handles notebook rendering. Notebooks are executed once and
saved with outputs; myst-nb renders the saved outputs without re-executing
(`nb_execution_mode = "off"`). This means: the docs always reflect the
state of the notebook at commit time, and the build does not require Ansys.

**Zero-warning build is non-negotiable.** A docs build with warnings looks
unmaintained and erodes trust. Every warning is a signal to a potential user
that the project is not actively cared for.

---

## Docstring and API documentation philosophy

**NumPy-style docstrings** — the standard across NumPy, SciPy, pandas,
QuTiP, and scikit-learn. Hardware physicists who are already familiar with
the scientific Python stack will find this immediately readable.

Docstrings should explain:
- *What* the method does (one-line summary)
- *Why* non-obvious parameters have their specific constraints or units
- *What* is returned and in what units/shape

Docstrings should not explain *how* the method is implemented — that
belongs in comments in the code, and only when the implementation is
non-obvious.

**API pages** are auto-generated from docstrings via `automodule` in `docs/source/api/`.
The narrative `key_classes_reference.rst` page should contain *explanation and
context*, not duplicated API tables. If a class appears in both the narrative
page and the API page, use a cross-reference link (`:class:`...``) in the
narrative page, not a redundant `.. autoclass::`.

---

## Release and versioning philosophy

pyEPR follows semantic versioning loosely:
- Patch releases (`0.9.x`) — bug fixes, documentation, dependency compatibility.
  Should not break any user code.
- Minor releases (`0.x.0`) — new features, potentially deprecating old patterns.
  Deprecated items stay for at least one minor version.
- Major releases (`x.0.0`) — reserved for significant architectural changes.
  The current codebase has not had one.

**Release timing:** pyEPR does not have a fixed release cadence. Releases are
driven by meaningful accumulation of changes: a significant new feature,
a batch of bug fixes, or a documentation overhaul. Small releases are fine;
releasing with only a version bump to pick up trivial changes is not useful.

**PyPI publish is automated** via OIDC Trusted Publishing. No API token is
needed. Create the GitHub Release with the correct tag; the workflow fires
automatically. Verify on PyPI after ~5 minutes.

---

## What the community cares about

Based on issues, questions, and usage patterns:

1. **AEDT version compatibility** — the most common issue. New AEDT releases
   change API strings, design creation behaviour, and scripting call signatures.
   Users upgrade AEDT and pyEPR breaks. This is the highest-friction pain point.

2. **The junction setup** — misconfigured junction rectangles and polylines are
   the most common user error. Clear error messages when junction parameters are
   wrong save enormous amounts of support time.

3. **The no-HFSS workflow** — constantly requested by students and theorists.
   Tutorial 6 is the answer; keep it working.

4. **Import stability** — researchers copy-paste pyEPR code into their own
   analysis scripts and expect it to continue working across pyEPR updates.
   The deprecated aliases (`pyEPR_Analysis`, etc.) exist for exactly this reason.

5. **Documentation** — the most consistent feedback on any scientific Python
   project is "better docs." The PyData theme migration, the tutorial gallery,
   and the zero-warning build are direct responses to this.
