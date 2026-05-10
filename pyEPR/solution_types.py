"""
pyEPR.solution_types
====================
Canonical names and alias sets for Ansys HFSS / Q3D solution types.

Background
----------
Ansys AEDT 2021.2 renamed the strings returned by ``GetSolutionType()``:

    Pre-2021.2 (canonical)      AEDT 2021.2+
    ─────────────────────────────────────────────────────────
    "DrivenModal"               "HFSS Modal Network"
    "DrivenTerminal"            "HFSS Terminal Network"
    —                           "HFSS Hybrid Modal Network"
    —                           "HFSS Hybrid Terminal Network"

pyEPR normalises these to the pre-2021.2 canonical strings in
``HfssDesign.__init__`` (PR #176) using substring matching (same approach
as PyAEDT), so all downstream pyEPR code always receives ``"DrivenModal"``,
``"DrivenTerminal"``, etc.

This module exposes the alias sets and helper functions publicly so that:

1. pyEPR's own ``get_setup()`` / ``connect_setup()`` can import the sets
   instead of duplicating them as inline tuples.
2. qiskit-metal (and other consumers) can import from here rather than
   maintaining a parallel copy of the alias sets.

Note for qiskit-metal maintainers
----------------------------------
Because pyEPR normalises at read time, ``design.solution_type`` will
always yield the *canonical* name (e.g. ``"DrivenModal"``), never the
raw AEDT-2021.2 string.  Metal's own frozenset layer is still correct as
defence-in-depth for code paths that call ``GetSolutionType()`` directly
without going through pyEPR's ``HfssDesign`` wrapper.  See pyEPR PR #176
and the qiskit-metal branch ``claude/sync-qtip-pyepr-versions-tEJd1``.
"""

# ── Canonical names (pre-AEDT-2021.2; used as the normalised output) ──────────

EIGENMODE: str = "Eigenmode"
DRIVEN_MODAL: str = "DrivenModal"
DRIVEN_TERMINAL: str = "DrivenTerminal"
Q3D: str = "Q3D"

# ── Alias sets ─────────────────────────────────────────────────────────────────
# Each set lists *every* string GetSolutionType() is known to emit for that
# family.  When HFSS introduces a new alias, add it here; normalize() and
# the is_*() predicates are automatically updated.

EIGENMODE_NAMES: frozenset = frozenset({
    "Eigenmode",
})

DRIVEN_MODAL_NAMES: frozenset = frozenset({
    "DrivenModal",                # pre-AEDT 2021.2 canonical
    "HFSS Modal Network",         # AEDT 2021.2+
    "HFSS Hybrid Modal Network",  # hybrid variant (some builds)
})

DRIVEN_TERMINAL_NAMES: frozenset = frozenset({
    "DrivenTerminal",                # pre-AEDT 2021.2 canonical
    "HFSS Terminal Network",         # AEDT 2021.2+
    "HFSS Hybrid Terminal Network",  # hybrid variant (some builds)
})

Q3D_NAMES: frozenset = frozenset({
    "Q3D",
})

# ── Predicates ─────────────────────────────────────────────────────────────────


def is_eigenmode(s: str) -> bool:
    """Return True if *s* is any known eigenmode solution-type string."""
    return s in EIGENMODE_NAMES


def is_drivenmodal(s: str) -> bool:
    """Return True if *s* is any known driven-modal solution-type string."""
    return s in DRIVEN_MODAL_NAMES


def is_driventerminal(s: str) -> bool:
    """Return True if *s* is any known driven-terminal solution-type string."""
    return s in DRIVEN_TERMINAL_NAMES


def is_q3d(s: str) -> bool:
    """Return True if *s* is any known Q3D solution-type string."""
    return s in Q3D_NAMES


def canonical_kind(s: str):
    """Map *s* to one of ``'eigenmode'``, ``'drivenmodal'``, ``'driventerminal'``,
    ``'q3d'``, or ``None`` (unrecognised).

    Comparisons are case-sensitive; HFSS always emits CamelCase.
    """
    if s in EIGENMODE_NAMES:
        return "eigenmode"
    if s in DRIVEN_MODAL_NAMES:
        return "drivenmodal"
    if s in DRIVEN_TERMINAL_NAMES:
        return "driventerminal"
    if s in Q3D_NAMES:
        return "q3d"
    return None


def normalize(raw: str) -> str:
    """Convert a raw ``GetSolutionType()`` string to the canonical pre-2021.2 form.

    Uses substring matching so that future HFSS aliases that still contain
    ``"Modal"`` / ``"Terminal"`` / ``"Eigenmode"`` are handled automatically
    (same strategy as PyAEDT).

    Unknown types (Q3D, SBR+, Transient, custom) are returned unchanged.
    Comparisons are case-sensitive.
    """
    if "Modal" in raw:
        return DRIVEN_MODAL
    if "Terminal" in raw:
        return DRIVEN_TERMINAL
    if "Eigenmode" in raw:
        return EIGENMODE
    return raw
