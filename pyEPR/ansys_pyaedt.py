"""
PyAEDT-based HFSS interface for pyEPR — a modern alternative to :mod:`pyEPR.ansys`.

pyEPR drives HFSS through the hand-rolled COM wrapper in :mod:`pyEPR.ansys`. This
module performs the *same* Energy-Participation-Ratio field extraction through
**PyAEDT** (``ansys.aedt.core``), Ansys's official, maintained Python API —
**entirely over gRPC, with no COM at all**:

* project / design / variation / frequency handling over PyAEDT's gRPC transport
  (``Hfss(...)``), including *attach-to-a-running-session*;
* eigenmode normalization via ``Solutions.EditSources`` (gRPC);
* the junction **line** voltage and the material-weighted **volume** energy
  integral ``Re∫ E·ε·E* dV`` over the gRPC field calculator
  (``hfss.ofieldsreporter``), read back with ``CalculatorWrite`` to a file.

The one subtlety that makes this work COM-free: results are read back with
``CalculatorWrite`` (write to a ``.fld`` file), **not** ``ClcEval`` /
``GetTopEntryValue`` — that stateful read-back round-trip is what does not survive
gRPC. ``ClcMaterial`` (the permittivity multiply), ``EnterVol``, ``EnterLine`` and
``Integrate`` all run fine over gRPC. Validated digit-for-digit against the COM
path (``p_mj = 0.9755`` on the demo transmon).

The extracted participations ``p_mj``, signs ``s_mj``, eigenfrequencies and
``L_j`` feed pyEPR's own physics unchanged
(:meth:`pyEPR.calcs.basic.CalcsBasic.epr_to_zpf`,
:func:`pyEPR.calcs.back_box_numeric.epr_numerical_diagonalization`,
:class:`pyEPR.QuantumAnalysis`).

Junctions are declared exactly as in pyEPR — via ``ProjectInfo.junctions`` — using
the keys ``Lj_variable`` (or ``Lj_henries``), ``line``, and an optional ``rect``
or ``solid``.

Requires
--------
* ``pyaedt`` (PyAEDT; provides the ``ansys.aedt.core`` namespace). No
  ``pywin32``/COM needed.

Example
-------
>>> import pyEPR as epr
>>> from pyEPR.ansys_pyaedt import PyaedtDistributedAnalysis
>>> pinfo = epr.ProjectInfo(project_path=r"C:/proj", project_name="Demo",
...                         design_name="Transmon", setup_name="Setup1",
...                         do_connect=False)               # PyAEDT does the connecting
>>> pinfo.junctions['j1'] = {'Lj_variable': 'Lj', 'line': 'jj_line', 'rect': 'jj_rect'}
>>> eprd = PyaedtDistributedAnalysis(pinfo, aedt_version="2026.1")
>>> eprd.do_EPR_analysis()
>>> f_ND, chi_ND = eprd.analyze()        # uses pyEPR's epr_numerical_diagonalization
>>> print(chi_ND.diagonal())             # anharmonicities (MHz)
"""
from __future__ import annotations

import math
import os
import tempfile
from typing import Dict, List, Optional, Tuple

import numpy as np

__all__ = ["PyaedtDistributedAnalysis", "compute_p_mj"]


# ---------------------------------------------------------------------------
#  PyAEDT imported lazily so `import pyEPR` never hard-requires it.
# ---------------------------------------------------------------------------
def _require_pyaedt():
    try:
        from ansys.aedt.core import Hfss  # noqa: F401  # pylint: disable=import-error
        from ansys.aedt.core.generic.general_methods import active_sessions  # noqa: F401  # pylint: disable=import-error
        return Hfss, active_sessions
    except Exception as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "PyAEDT is required for the ansys_pyaedt backend. "
            "Install it with `pip install pyaedt`."
        ) from exc


# ---------------------------------------------------------------------------
#  Unit parsing  (e.g. "12.9201nH" -> 1.29201e-08)
# ---------------------------------------------------------------------------
_L_UNITS = [("nh", 1e-9), ("uh", 1e-6), ("mh", 1e-3), ("kh", 1e3), ("h", 1.0)]


def _parse_henries(expr: str) -> float:
    """Parse an inductance string with units into SI Henries."""
    s = str(expr).strip().lower()
    for unit, scale in _L_UNITS:
        if s.endswith(unit):
            return float(s[: -len(unit)]) * scale
    return float(s)  # already unit-less / SI


# ---------------------------------------------------------------------------
#  p_mj — pyEPR's `line_voltage` recipe (pure arithmetic, no HFSS)
# ---------------------------------------------------------------------------
def compute_p_mj(
    *,
    V_peak: float,
    freq_hz: float,
    Lj_henries: float,
    U_E_raw: float,
    Cj_farads: float = 0.0,
    sum_Ucap_other: float = 0.0,
) -> Tuple[float, int]:
    """Junction participation ``p_mj`` and sign from a line voltage.

    Mirrors :meth:`pyEPR.DistributedAnalysis.calc_p_junction` (``line_voltage``)::

        ω      = 2π f
        I_peak = V_peak / (ω · L_J)
        U_J_ind = ½ · L_J · I_peak²
        U_J_cap = ½ · C_J · V_peak²
        U_norm  = U_E_raw/2 + Σ_j U_J_cap
        p_mj    = U_J_ind / U_norm,   s_mj = +1 if V_peak > 0 else -1

    ``U_E_raw`` is the raw electric integral ``Re∫E·εE* dV`` for the whole mode
    (= 4·U_E), exactly what pyEPR's ``calc_energy_electric`` returns and halves.

    Returns ``(p_mj, sign)``.
    """
    if Lj_henries <= 0 or freq_hz <= 0:
        raise ValueError("Lj_henries and freq_hz must be > 0")
    omega = 2.0 * math.pi * freq_hz
    I_peak = V_peak / (omega * Lj_henries)
    U_J_ind = 0.5 * Lj_henries * I_peak ** 2
    U_J_cap = 0.5 * Cj_farads * V_peak ** 2
    U_norm = (U_E_raw / 2.0) + U_J_cap + sum_Ucap_other
    if U_norm <= 0:
        raise ValueError(
            f"non-positive U_norm ({U_norm}); the mode has no electric energy "
            f"(U_E_raw={U_E_raw}). Is the mode solved and selected?"
        )
    return U_J_ind / U_norm, (1 if V_peak > 0 else -1)


# ---------------------------------------------------------------------------
#  gRPC field calculator (no COM)
# ---------------------------------------------------------------------------
class _GrpcFieldCalc:
    """Field calculator + eigenmode selection over PyAEDT's **gRPC** handle.

    No COM. Two design points make this work where the naive approach fails:

    1. **Normalize the eigenmode with ``EditSources``** (the ``Solutions`` module
       call, which *does* run over gRPC) before reading any field. This gives the
       line voltage and the energy a single, self-consistent field scaling — the
       coupling that otherwise makes ``p_mj`` come out wrong.
    2. **Read results back with ``CalculatorWrite``** (write to a ``.fld`` file),
       NOT ``ClcEval``/``GetTopEntryValue`` — that stateful read-back is the one
       thing that doesn't survive gRPC. ``ClcMaterial`` / ``EnterVol`` /
       ``EnterLine`` / ``Integrate`` themselves are all fine over gRPC.
    """

    def __init__(self, hfss, setup_solution: str):
        self._ofr = hfss.ofieldsreporter
        self._sol = hfss.odesign.GetModule("Solutions")
        self._setup = setup_solution
        self._wdir = getattr(hfss, "working_directory", None) or tempfile.gettempdir()

    def set_mode(self, mode_index_0based: int, phase: float = 0.0) -> None:
        """Normalize fields to one eigenmode via EditSources (over gRPC)."""
        n = mode_index_0based + 1            # HFSS is 1-based
        n_modes = max(10, n)
        mags = ["1" if i + 1 == n else "0" for i in range(n_modes)]
        phs = [str(phase) if i + 1 == n else "0" for i in range(n_modes)]
        self._sol.EditSources(
            [
                ["FieldType:=", "EigenPeakElectricField"],
                ["Name:=", "Modes", "Magnitudes:=", mags, "Phases:=", phs],
            ]
        )

    def _evaluate(self, build) -> float:
        """Build a calculator stack and read it back over gRPC via CalculatorWrite."""
        ofr = self._ofr
        ofr.CalcStack("clear")
        build(ofr)
        path = os.path.join(self._wdir, "epr_calc.fld")
        if os.path.isfile(path):
            os.remove(path)
        ofr.CalculatorWrite(path, ["Solution:=", self._setup], ["Phase:=", "0deg"])
        with open(path) as fh:
            val = float(fh.readlines()[-1].strip())
        try:
            os.remove(path)
        except OSError:
            pass
        ofr.CalcStack("clear")
        return val

    def energy_electric(self, obj: str = "AllObjects") -> float:
        """Raw electric integral ``Re∫ E·ε·E* dV`` (= 4·U_E) over ``obj``, over gRPC."""
        def build(o):
            o.EnterQty("E")
            o.ClcMaterial("Permittivity (epsi)", "mult")   # material multiply over gRPC
            o.EnterQty("E")
            o.CalcOp("Conj")
            o.CalcOp("Dot")
            o.CalcOp("Real")
            o.EnterVol(obj)
            o.CalcOp("Integrate")
        return self._evaluate(build)

    def line_voltage(self, line: str) -> float:
        """Signed peak line voltage ``sign(Re)·√(Re²+Im²)`` along ``line``, over gRPC."""
        def part(p):
            def build(o):
                o.EnterQty("E")
                o.CalcOp(p)                 # "Real" or "Imag"
                o.EnterLine(line)
                o.CalcOp("Tangent")
                o.CalcOp("Dot")
                o.EnterLine(line)
                o.CalcOp("Integrate")
            return build
        re = self._evaluate(part("Real"))
        im = self._evaluate(part("Imag"))
        mag = math.sqrt(re * re + im * im)
        return mag if re >= 0 else -mag


# ---------------------------------------------------------------------------
#  Session targeting (multi-session safe)
# ---------------------------------------------------------------------------
def _owning_session_from_lock(project_path: str, grpc: Dict[int, int]):
    """``(pid, port)`` of the running gRPC session that holds this project's lock.

    AEDT writes the owning ``DesktopProcessID`` into ``<project>.aedt.lock``;
    attaching there means we use the session the user already has the project open
    in, instead of an arbitrary (possibly wrong-version) instance.
    """
    lock = str(project_path) + ".lock"
    if not os.path.isfile(lock):
        return None
    pid = None
    try:
        for ln in open(lock, errors="ignore"):
            if ln.strip().startswith("DesktopProcessID="):
                pid = int(ln.strip().split("=", 1)[1])
                break
    except (OSError, ValueError):
        return None
    return (pid, grpc[pid]) if pid in grpc else None


# ---------------------------------------------------------------------------
#  Main analysis object
# ---------------------------------------------------------------------------
class PyaedtDistributedAnalysis:
    """EPR field extraction through PyAEDT (pure gRPC) — the PyAEDT analogue of
    :class:`pyEPR.DistributedAnalysis`.

    After :meth:`do_EPR_analysis` the results live on the object as plain NumPy
    arrays: ``freqs_GHz`` (M,), ``Ljs`` (J,), ``PJ`` (M by J ``p_mj``) and
    ``SJ`` (M by J signs) — the inputs pyEPR's diagonalizer needs.

    Parameters
    ----------
    project_info : pyEPR.ProjectInfo
        Construct it with ``do_connect=False`` — PyAEDT does the connecting.
        Junctions are read from ``project_info.junctions`` using the keys
        ``Lj_variable`` (or ``Lj_henries``), ``line`` and, optionally,
        ``rect`` / ``solid``.
    aedt_version : str
        AEDT version string, e.g. ``"2026.1"``.
    non_graphical : bool
        Forwarded to PyAEDT only when no running session is found.
    new_desktop : bool
        Forwarded to PyAEDT only when no running session is found; when a
        session is already running this object attaches to it instead.
    """

    def __init__(
        self,
        project_info,
        aedt_version: str = "2026.1",
        non_graphical: bool = False,
        new_desktop: bool = False,
    ):
        self.pinfo = project_info
        self.aedt_version = aedt_version
        self.non_graphical = non_graphical
        self.new_desktop = new_desktop

        self._hfss = None
        self.junction_names: List[str] = list(getattr(project_info, "junctions", {}) or {})

        # Results (populated by do_EPR_analysis)
        self.freqs_GHz: Optional[np.ndarray] = None
        self.Ljs: Optional[np.ndarray] = None
        self.PJ: Optional[np.ndarray] = None
        self.SJ: Optional[np.ndarray] = None
        self.results: Dict = {}

    # ---- connection --------------------------------------------------------
    @property
    def _project_path(self) -> str:
        """Full ``.aedt`` path from the ProjectInfo (``project_path`` + name)."""
        p = self.pinfo.project_path
        name = self.pinfo.project_name
        if p and name and not str(p).lower().endswith(".aedt"):
            return os.path.join(str(p), f"{name}.aedt")
        return str(p)

    def connect(self):
        """Attach to (or open) the project through PyAEDT — pure gRPC."""
        Hfss, active_sessions = _require_pyaedt()
        project_path = self._project_path

        sessions = {}
        try:
            sessions = active_sessions(
                version=self.aedt_version, student_version=False,
                non_graphical=False) or {}
        except Exception:
            sessions = {}
        grpc = {pid: port for pid, port in sessions.items() if port and port > 0}

        if grpc:
            owner = _owning_session_from_lock(project_path, grpc)
            pid, port = owner if owner else (next(iter(grpc)), grpc[next(iter(grpc))])
            self._hfss = Hfss(
                project=project_path, design=self.pinfo.design_name,
                version=self.aedt_version, port=port, aedt_process_id=pid,
                new_desktop=False, non_graphical=False)
        else:
            self._hfss = Hfss(
                project=project_path, design=self.pinfo.design_name,
                version=self.aedt_version, non_graphical=self.non_graphical,
                new_desktop=self.new_desktop)
        return self

    def disconnect(self):
        """Release the PyAEDT handle without closing a project we didn't open."""
        if self._hfss is not None:
            try:
                self._hfss.release_desktop(close_projects=False,
                                           close_desktop=self.new_desktop)
            except Exception:
                pass
            self._hfss = None

    def __enter__(self):
        return self.connect()

    def __exit__(self, *exc):
        self.disconnect()

    # ---- extraction --------------------------------------------------------
    def _resolve_Lj(self, jdict) -> float:
        """Junction inductance in Henries from a junction dict.

        Uses ``Lj_henries`` if given; otherwise reads the HFSS design variable
        named by ``Lj_variable`` (e.g. ``"12.9201nH"``) over gRPC and parses it.
        """
        if jdict.get("Lj_henries") is not None:
            return float(jdict["Lj_henries"])
        var = jdict.get("Lj_variable")
        if not var:
            raise ValueError("junction needs 'Lj_variable' or 'Lj_henries'")
        return _parse_henries(self._hfss.odesign.GetVariableValue(var))

    def _eigenfrequencies(self, setup_solution: str, variation: str = "") -> np.ndarray:
        """Eigenmode frequencies (GHz) via ExportEigenmodes on the gRPC Solutions module."""
        fd, path = tempfile.mkstemp(suffix=".txt", prefix="pyaedt_eig_")
        os.close(fd)
        try:
            self._hfss.odesign.GetModule("Solutions").ExportEigenmodes(
                setup_solution, variation, path)
            data = np.genfromtxt(path, dtype="str")
            if data.ndim == 1:
                data = np.array([data])
            return np.array([float(row[1]) for row in data])   # column 1 = GHz
        finally:
            try:
                os.remove(path)
            except OSError:
                pass

    def do_EPR_analysis(self, variation: str = "", modes: Optional[List[int]] = None):
        """Extract ``p_mj``, signs, frequencies and ``L_j`` for every mode — pure gRPC.

        Populates ``self.freqs_GHz``, ``self.Ljs``, ``self.PJ`` (M×J), ``self.SJ``
        (M×J). Returns ``self``.
        """
        if self._hfss is None:
            self.connect()
        setup_solution = f"{self.pinfo.setup_name} : LastAdaptive"
        fc = _GrpcFieldCalc(self._hfss, setup_solution)

        freqs = self._eigenfrequencies(setup_solution, variation)
        if modes is None:
            modes = list(range(len(freqs)))
        freqs = freqs[modes]
        M = len(modes)

        jnames = self.junction_names
        J = len(jnames)
        jdicts = [self.pinfo.junctions[n] for n in jnames]
        Ljs = np.array([self._resolve_Lj(jd) for jd in jdicts]) if J else np.array([])

        PJ = np.zeros((M, J))
        SJ = np.ones((M, J), dtype=int)

        for mi, mode in enumerate(modes):
            fc.set_mode(mode)                          # EditSources normalization (gRPC)
            U_E_raw = fc.energy_electric("AllObjects")
            Ucaps, Vs = [], []
            for jd in jdicts:
                V = fc.line_voltage(jd["line"])
                Vs.append(V)
                Ucaps.append(0.5 * float(jd.get("Cj_farads", 0.0) or 0.0) * V * V)
            for ji, jd in enumerate(jdicts):
                p, s = compute_p_mj(
                    V_peak=Vs[ji], freq_hz=float(freqs[mi]) * 1e9,
                    Lj_henries=float(Ljs[ji]), U_E_raw=U_E_raw,
                    Cj_farads=float(jd.get("Cj_farads", 0.0) or 0.0),
                    sum_Ucap_other=sum(Ucaps) - Ucaps[ji],
                )
                PJ[mi, ji] = p
                SJ[mi, ji] = s

        self.freqs_GHz, self.Ljs, self.PJ, self.SJ = freqs, Ljs, PJ, SJ
        self.results = {
            "freqs_hfss_GHz": freqs, "Ljs": Ljs, "Pm": PJ, "Sm": SJ,
            "junctions": jnames, "modes": modes,
        }
        return self

    # ---- hand off to pyEPR's physics --------------------------------------
    def epr_zpf(self) -> np.ndarray:
        """Reduced zero-point flux ``ϕ_zpf`` (M×J) via pyEPR's ``epr_to_zpf``."""
        from .calcs.basic import CalcsBasic

        if self.PJ is None:
            raise RuntimeError("call do_EPR_analysis() first")
        H = 6.62607015e-34
        PHI0_RED = 3.29105976e-16
        EJ_GHz = np.array([(PHI0_RED ** 2 / (Lj * H)) / 1e9 for Lj in self.Ljs])
        return CalcsBasic.epr_to_zpf(
            self.PJ, self.SJ, np.diag(self.freqs_GHz), np.diag(EJ_GHz))

    def analyze(self, cos_trunc: int = 8, fock_trunc: int = 7):
        """Numerically diagonalize via pyEPR's ``epr_numerical_diagonalization``.

        Returns ``(f_ND, chi_ND)`` — dressed frequencies (GHz) and the cross-Kerr /
        anharmonicity matrix (MHz, diagonal = anharmonicity).
        """
        from .calcs.back_box_numeric import epr_numerical_diagonalization

        phi = self.epr_zpf()
        return epr_numerical_diagonalization(
            self.freqs_GHz, np.asarray(self.Ljs), phi,
            cos_trunc=cos_trunc, fock_trunc=fock_trunc)
