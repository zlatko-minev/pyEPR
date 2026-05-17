Welcome to pyEPR :beers:! &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(see [arXiv:2010.00620](https://arxiv.org/abs/2010.00620))
===================
### Automated Python module for the design and quantization of Josephson quantum circuits

[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.png?v=103)](https://github.com/zlatko-minev/pyEPR)
[![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/zlatko-minev/pyEPR)
[![Star pyEPR](https://img.shields.io/badge/⭐-Star%20pyEPR-blue?style=flat)](https://github.com/zlatko-minev/pyEPR/stargazers)
[![Fork pyEPR](https://img.shields.io/badge/🍴-Fork%20pyEPR-blue?style=flat)](https://github.com/zlatko-minev/pyEPR/fork)
<br>
[![PyPI version](https://badge.fury.io/py/pyEPR-quantum.svg)](https://badge.fury.io/py/pyEPR-quantum)
[![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)
[![CI](https://github.com/zlatko-minev/pyEPR/actions/workflows/ci.yaml/badge.svg)](https://github.com/zlatko-minev/pyEPR/actions/workflows/ci.yaml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zlatko-minev/pyEPR/master?filepath=_tutorial_notebooks%2FTutorial%206.%20EPR%20without%20HFSS%20%E2%80%94%20purely%20numerical%20workflow.ipynb)


## What is pyEPR?

**pyEPR** implements the *energy-participation ratio* (EPR) method for the design and quantization of superconducting quantum circuits. The EPR framework, introduced in [Minev *et al.*, npj Quantum Information (2021)](https://www.nature.com/articles/s41534-021-00461-8) ([arXiv:2010.00620](https://arxiv.org/abs/2010.00620)), bridges classical electromagnetic simulation and quantum circuit theory:

1. **Classical EM** — simulate the linearized circuit in Ansys HFSS to obtain eigenmode frequencies and field distributions.
2. **EPR extraction** — compute the energy participation ratios $p_{mj}$ of each Josephson junction $j$ in each mode $m$, and the zero-point phase fluctuations $\varphi_{\rm zpf}^{(mj)} = \sqrt{p_{mj}\,\hbar\omega_m / (2E_J)}$.
3. **Quantum Hamiltonian** — diagonalize $H = \sum_m \omega_m a_m^\dagger a_m - \sum_j E_J [\cos(\hat\varphi_j) - 1 + \hat\varphi_j^2/2]$ numerically to get dressed frequencies and the dispersive-shift matrix $\chi$.

The method works for any number of modes and junctions, handles strongly anharmonic circuits (fluxonium, $\varphi_{\rm zpf} \gtrsim 1$), and requires no manual circuit diagram — only the 3D geometry in HFSS.

### Documentation

[Read the docs here.](https://pyepr-docs.readthedocs.io)
<br>

## How to cite

If you use pyEPR in your research, please cite:

* **EPR method paper** (primary reference for the method):
  Z. K. Minev, Z. Leghtas, S. O. Mundhada, L. Christakis, I. M. Pop, M. H. Devoret,
  *Energy-participation quantization of Josephson circuits*,
  [npj Quantum Information **7**, 131 (2021)](https://doi.org/10.1038/s41534-021-00461-8) · [arXiv:2010.00620](https://arxiv.org/abs/2010.00620)

* **Software** (Zenodo DOI for the specific version):
  [![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)

Use this [BibTeX file](https://github.com/zlatko-minev/pyEPR/blob/master/pyEPR.bib) for both.

Related references:
* Z. K. Minev, Ph.D. Dissertation, Yale University (2018), Chapter 4. ([arXiv:1902.10355](https://arxiv.org/abs/1902.10355))
* A. Petrescu, C. T. Hann, Z. K. Minev *et al.*, EPR for very anharmonic circuits. ([arXiv:2411.15039](https://arxiv.org/abs/2411.15039))

## Who uses pyEPR?
* Yale University, Michel Devoret lab [QLab](https://qulab.eng.yale.edu/), CT, USA
* Yale University, Rob Schoelkopf lab [RSL](https://rsl.yale.edu/), CT, USA
* [IBM Quantum](https://www.ibm.com/quantum-computing/) and IBM's Qiskit Metal
* [QUANTIC](https://team.inria.fr/quantic/people.html#) (QUANTUM INFORMATION CIRCUITS), PARIS — INRIA, ENS, MINES PARISTECH, UPMC, CNRS. Groups of Zaki Leghtas and team. France
* [Quantum Circuit Group](http://www.physinfo.fr/) Benjamin Huard, Ecole Normale Supérieure de Lyon, France
* Emanuel Flurin, CEA Saclay, France
* Ioan Pop group, KIT Physikalisches Institut, Germany
* UC Berkeley, [Quantum Nanoelectronics Laboratory](https://physics.berkeley.edu/quantum-nanoelectronics-laboratory), Irfan Siddiqi, CA, USA
* [Quantum Circuits, Inc.](https://quantumcircuits.com/), CT, USA
* [Seeqc](https://seeqc.com/) (spin-out of Hypres) Digital Quantum Computing, USA
* Serge [Rosenblum Lab](https://www.weizmann.ac.il/condmat/rosenblum/) quantum circuits group, Weizmann Institute, Israel
* University of Oxford — [Leek Lab](https://leeklab.org/), UK
* Britton [Plourde Lab](https://bplourde.expressions.syr.edu/), Syracuse University, USA
* Javad [Shabani Lab](https://wp.nyu.edu/shabanilab/) Quantum Materials & Devices, NYU, NY, USA
* UChicago Dave Schuster Lab, USA
* SQC lab — Shay Hacohen-Gourgy, Israel
* Lawrence Berkeley National Lab
* Colorado School of Mines, USA
* IPQC, SJTU, Shanghai, China
* Bhabha Atomic Research Centre, India
* Quantum Computing UK
* Alice&Bob, France
* Centre for Quantum Technologies / Qcrew
* Quantum Device Lab ETHZ — Andreas Wallraff
* Bleximo
* ... and many more! (Please e-mail `zlatko.minev@aya.yale.edu` with updates.)


<br>

# Contents:
* [Start here: Using `pyEPR`](#start-here-using-pyepr)
* [Tutorial Notebooks](#tutorial-notebooks)
* [Video Tutorials](#pyepr-video-tutorials)
* [Setup and Installation](#installation-and-setup-of-pyepr)
* [HFSS Project Setup for `pyEPR`](#hfss-project-setup-for-pyepr)
* [Troubleshooting `pyEPR`](#troubleshooting-pyepr)
* [Authors and Contributors](#authors-and-contributors)


![Intro image](imgs/read_me_0.png 'Intro image')

# Start here: Using `pyEPR`

1. **Install** — see [Installation and setup](#installation-and-setup-of-pyepr) below. The fastest path: `pip install pyEPR-quantum`.
2. **Tutorials** — work through the [Jupyter notebook tutorials](#tutorial-notebooks) below. No HFSS licence? Start with Tutorial 6.
3. **Read the docs** — [pyepr-docs.readthedocs.io](https://pyepr-docs.readthedocs.io) for the API reference and detailed guides.
4. **Cite `pyEPR`** — [arXiv:2010.00620](https://arxiv.org/abs/2010.00620) [![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)

## Quickstart — no Ansys required

No HFSS licence? You can run the full quantum Hamiltonian diagonalization from scratch in three steps. The only requirements are `numpy`, `qutip`, and `matplotlib` (all installed with `pip install pyEPR-quantum`).

```python
import numpy as np
from pyEPR.calcs.back_box_numeric import epr_numerical_diagonalization

# Standard transmon: E_J/h = 20 GHz, E_C/h = 300 MHz
# Plasma frequency f_p = sqrt(8 E_J E_C)/h ≈ 6.93 GHz
freqs   = np.array([6.928])     # GHz  (linearised plasma frequency)
Ljs     = np.array([8.2e-9])    # H    (Josephson inductance L_J = (Φ₀/2π)² / E_J)
phi_zpf = np.array([[0.416]])   # dimensionless  (n_modes × n_junctions)

# Diagonalize — returns dressed frequencies in Hz and χ matrix in MHz
f_dressed, chi_matrix = epr_numerical_diagonalization(
    freqs, Ljs, phi_zpf, cos_trunc=8, fock_trunc=15
)
print(f"Qubit frequency : {f_dressed[0].real / 1e9:.3f} GHz")
print(f"Anharmonicity   : {chi_matrix[0,0].real:.0f} MHz")
# → Qubit frequency : 6.614 GHz
# → Anharmonicity   : 337 MHz
```

For a multi-mode transmon + resonator, fluxonium, or a custom junction potential, see **[Tutorial 6](https://mybinder.org/v2/gh/zlatko-minev/pyEPR/master?filepath=_tutorial_notebooks%2FTutorial%206.%20EPR%20without%20HFSS%20%E2%80%94%20purely%20numerical%20workflow.ipynb)** (runs in-browser via Binder, no install needed).


#### Start-up example

The following code illustrates how to perform a complete EPR analysis of a simple two-qubit, one-cavity device in just a few lines. In the HFSS file, first specify the non-linear junction rectangles and variables (see [HFSS Project Setup](#hfss-project-setup-for-pyepr) below). All operations in the eigenmode analysis and Hamiltonian computation are fully automated.

```python
import pyEPR as epr

# 1. Connect to Ansys HFSS and load your design
pinfo = epr.ProjectInfo(project_path = r'C:\sim_folder',
                        project_name = r'cavity_with_two_qubits',
                        design_name  = r'Alice_Bob')

# 2a. Specify Josephson junctions (non-linear elements)
pinfo.junctions['jAlice'] = {'Lj_variable':'Lj_alice', 'rect':'rect_alice',
                              'line': 'line_alice', 'Cj_variable':'Cj_alice'}
pinfo.junctions['jBob']   = {'Lj_variable':'Lj_bob',   'rect':'rect_bob',
                              'line': 'line_bob',   'Cj_variable':'Cj_bob'}
pinfo.validate_junction_info()

# 2b. (Optional) Dissipative elements for loss calculation
pinfo.dissipative['dielectrics_bulk']    = ['si_substrate']
pinfo.dissipative['dielectric_surfaces'] = ['interface1', 'interface2']

# 3. Classical EM analysis — extract EPR participation ratios φ_zpf
eprd = epr.DistributedAnalysis(pinfo)
eprd.hfss_report_full_convergence()   # check HFSS convergence
eprd.do_EPR_analysis()                # compute p_mj and φ_zpf for all modes/junctions

# 4. Quantum Hamiltonian analysis — dress frequencies and compute χ matrix
epra = epr.QuantumAnalysis(eprd.data_filename)
epra.analyze_all_variations(cos_trunc=8, fock_trunc=15)

# 5. Report results
swp_variable = 'Lj_alice'
epra.plot_hamiltonian_results(swp_variable=swp_variable)
epra.report_results(swp_variable=swp_variable, numeric=True)
```

# Tutorial Notebooks

The tutorials are Jupyter notebooks in the [`_tutorial_notebooks/`](https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks) folder. They build on each other — work through them in order for the best experience.

| # | Title | Requires HFSS? | Topics |
|---|---|---|---|
| 1 | [Startup example](https://github.com/zlatko-minev/pyEPR/blob/master/_tutorial_notebooks/Tutorial%201.%20%20Startup%20example.ipynb) | Yes | Full end-to-end workflow: connect to HFSS, define junctions, run EPR, get χ matrix |
| 2 | [Dielectric loss EPR](https://github.com/zlatko-minev/pyEPR/blob/master/_tutorial_notebooks/Tutorial%202.%20%20Field%20calculations%20-%20dielectric%20energy%20participation%20ratios%20(EPRs).ipynb) | Yes | Dielectric energy participation ratios, loss rates, HFSS fields calculator |
| 3 | [Circuit QED parameters](https://github.com/zlatko-minev/pyEPR/blob/master/_tutorial_notebooks/Tutorial%203.%20%20toolbox_circuits.ipynb) | No | Josephson junction physics, E_J, E_C, L_J, I_c conversions, transmon model |
| 4 | [Parametric sweeps](https://github.com/zlatko-minev/pyEPR/blob/master/_tutorial_notebooks/Tutorial%204.%20Parametric%20sweep%20options.ipynb) | Yes | All HFSS Optimetrics sweep types: linear, log, file-based |
| 5 | [Generic junction potential & fluxonium](https://github.com/zlatko-minev/pyEPR/blob/master/_tutorial_notebooks/Tutorial%205.%20Generic%20junction%20potential%20and%20fluxonium%20EPR.ipynb) | No | Exact cosine for fluxonium, custom potentials, asymmetric SQUIDs |
| 6 | [EPR without HFSS](https://github.com/zlatko-minev/pyEPR/blob/master/_tutorial_notebooks/Tutorial%206.%20EPR%20without%20HFSS%20%E2%80%94%20purely%20numerical%20workflow.ipynb) | No | Purely numerical workflow: supply freqs, Ljs, φ_zpf directly |

> **No HFSS?** Start with Tutorials 3, 5, and 6, which are self-contained and require only `pip install pyEPR-quantum`.

# `pyEPR` Video Tutorials <img src="https://developers.google.com/site-assets/logo-youtube.svg" height=30>
<div style="overflow:auto;">
<table style="">
  <tr>
    <th>
	<a href="https://www.youtube.com/watch?v=fSRYvD-ITnQ&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=1">
	Tutorial 1 - Overview
	<br>
	<img src="https://img.youtube.com/vi/fSRYvD-ITnQ/0.jpg" alt="pyEPR Tutorial 1 - Overview" width=250>
	</a>
    </th>
    <th>
	<a href="https://www.youtube.com/watch?v=ZTi1pb6wSbE&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=2">
	Tutorial 2 - Setup of Conda & Git
	<br>
	<img src="https://img.youtube.com/vi/ZTi1pb6wSbE/0.jpg" alt="pyEPR Tutorial 2 - Setup of Conda & Git" width=250>
	</a>
    </th>
    <th>
	<a href="https://www.youtube.com/watch?v=L79nlXY2w4s&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=3">
	Tutorial 3 - Setup of Packages & Config
	<br>
	<img src="https://img.youtube.com/vi/L79nlXY2w4s/0.jpg" alt="pyEPR Tutorial 3 - Setup of Packages & Config" width=250>
	</a>
    </th>
  </tr>
</table>
</div>

# Installation and setup of `pyEPR`
-------------

**Requires Python 3.9–3.12.**  All dependencies are installed automatically.

### Quick install

**uv** (recommended — fast, modern):
```sh
uv pip install pyEPR-quantum
```

**pip**:
```sh
pip install pyEPR-quantum
```

**conda** (conda-forge channel):
```sh
conda install -c conda-forge pyepr-quantum
```

> The conda-forge package name is `pyepr-quantum` (lower-case); the PyPI name is `pyEPR-quantum`.
> Either way you import it as `import pyEPR as epr`.

### Development / editable install

```sh
git clone https://github.com/zlatko-minev/pyEPR.git
cd pyEPR
pip install -e ".[test]"
pytest          # runs all tests that don't need a live HFSS session
```

### Platform support

pyEPR has two layers with different platform requirements:

| Feature | Windows | macOS | Linux |
|---|---|---|---|
| EPR / quantum analysis (`DistributedAnalysis`, `QuantumAnalysis`) | ✅ | ✅ | ✅ |
| Ansys HFSS COM interface (`HfssDesign`, `ansys.py`) | ✅ | ⚠️ limited | ⚠️ limited |

**EPR and quantum analysis** are pure-Python and work on all platforms — no Ansys installation required for post-processing or for the fully numerical workflow (Tutorials 3, 5, 6).

**The HFSS COM interface** (`ansys.py`) was written for Windows, where HFSS exposes automation via `pythoncom`/`win32com`. On macOS/Linux you can still reach a remote Windows HFSS instance over a network COM bridge.

> **Note on PyAEDT:** [PyAEDT](https://github.com/ansys/pyaedt) is Ansys's official cross-platform Python scripting library for AEDT. pyEPR predates it and focuses on the quantum EPR quantization workflow that PyAEDT does not cover. They are complementary: use PyAEDT for geometry/mesh/solve scripting, pyEPR for EPR-based Hamiltonian extraction.

### Optional configuration

Edit `pyEPR/_config_user.py` to set your data-save directory and logging level.
To keep local changes out of git:
```sh
git update-index --skip-worktree pyEPR/_config_user.py
```

#### Legacy users
pyEPR was significantly refactored in v0.8 (2020). Key classes were renamed — see the tutorials and docs. If you cannot upgrade yet, use the stable v0.7 tag.


# HFSS Project Setup for `pyEPR`
-------------
#### Eigenmode Design — How to set up junctions

The EPR method requires each Josephson junction to be modelled as a **lumped RLC boundary** on a rectangle in HFSS, together with a **polyline** spanning the junction that defines the current direction. Follow these steps:

 1. **Define circuit geometry and electromagnetic boundary conditions.**
    1. *Junction rectangle and BC:* Create a rectangle for each Josephson junction and give it a descriptive name, e.g., `rect_alice` for a qubit named Alice. A 50 × 100 μm rectangle is typical, though smaller sizes work fine. Assign a `Lumped RLC` boundary condition on this rectangle surface, with an inductance value given by a local variable (e.g., `Lj_alice`). The variable name is passed to pyEPR.
    2. *Current-direction polyline:* Over each junction rectangle, draw a model `polyline` spanning the full length of the rectangle to define the current-flow direction. Best practice: define it in an object coordinate system on the junction rectangle so they move together when the geometry is modified. The polyline name is passed to pyEPR.
    3. *Note:* The linearised Josephson inductance used in HFSS is $L_J = (\Phi_0/2\pi)^2 / E_J$. This is the small-signal inductance that appears in the classical eigenmode problem; the full nonlinear cosine potential is restored in the quantum Hamiltonian step.
 2. **Meshing.**
    1. Lightly mesh the thin-film metal boundary condition.
    2. Lightly mesh the junction rectangles.
 3. **Simulation setup.**
    1. We recommend `mixed order` solutions for best accuracy.
    2. Make sure "Save Fields" is enabled for each variation if running a parametric sweep (see Tutorial 4).

<p align="center">
  <img width="50%" height="50%" src="imgs/xmon-example.gif">
</p>


# Troubleshooting pyEPR
---------------------
###### First run: pint error: system='mks' unknown.
Please update to pint version newer than 0.7.2:
```
pip install pint --upgrade
```

###### QuTiP installation
QuTiP is a required dependency and is installed automatically with `pip install pyEPR-quantum`.
If you need to install it separately:
```sh
pip install qutip          # PyPI (qutip >= 5.0 required)
conda install -c conda-forge qutip   # conda-forge
```

###### `AttributeError: module 'numpy' has no attribute '__config__'` on importing qutip
Update your numpy installation:
```
conda update qutip
conda update numpy
```

###### COM Error on opening HFSS
Check the project and design file names carefully. Make sure the file path does not contain apostrophes or other special characters (e.g., `C:\Minev's PC\my:Project`). Check that HFSS has not popped up an error dialogue such as "File locked." Manually open HFSS and the file.

###### COM error on calculation of expression
Either HFSS has popped an error dialog, frozen up, or a name is misspelled. Verify all object and variable names match exactly.

###### HFSS refuses to close
If your script terminates improperly, this can happen. pyEPR tries to catch termination events and handle them. Call `hfss.release()` when finished. Use the Task Manager (Activity Monitor on macOS) to kill HFSS if necessary.

###### Parametric sweep error — missing field solutions
When running a parametric sweep in HFSS, make sure you save the fields for each variation before running pyEPR. Right-click on your ParametricSetup → Properties → Options → "Save Fields and Mesh". See Tutorial 4 for details.

###### `ValueError: cannot set WRITEABLE flag to True of this array`
This occurs with numpy 1.16 and certain HDF files. Upgrade numpy or downgrade to 1.15.4.

# Authors and Contributors
* *Authors:* [Zlatko Minev](https://www.zlatko-minev.com/) & [Zaki Leghtas](http://cas.ensmp.fr/~leghtas/), with contributions from many friends and colleagues.
* 2015 – present.
* Contributors: [Phil Reinhold](https://github.com/PhilReinhold), Lysander Christakis, [Devin Cody](https://github.com/devincody), Zachary Parrott, and many others.
* Original versions of pyHFSS.py and pyNumericalDiagonalization.py contributed by [Phil Reinhold](https://github.com/PhilReinhold) — excellent original [repo](https://github.com/PhilReinhold/pyHFSS).
* Terms of use: Use freely and kindly cite the paper ([arXiv:2010.00620](https://arxiv.org/abs/2010.00620)) and/or this package.
* How can I contribute? Contact [Z. Minev](https://www.zlatko-minev.com/) or [Z. Leghtas](http://cas.ensmp.fr/~leghtas/).  [![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)

## How do I cite `pyEPR`?
[![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)

Use this [BibTeX file](https://github.com/zlatko-minev/pyEPR/blob/master/pyEPR.bib) for `pyEPR` and cite the energy-participation-ratio paper [arXiv:2010.00620](https://arxiv.org/abs/2010.00620) for the method.


[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/zlatko-minev/pyEPR/graphs/commit-activity)

[![Twitter](https://github.frapsoft.com/social/twitter.png)](https://twitter.com/zlatko_minev)
[![Github](https://github.frapsoft.com/social/github.png)](https://github.com/zlatko-minev/)
