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


### Documentation

[Read the docs here.](https://pyepr-docs.readthedocs.io)
<br>

## Scientific work:
* Minev, Z. K., Leghtas, Z., Mudhada, S. O., Reinhold, P., Diringer, A., & Devoret, M. H. (2018). [pyEPR: The energy-participation-ratio (EPR) open-source framework for quantum device design.](https://github.com/zlatko-minev/pyEPR/blob/master/pyEPR.bib)  [![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)
* Minev, Z. K., Leghtas, Z., Mundhada, S. O., Christakis, L., Pop, I. M., & Devoret, M. H. (2020). Energy-participation quantization of Josephson circuits. ArXiv. Retrieved from http://arxiv.org/abs/2010.00620 (2020)
* Z.K. Minev, Ph.D. Dissertation, Yale University (2018), Chapter 4. ([arXiv:1902.10355](https://arxiv.org/abs/1902.10355))  (2018) 

## pyEPR Working group meeting -- Planning for the future of pyEPR

* Please sign-up here: https://github.com/zlatko-minev/pyEPR/issues/45 or [directly here](https://docs.google.com/forms/d/e/1FAIpQLScd3WyfzDS47D0WB9skkSPQAXCnKLf7JMxsZ7BnMwK0LjE3Sw/viewform?usp=sf_link) :bangbang: :beers:
- See [pyEPR wiki](https://github.com/zlatko-minev/pyEPR/wiki) for notes from first meeting.
- We will schedule a follow-up meeting in 1-2 mo.

<br>

## Who uses pyEPR?
* Yale University, Michel Devoret lab [QLab](https://qulab.eng.yale.edu/), CT, USA
* Yale University, Rob Schoelkopf lab [RSL](https://rsl.yale.edu/), CT, USA
* [IBM Quantum](https://www.ibm.com/quantum-computing/) and IBM's Qiskit Metal
* [QUANTIC](https://team.inria.fr/quantic/people.html#) (QUANTUM INFORMATION CIRCUITS), PARISINRIA, ENS, MINES PARISTECH, UPMC, CNRS. Groups of Zaki Leghtas and team. France
* [Quantum Circuit Group](http://www.physinfo.fr/) Benjamin Huard, Ecole Normale Supérieure de Lyon, France
* Emanuel Flurin, CEA Saclay, France
* Ioan Pop group, KIT Physikalisches Institut, Germany 
* UC Berkeley, [Quantum Nanoelectronics Laboratory](https://physics.berkeley.edu/quantum-nanoelectronics-laboratory), Irfan Siddiqi, CA, USA
* [Quantum Circuits, Inc.](https://quantumcircuits.com/), CT, USA
* [Seeqc](https://seeqc.com/) (spin-out of Hypres) Digital Quantum Computing, USA
* Serge [Rosenblum Lab] quantum circuits group (https://www.weizmann.ac.il/condmat/rosenblum/) in the Weizmann Instatue, Israel
* University of Oxford - LeekLab - Peter [Leek Lab](https://leeklab.org/), UK
* Britton [Plourde Lab](https://bplourde.expressions.syr.edu/), Syracuse University, USA
* Javad [Shabani Lab](https://wp.nyu.edu/shabanilab/) Quantum Materials & Devices, NYU, NY, USA
* UChicago Dave Schuster Lab, USA
* SQC lab - Shay Hacohen Gourgy, Israel
* Lawrence Berkeley National Lab
* Colorado School of Mines, USA
* Syracuse University, USA
* IPQC, SJTU, Shanghai, China
* Bhabha Atomic Research Centre, India
* Quantum Computing UK
* Alice&Bob, France
* Centre for Quantum Technologies / Qcrew
* Quantum Device Lab ETHZ; Andreas Wallraff
* Bleximo
* ... and many more! (Please e-mail `zlatko.minev@aya.yale.edu` with updates.)


<br>

# Contents:
* [Start here: Using `pyEPR`](#start-here-using-pyepr)
* [Video Tutorials](#pyepr-video-tutorials)
* [Setup and Installation](#installation-and-setup-of-pyepr)
* [HFSS Project Setup for `pyEPR`](#hfss-project-setup-for-pyepr)
* [Troubleshooting `pyEPR`](#troubleshooting-pyepr)
* [Authors and Contributors](#authors-and-contributors)


![Intro image](imgs/read_me_0.png 'Intro image')

# Start here: Using `pyEPR`

1. **Install** — see [Installation and setup](#installation-and-setup-of-pyepr) below. The fastest path: `pip install pyEPR-quantum`.
2. **Tutorials** — work through the [Jupyter notebook tutorials](https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks) to learn the full workflow.
3. **Read the docs** — [pyepr-docs.readthedocs.io](https://pyepr-docs.readthedocs.io) for the API reference and detailed guides.
4. **Cite `pyEPR`** — [arXiv:2010.00620](https://arxiv.org/abs/2010.00620) / [arXiv:1902.10355](https://arxiv.org/abs/1902.10355) [![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)



#### Start-up example

[Jupyter notebook tutorials](https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks)

The following code illustrates how to perform a complete analysis of a simple two-qubit, one-cavity device in just a few lines of code with `pyEPR`.  In the HFSS file, before running the script, first specify the non-linear junction rectangles and variables (see Sec. pyEPR Project Setup in HFSS). All operations in the eigen analysis and Hamiltonian computation are fully automated. The results are saved, printed, and succinctly plotted.


```python
# Load pyEPR. See the tutorial notebooks!
import pyEPR as epr

# 1. Connect to your Ansys, and load your design
pinfo = epr.ProjectInfo(project_path = r'C:\sim_folder',
                        project_name = r'cavity_with_two_qubits',
                        design_name  = r'Alice_Bob')


# 2a. Non-linear (Josephson) junctions
pinfo.junctions['jAlice'] = {'Lj_variable':'Lj_alice', 'rect':'rect_alice', 'line': 'line_alice', 'Cj_variable':'Cj_alice'}
pinfo.junctions['jBob']   = {'Lj_variable':'Lj_bob',   'rect':'rect_bob',   'line': 'line_bob', 'Cj_variable':'Cj_bob'}
pinfo.validate_junction_info() # Check that valid names of variables and objects have been supplied.

# 2b. Dissipative elements: specify
pinfo.dissipative['dielectrics_bulk']    = ['si_substrate', 'dielectric_object2'] # supply names of hfss objects
pinfo.dissipative['dielectric_surfaces'] = ['interface1', 'interface2']
# Alternatively, these could be specified in ProjectInfo with
# pinfo = epr.ProjectInfo(..., dielectrics_bulk = ['si_substrate', 'dielectric_object2'])

# 3.  Perform microwave analysis on eigenmode solutions
eprd = epr.DistributedAnalysis(pinfo)
if 1: # automatic reports
  eprd.quick_plot_frequencies(swp_var) # plot the solved frequencies before the analysis
  eprd.hfss_report_full_convergence() # report convergence
eprd.do_EPR_analysis()

# 4a.  Perform Hamiltonian spectrum post-analysis, building on mw solutions using EPR
epra = epr.QuantumAnalysis(eprd.data_filename)
epra.analyze_all_variations(cos_trunc = 8, fock_trunc = 7)

# 4b. Report solved results
swp_variable = 'Lj_alice' # suppose we swept an optimetric analysis vs. inductance Lj_alice
epra.plot_hamiltonian_results(swp_variable=swp_variable)
epra.report_results(swp_variable=swp_variable, numeric=True)
epra.quick_plot_mode(0,0,1,numeric=True, swp_variable=swp_variable)
```

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
<!--
 [![pyEPR Tutorial 1 - Overview](https://img.youtube.com/vi/fSRYvD-ITnQ/0.jpg)](https://www.youtube.com/watch?v=fSRYvD-ITnQ&list=PLnak_fVcHp17tydgFosNtetDNjKbEaXtv&index=1) -->

[Jupyter notebook tutorials](https://github.com/zlatko-minev/pyEPR/tree/master/_tutorial_notebooks)

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
pip install -e “.[test]”
pytest          # runs all tests that don’t need a live HFSS session
```

### Platform support

pyEPR has two layers with different platform requirements:

| Feature | Windows | macOS | Linux |
|---|---|---|---|
| EPR / quantum analysis (`DistributedAnalysis`, `QuantumAnalysis`) | ✅ | ✅ | ✅ |
| Ansys HFSS COM interface (`HfssDesign`, `ansys.py`) | ✅ | ⚠️ limited | ⚠️ limited |

**EPR and quantum analysis** are pure-Python and work on all platforms — no Ansys installation required for post-processing.

**The HFSS COM interface** (`ansys.py`) was written for Windows, where HFSS exposes automation via `pythoncom`/`win32com`.  On macOS/Linux you can still reach a remote Windows HFSS instance over a network COM bridge.

> **Note on PyAEDT:** [PyAEDT](https://github.com/ansys/pyaedt) is Ansys’s official cross-platform Python scripting library for AEDT.  pyEPR predates it and focuses on the quantum EPR quantization workflow that PyAEDT does not cover.  They are complementary: use PyAEDT for geometry/mesh/solve scripting, pyEPR for EPR-based Hamiltonian extraction.

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
#### Eigenmode Design --- How to set up junctions
You may find an advised work flow and some setup tips here.

 1. Define circuit geometry & electromagnetic boundary condition (BC).
   1. Junction rectangles and BC: Create a rectangle for each Josephson junction and give it a good name; e.g., `jAlice` for a qubit named Alice. We recommend 50 x 100 um rectangle for a simple simulation, although orders of magnitude smaller rectangles work as well. Note the length of this junction, you will supply it to pyEPR. Assign a `Lumped RLC` BC on this rectangle surface, with an inductance value given by a local variable, `Lj1` for instance. The name of this variable will also be supplied to the pyEPR.
   2. Over each junction rectangle draw a model `polyline` to define give a sense of the junction current-flow direction. This line should spans the length of the full junction rectangle. Define it using an object coordinate system on the junction rectangle (so that they move together when the geometry is altered). The name of this line will be supplied to the pyEPR module.
 2. Meshing.
   1. Lightly mesh the thin-film metal BC. Lightly mesh the junction rectangles.
 3. Simulation setup
   1. We recommend `mixed order` solutions.

<p align="center">
  <img width="50%" height="50%" src="imgs/xmon-example.gif">
</p>


# Troubleshooting pyEPR
---------------------
###### First run: pint error: system='mks' unknown.
Please update to pint version newer than 0.7.2. You may use
```
pip install pint --upgrade
```

###### No attribute `StringIO` during do_EPR_analysis()

`AttributeError: module 'pandas.compat' has no attribute 'StringIO'`

Caleb pointed this out, see [here](https://stackoverflow.com/questions/58372475/attributeerror-module-pandas-compat-has-no-attribute-iteritems) and here for [solution](https://github.com/zlatko-minev/pyEPR/issues/21). You need to change the pandas version. [pyEPR to be upgraded]

This was solved in [this commit](https://github.com/DanielCohenHillel/pyEPR/commit/fd2b5897d6f819681b8a605734cfed855c002df6). Try to update your pyEPR version to the current master.


###### When importing qutip an error occurs `AttributeError: module 'numpy' has no attribute '__config__'`
You probably have to update your numpy installation. For me, the following bash sequence worked:
```
conda update qutip
conda update numpy
```

###### QuTiP installation
QuTiP is a required dependency and is installed automatically with `pip install pyEPR-quantum`.
If you need to install it separately:
```sh
pip install qutip          # PyPI (qutip >= 5.0 required)
conda install -c conda-forge qutip   # conda-forge
```


###### COM Error on opening HFSS
Check the project and design file names carefully. Make sure that the file-path doesn't have apostrophes or other bad characters, such as in C:\\Minev's PC\\my:Project.  Check that HFSS hasn't popped up an error dialogue, such as "File locked." Manually open HFSS and the file.

###### COM error on calculation of expression
Either HFSS popped an error dialog, froze up, or you miss-typed the name of something.

###### HFSS refuses to close
If your script terminates improperly, this can happen. pyHFSS tries to catch termination events and handle them. Your safety should be guaranteed however, if you call `hfss.release()` when you have finished. Use the Task-manager (Activity Monitor on MAC) to kill HFSS if you want.

###### Parametric Sweep Error
When running a parametric sweep in HFSS, make sure you are actually saving the fields for each variation before running pyEPR. This can be done by right-clicking on your ParametricSetup -> properties -> options -> "Save Fields and Mesh".

###### Spyder pops up command window cmd with tput.exe executed
This problem is due to pandas 0.20.1, update to 0.20.3 or better solves this issue.
<br>

###### `ValueError: cannot set WRITEABLE flag to True of this array`
This error happens when trying to read in an hdf file with numpy version 1.16, see [git issue here](https://github.com/numpy/numpy/issues/12791). A solution is to downgrade numpy to 1.15.4 or upgrade to newer versions of hdf and numpy.

# Authors and Contributors
* _Authors:_ [Zlatko Minev](https://www.zlatko-minev.com/) & [Zaki Leghtas](http://cas.ensmp.fr/~leghtas/), with contributions from many friends and colleagues. ([arXiv:2010.00620](https://arxiv.org/abs/2010.00620))
* 2015 - present.
* Contributors: [Phil Rheinhold](https://github.com/PhilReinhold), Lysander Christakis, [Devin Cody](https://github.com/devincody), ...
Original versions of pyHFSS.py and pyNumericalDiagonalization.py contributed by [Phil Rheinhold](https://github.com/PhilReinhold), excellent original [repo](https://github.com/PhilReinhold/pyHFSS).
* Terms of use: Use freely and kindly cite the paper (arXiv link to be posted here) and/or this package.
* How can I contribute? Contact [Z. Minev](https://www.zlatko-minev.com/) or [Z. Leghtas](http://cas.ensmp.fr/~leghtas/).  [![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)

## How do I cite `pyEPR`?
 [![DOI](https://zenodo.org/badge/101073856.svg)](https://zenodo.org/badge/latestdoi/101073856)
Use this [bibtex](https://github.com/zlatko-minev/pyEPR/blob/master/pyEPR.bib) for `pyEPR` and for the method use the energy-participation-ratio paper [arXiv:2010.00620](https://arxiv.org/abs/2010.00620). 


[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/zlatko-minev/pyEPR/graphs/commit-activity)

[![Twitter](https://github.frapsoft.com/social/twitter.png)](https://twitter.com/zlatko_minev)
[![Github](https://github.frapsoft.com/social/github.png)](https://github.com/zlatko-minev/)
