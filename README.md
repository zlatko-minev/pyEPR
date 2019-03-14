Welcome to pyEPR :beers:!
===================
[![Open Source Love](https://badges.frapsoft.com/os/v1/open-source.png?v=103)](https://github.com/zlatko-minev/pyEPR)
[![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/zlatko-minev/pyEPR)
[![star this repo](http://githubbadges.com/star.svg?user=zlatko-minev&repo=pyEPR&style=flat)](https://github.com/zlatko-minev/pyEPR)
[![fork this repo](http://githubbadges.com/fork.svg?user=zlatko-minev&repo=pyEPR&style=flat)](https://github.com/zlatko-minev/pyEPR/fork)

### Automated Python module for the design and quantization of Josephson quantum circuit

Superconducting circuits incorporating non-linear devices, such as Josephson junctions and nanowires, are among the leading platforms for emerging quantum technologies. Promising applications require designing and optimizing circuits with ever-increasing complexity and controlling their dissipative and Hamiltonian parameters to several significant digits. Therefore, there is a growing need for a systematic, simple, and robust approach for precise circuit design, extensible to increased complexity. The energy-participation ratio (EPR) approach presents such an approach to unify the design of dissipation and Hamiltonians around a single concept — the energy participation, a number between zero and one — in a single-step electromagnetic simulation. This markedly reduces the required number of simulations and allows for robust extension to complex systems. The approach is general purpose, derived ab initio, and valid for arbitrary non-linear devices and circuit architectures. Experimental results on a variety of circuit quantum electrodynamics (cQED) devices and architectures, 3D and flip-chip (2.5D), have been demonstrated to exhibit ten percent to percent-level agreement for non-linear coupling and modal Hamiltonian parameters over five-orders of magnitude and across a dozen samples. Here, in this package, all routines of the EPR approach are fully automated. 



### References
* Z.K. Minev, Ph.D. Dissertation, Yale University (2018), see Chapter 4. [arXiv:1902.10355](https://arxiv.org/abs/1902.10355) 
* Z.K. Minev, Z. Leghtas, _et al._ (to appear soon on arXiv) (2019)
<br>

### Showcase

![Intro image](read_me_0.png 'Intro image')

<br><br>

# Contents: 
* [Start here: Using `pyEPR`](#start-here-using-pyepr)
* [Setup and Installation](#installation-and-setup-of-pyepr)
* [HFSS Project Setup for `pyEPR`](#hfss-project-setup-for-pyepr)
* [Troubleshooting `pyEPR`](#troubleshooting-pyepr)
* [Authors and Contributors](#authors-and-contributors)



# Start here: Using `pyEPR`

1. **Fork**  :fork_and_knife: the [``pyEPR top-level repository`` ](https://github.com/zlatko-minev/pyEPR) on GitHub. ([How to fork a GitHub repo?](https://help.github.com/en/articles/fork-a-repo)). Nice to **star** :star: the pyEPR master repository too.
2. **Clone** :point_down: your forked repository locally. ([How to clone a GitHub repo?](https://help.github.com/en/articles/cloning-a-repository)). Setup the `pyEPR` python code by following [Installation and Python Setup](#installation-of-pyepr).
3. **Run** the `pyEPR` script and interface to analyze and optimize your quantum circuits. Examples provided below, more to come soon. 
4. **Cite us** :clap: and enjoy :birthday:! [arXiv:1902.10355](https://arxiv.org/abs/1902.10355) 

#### Start-up example 
The following code illustrates how to perform a complete analysis of a simple two-qubit, one-cavity device in just a few lines of code with `pyEPR`.  In the HFSS file, before running the script, first specify the non-linear junction rectangles and variables (see Sec. pyEPR Project Setup in HFSS). All operations in the eigen analysis and Hamiltonian computation are fully automated. The results are saved, printed, and succinctly plotted. 

```python
from pyEPR import *

# 1.  Project and design. Open link to HFSS controls.
project_info = Project_Info('c:/sims', 
			    project_name = 'two_qubit_one_cavity', # Project file name (string). "None" will get the current active one.
			    design_name  = 'Alice_Bob'             # Design name (string). "None" will get the current active one.
			    )

# 2a. Junctions. Specify junctions in HFSS model
project_info.junctions['jAlice'] = {'Lj_variable':'LJAlice', 'rect':'qubitAlice', 'line': 'alice_line', 'length':0.0001}
project_info.junctions['jBob']   = {'Lj_variable':'LJBob',   'rect':'qubitBob',   'line': 'bob_line',   'length':0.0001}

# 2b. Dissipative elements.
project_info.dissipative.dielectrics_bulk    = ['si_substrate']    # supply names here, there are more options in project_info.dissipative.
project_info.dissipative.dielectric_surfaces = ['interface']   

# 3.  Run analysis
epr_hfss = pyEPR_HFSS(project_info)
epr_hfss.do_EPR_analysis()

# 4.  Hamiltonian analysis
epr      = pyEPR_Analysis(epr_hfss.data_filename)
epr.analyze_all_variations(cos_trunc = 8, fock_trunc = 7)
epr.plot_Hresults()
```

# Installation and setup of `pyEPR`
-------------
If you are starting from scratch, follow the installation guide below to setup a Python 2.7 or 3.x environment ant to fork this repository. To keep up to date with this git, you can use SourceTree, a git-gui manager.

**Recommended procedure.**   <br />

 1. Install a Python 2.7 or 3.x environment.
   * We recommend [Anaconda CE](https://www.continuum.io/downloads), and have tested this installation procedure with Anaconda v4.4 64-bit edition on Windows. Other environments, such as Python XY, or 32 bit have been used in the past, and should work too. We assume you will use Anaconda. First, install Anaconda in "C:\Anaconda2". Note if you are using Python 3, remember to change 2 to 3 in the filenames.
   * Set the Windows System PATH variable. In Control Panel, search for Environment Variables (in System), and open it. In the section System Variables, find the PATH environment variable and select it. Click Edit.  Place`C:\Anaconda2;C:\Anaconda2\Scripts;C:\Anaconda2\Library\bin;` at the beginning. If you have a previous Python installation this step is *very* important, especially to compile the qutip module. You may verity your path using the following command in the Command Prompt (terminal):
      ```sh
      $ echo %PATH%
      ```

 2. Install the required package [pint](http://pint.readthedocs.io/en/latest/). In a terminal window
 ```sh
 conda install -c conda-forge pint
 ```
   * We also use [pandas](http://pandas.pydata.org/). However, as of Aug. 2017, Anaconda includes the pandas package by default, so you do not need to install it manually.
 3. Fork this pyEPR repository on GitHub with your GitHub account. You may clone the fork to your PC and manage it using the [SourceTree](https://www.sourcetreeapp.com/) git-gui manager.
 4. Add the pyEPR repository folder to your python search path.
 5. Edit pyEPR module `config.py`  to set your data-saving directory and other parameters of interest.
 6. **ENJOY! **  :+1:

 
#### Note for Mac/Linux.
Follow the same instructions above. You shouldn't have to install mingw or modify distutils.cfg, since your distribution should come with gcc as the default compiler.

#### Optional package installation
You may also choose to install the optional qutip package for some advanced numerical analysis of the Hamiltonian.
We use [Qutip](http://qutip.org/) to handle quantum objects. Follow the instruction on their website. As of Aug. 2017, qutip is part of conda, and you can use
```sh
conda install qutip
```
If this doesn't work, try  installing from conda forge
```sh
conda install -c conda-forge qutip
```
### Alternative, manual install 
If you wish to install manually, follow the following procedure. Some of this can get a bit tricky at times.
First, you need to install a C compiler, since qutip uses Cython. If you dont have VS9, gcc, or mingw installed, the following works:
```sh
pip install -i https://pypi.anaconda.org/carlkl/simple mingwpy
```
Let anaconda know to use this compiler by creating the file `C:\Anaconda2\Lib\distutils\distutils.cfg` with the following content
```
[build]
compiler = mingw32
[build_ext]
compiler = mingw32
```
Next, let's install qutip. You can choose to use conda intall or pip install, or pull from the git directly  as done here:
```sh
conda install git
pip install git+https://github.com/qutip/qutip.git
```


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


# Troubleshooting pyEPR
---------------------
###### First run: pint error: system='mks' unknown.
Please update to pint version newer than 0.7.2. You may use
```
pip install pint --upgrade
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

# Authors and Contributors
* _Authors:_ [Zlatko Minev](https://www.zlatko-minev.com/) & [Zaki Leghtas](http://cas.ensmp.fr/~leghtas/), with contributions from many friends and colleagues.
* 2015 - present. 
* Contributors: [Phil Rheinhold](https://github.com/PhilReinhold), Lysander Christakis, [Devin Cody](https://github.com/devincody), ...
Original versions of pyHFSS.py and pyNumericalDiagonalization.py contributed by [Phil Rheinhold](https://github.com/PhilReinhold), excellent original [repo](https://github.com/PhilReinhold/pyHFSS).
* Terms of use: Use freely and kindly cite the paper (arXiv link to be posted here) and/or this package.
* How can I contribute? Contact [Z. Minev](zlatko-minev.com) or [Z. Leghtas](http://cas.ensmp.fr/~leghtas/).

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/zlatko-minev/pyEPR/graphs/commit-activity)

[![Twitter](https://github.frapsoft.com/social/twitter.png)](https://twitter.com/zlatko_minev)
[![Github](https://github.frapsoft.com/social/github.png)](https://github.com/zlatko-minev/)
