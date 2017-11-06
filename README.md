Welcome to pyEPR!
===================

### The energy-participation ratio (EPR) approach to the design of quantum Josephson circuits <br> — the Python package.

##### By Zlatko Minev & Zaki Leghtas.

* Join the project!  Contact Zlatko or Zaki if you wish to contribute
* Timeframe: 2015 - present
* Terms of use: Use freely and kindly cite the paper (arXiv link to be posted here) or this package.
* Contributors:  [Zlatko Minev](https://github.com/zlatkom), [Zaki Leghtas](https://github.com/leghtas/), [Phil Rheinhold](https://github.com/PhilReinhold), Lysander Christakis, [Devin Cody](https://github.com/devincody), ...
* pyHFSS and pyNumericalDiagonalization were contributed by [Phil Rheinhold](https://github.com/PhilReinhold). For his excellent, original pyHFSS interface between python and HFSS see [pyHFSS](https://github.com/PhilReinhold/pyHFSS)

# Features
---------------------
TBA


# Installation of pyEPR
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

####Optional package installation
You may also choose to install the optional qutip package for some advanced numerical analysis of the Hamiltonian.
We use [Qutip](http://qutip.org/) to handle quantum objects. Follow the instruction on their website. As of Aug. 2017, qutip is part of conda, and you can use
```sh
	conda install qutip
```
If this doesn't work, try  installing from conda forge
```sh
	conda install -c conda-forge qutip
```

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

# pyEPR Project Setup in HFSS
-------------
#### Eigenmode Design --- How to set up junctions
You may find an advised workflow and some setup tips here.

 1. Define circuit geometry & electromagnetic boundary condition (BC).
   1. Junction rectangles and BC: Create a rectangle for each Josephson junction and give it a good name; e.g., `jAlice` for a qubit named Alice. We recommend 50 x 100 um rectangle for a simple simulation, although orders of magnitude smaller rectangles work as well. Note the length of this junction, you will supply it to pyEPR. Assign a `Lumped RLC` BC on this rectangle surface, with an inductance value given by a local variable, `Lj1` for instance. The name of this variable will also be supplied to the pyEPR.
   2. Over each junction rectangle draw a model `polyline` to define give a sense of the junction current-flow direction. This line should spans the length of the full junction rectangle. Define it using an object coordinate system on the junction rectangle (so that they move together when the geometry is altered). The name of this line will be supplied to the pyEPR module.
 2. Meshing.
   1. Lightly mesh the thin-film metal BC. Lightly mesh the junction rectangles.
 3. Simulation setup
   1. We recommend `mixed order` solutions.


# Troubleshooting pyEPR
---------------------
###### First run: pint error: system='mks' unkown.
Please update to pint version newer than 0.7.2. You may use
```
    pip install pint --upgrade
```

###### COM Error on opening HFSS
Check the project and design file names carefully. Make sure that the file-path doesn't have apostrophes or other bad characters, such as in C:\\Minev's PC\\my:Project.  Check that HFSS hasn't popped up an error dialogue, such as "File locked." Manually open HFSS and the file.

###### COM error on calculation of expression
Either HFSS popped an error dialogue, froze up, or you miss-typed the name of something.

###### HFSS refuses to close
If your script terminates improperly, this can happen. pyHFSS tries to catch termination events and handle them. Your safety should be guaranteed however, if you call `hfss.release()` when you have finished. Use the Task-manager (Activity Monitor on MAC) to kill HFSS if you want.

###### Parametric Sweep Error
When running a parametric sweep in HFSS, make sure you are actually saving the fields for each variation before running pyEPR. This can be done by right-clicking on your ParametricSetup -> properties -> options -> "Save Fields and Mesh".

###### Spyder pops up cmd window with tput.exe
This problem is due to pandas 0.20.1, update to 0.20.3 or better solves this issue.