Core module of the pyEPR package
===================
Zlatko Minev & Zaki Leghtas

When using this packed directly, the directory containing this folder should be added to the Python search path.

##### core.py
Contains the core analysis and run functions.

##### toolbox
Folder that contains key and utility modues used in pyEPR.
- conversions:  conversion functions for Lj to Ej, Cj to Ej, and such.
- plotting: useful in visualizaiton and analysis.
- constants: numerical constants
- pythonic:  useful pythonic functions
- calcs_basic: epr to zpf, and other general functions
- calcs_transmon: transmon related functions

##### hfss.py
Interface functions with Ansoft HFSS.
Contributed by Phil Rheinhold.  Originally part of [pyHFSS](https://github.com/PhilReinhold/pyHFSS).
Updated and modified by Zlatko Minev & Zaki Leghtas.

##### numeic_diag.py
Internal use only. For numerical diagonalizaiton.
Written by Phil Rheinhold.
Updated by Zlatko Minev & Lysander Christakis.
This file is tricky, use caution to modify.