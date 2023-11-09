Core module of the pyEPR package
===================
Zlatko Minev & Zaki Leghtas

When using this packed directly, the directory containing this folder should be added to the Python search path.

##### _config_user.py
User specified configuration file. Overwrites `_config_default.py` dictionary.
A user should not edit `_config_default.py` directly. A user should overwrite values in `_config_user.py`

##### core.py
Contains the core analysis and run functions.

##### toolbox
Module  that contains key and utility modules used in pyEPR.
- plotting: useful in visualization and analysis.
- pythonic:  useful pythonic functions
- report: used to plot reports

##### calcs
Module that contains calculation useful sub-modules
- constants: numerical constants
- basic: epr to zpf, and other general functions
- convert:  conversion functions for Lj to Ej, Cj to Ej, and such.
- transmon: transmon related functions

##### ansys.py
Interface functions with Ansys. (A long time ago, this used to be Ansoft HFSS).
Contributed by Phil Rheinhold. Originally part of [pyHFSS](https://github.com/PhilReinhold/pyHFSS).
Updated and modified by Zlatko Minev & Zaki Leghtas.

##### numeic_diag.py
Internal use only. For numerical diagonalization.
Written by Phil Rheinhold.
Updated by Zlatko Minev & Lysander Christakis.
This file is tricky, use caution to modify.