# This file is part of pyEPR: Energy participation ratio (EPR) design of
# quantum circuits in python
#
#    Copyright (c) 2015-2020 and later, Zlatko K. Minev and Zaki Leghtas
#    All rights reserved.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions are
#    met:
#
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#
#    2. Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#    3. Neither the name of the pyEPR nor the names
#       of its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################
"""
**pyEPR**

Automated Python module for the design and quantization of Josephson quantum circuits

Abstract: Superconducting circuits incorporating non-linear devices, such as Josephson
junctions and nanowires, are among the leading platforms for emerging quantum technologies.
Promising applications require designing and optimizing circuits with ever-increasing
complexity and controlling their dissipative and Hamiltonian parameters to several
significant digits. Therefore, there is a growing need for a systematic, simple, and robust
approach for precise circuit design, extensible to increased complexity.
The energy-participation ratio (EPR) approach presents such an approach to unify the design
of dissipation and Hamiltonians around a single concept — the energy participation, a number
between zero and one — in a single-step electromagnetic simulation. This markedly reduces
the required number of simulations and allows for robust extension to complex systems.
The approach is general purpose, derived ab initio, and valid for arbitrary non-linear
devices and circuit architectures. Experimental results on a variety of circuit quantum
electrodynamics (cQED) devices and architectures, 3D and flip-chip (2.5D), have been
demonstrated to exhibit ten percent to percent-level agreement for non-linear coupling
and modal Hamiltonian parameters over five-orders of magnitude and across a dozen samples.

Here, in this package, all routines of the EPR approach are fully automated.
An interface with ansys is provided.
Automated analysis of lumped and distributed circuits is provided.

@author: Zlatko Minev, Zaki Leghas, ... and the pyEPR team
@site: https://github.com/zlatko-minev/pyEPR
@license: "BSD-3-Clause"
@version: 0.8.5.6
@maintainer: Zlatko K. Minev and  Asaf Diringer
@email: zlatko.minev@aya.yale.edu
@url: https://github.com/zlatko-minev/pyEPR
@status: "Dev-Production"
"""
# pylint: disable= wrong-import-position, invalid-name

# Compatibility with python 2.7 and 3
#from __future__ import division, print_function, absolute_import

import logging
import warnings
from pathlib import Path

from addict import Dict

##############################################################################
# Python header

__author__ = "Zlatko Minev, Zaki Leghas, and the pyEPR team"
__copyright__ = "Copyright 2015-2020, pyEPR team"
__credits__ = [
    "Zlatko Minev", "Zaki Leghtas,", "Phil Rheinhold", "Asaf Diringer",
    "Will Livingston", "Steven Touzard"
]
__license__ = "BSD-3-Clause"
__version__ = "0.8.5.6"
__maintainer__ = "Zlatko K. Minev and  Asaf Diringer"
__email__ = "zlatko.minev@aya.yale.edu"
__url__ = r'https://github.com/zlatko-minev/pyEPR'
__status__ = "Dev-Production"

##############################################################################
# Config setup
from ._config_default import get_config

config = get_config()

##############################################################################
# Set up logging -- only on first loading of module, not on reloading.
logger = logging.getLogger('pyEPR')  # singleton
if not len(logger.handlers):
    from .toolbox._logging import set_up_logger
    set_up_logger(logger)
    del set_up_logger

##############################################################################
#
# Check that required packages are available. If not raise log warning.

try:
    import pandas as pd
    warnings.filterwarnings('ignore',
                            category=pd.io.pytables.PerformanceWarning)
    del pd
except (ImportError, ModuleNotFoundError):
    if config.internal.warn_missing_import:
        logger.warning("IMPORT WARNING: `pandas` python package not found. %s",
                       config.internal.error_msg_missing_import)

# Check for a few usually troublesome packages
if config.internal.warn_missing_import:

    # Check for qutip
    try:
        import qutip
        del qutip
    except (ImportError, ModuleNotFoundError):
        logger.warning(
            """IMPORT WARNING: `qutip` package not found.
        Numerical diagonalization will not work. Please install, e.g.:
            $ conda  install -c conda-forge qutip
        %s""", config.internal.error_msg_missing_import)

    try:
        import pythoncom
        del pythoncom
    except (ImportError, ModuleNotFoundError):
        logger.warning(
            """IMPORT WARNING:
        Python package 'pythoncom' could not be loaded
        It is used in communicating with HFSS on PCs. If you wish to do this, please set it up.
        For Linux, check the HFSS python linux files for the com module used. It is equivalent,
        and can be used just as well.
        %s""", config.internal.error_msg_missing_import)

    try:
        from win32com.client import Dispatch, CDispatch
        del Dispatch
        del CDispatch
    except (ImportError, ModuleNotFoundError):
        logger.warning(
            """IMPORT WARNING: Could not load from 'win32com.client'.
        The communication to hfss won't work. If you want to use it, you need to set it up.
        %s""", config.internal.error_msg_missing_import)

    try:
        import pint  # units
        del pint
    except (ImportError, ModuleNotFoundError):
        logger.error(
            """IMPORT ERROR:
        Python package 'pint' could not be loaded. It is used in communicating with HFSS. Try:
            $ conda install -c conda-forge pint \n%s""",
            config.internal.error_msg_missing_import)

# remove unused
del Path, warnings, logging

##############################################################################
# pyEPR convenience variable and function imports

from . import toolbox
from . import calcs
from . import ansys
from . import core

from .ansys import parse_units, parse_units_user, parse_entry
from .core import ProjectInfo, DistributedAnalysis, QuantumAnalysis,\
                  Project_Info, pyEPR_HFSSAnalysis, pyEPR_Analysis # names to be deprecated

__all__ = [
    'logger',
    'config',
    'toolbox',
    'calcs',
    'ansys',
    'core',
    'ProjectInfo',
    'DistributedAnalysis',
    'QuantumAnalysis',
    'Project_Info',
    'pyEPR_HFSSAnalysis',
    'pyEPR_Analysis',  # names to be deprecated
    'parse_units',
    'parse_units_user',
    'parse_entry'
]

# TODO: Add "about" method. Add to tutorial
