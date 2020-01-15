# This file is part of pyEPR: Energy participation ratio (EPR) design of quantum circuits in python
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

# pylint: disable = wrong-import-position

# Compatibility with python 2.7 and 3
from __future__ import division, print_function, absolute_import

import logging
import warnings
from pathlib import Path
from collections import OrderedDict

try:
    from attrdict import AttrDict as Dict
except (ImportError, ModuleNotFoundError):
    raise ImportError("""Please install python package `AttrDict`.
    AttrDict is in PyPI, so it can be installed directly
    (https://github.com/bcj/AttrDict) using:
        $ pip install attrdict""")

##############################################################################
# Config setup
from ._config_default import get_config
config = get_config()

x=5

##############################################################################
# Set up logging -- only on first loading of module, not on reloading.
logger = logging.getLogger('pyEPR')  # singleton
if not len(logger.handlers):
    from .toolbox._logging import set_up_logger
    set_up_logger(logger)
    del set_up_logger



##############################################################################
#
# CHECK that required packages are available. If not raise log warning.

try:
    import matplotlib as mpl
except (ImportError, ModuleNotFoundError):
    raise ImportError(f"""IMPORT WARNING: Could not find package `matplotlib`.
        Default plotting will not work unless you install it. Please install it.
        {config.internal.error_msg_missing_import}""")

try:
    import pandas as pd
    warnings.filterwarnings('ignore', category=pd.io.pytables.PerformanceWarning)
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
        logger.warning("""IMPORT WARNING: `qutip` package not found.
        Numerical diagonalization will not work. Please install, e.g.:
            $ conda  install -c conda-forge qutip
        %s""", config.internal.error_msg_missing_import)

    try:
        import pythoncom
        del pythoncom
    except (ImportError, ModuleNotFoundError):
        logger.warning("""IMPORT WARNING:
        Python package 'pythoncom' could not be loaded
        It is used in communicting with HFSS on PCs. If you wish to do this, please set it up.
        For Linux, check the HFSS python linux files for the com module used. It is equivalent,
        and can be used just as well.
        %s""", config.internal.error_msg_missing_import)

    try:
        from win32com.client import Dispatch, CDispatch
        del Dispatch, CDispatch
    except (ImportError, ModuleNotFoundError):
        logger.warning("""IMPORT WARNING: Could not load from 'win32com.client'.
        The communication to hfss won't work. If you want to use it, you need to set it up.
        %s""", config.internal.error_msg_missing_import)

    try:
        from pint import UnitRegistry  # units
        del UnitRegistry
    except (ImportError, ModuleNotFoundError):
        logger.error("""IMPORT ERROR:
        Python package 'pint' could not be loaded. It is used in communicting with HFSS. Try:
            $ conda install -c conda-forge pint \n%s""", config.internal.error_msg_missing_import)


##############################################################################
# pyEPR convenience variable and function imports

from . import toolbox
from . import calcs
from . import hfss
from . import core
from . import numeric_diag

from .toolbox.plotting import mpl_dpi
from .hfss import load_ansys_project, get_active_design, get_active_project,\
    HfssProject, CalcObject, parse_units, parse_units_user
from .hfss import release as hfss_release
from .core import Project_Info, pyEPR_HFSS, pyEPR_Analysis
