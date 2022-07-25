'''
pyEPR.ansys
    2014-present

Purpose:
    Handles Ansys interaction and control from version 2014 onward.
    Tested most extensively with V2016 and V2019R3.

@authors:
    Originally contributed by Phil Reinhold.
    Developed further by Zlatko Minev, Zaki Leghtas, and the pyEPR team.
    For the base version of hfss.py, see https://github.com/PhilReinhold/pyHFSS
'''

# Python 2.7 and 3 compatibility
from __future__ import (division, print_function)

from typing import List

import atexit
import os
import re
import signal
import tempfile
import time
import types
from collections.abc import Iterable
from copy import copy
from numbers import Number
from pathlib import Path

import numpy as np
import pandas as pd
from sympy.parsing import sympy_parser
import io

from . import logger

# Handle a  few usually troublesome to import packages, which the use may not have
# installed yet
try:
    import pythoncom
except (ImportError, ModuleNotFoundError):
    pass #raise NameError ("pythoncom module not installed. Please install.")

try:
    # TODO: Replace `win32com` with Linux compatible package.
    # See Ansys python files in IronPython internal.
    from win32com.client import Dispatch, CDispatch
except (ImportError, ModuleNotFoundError):
    pass #raise NameError ("win32com module not installed. Please install.")

try:
    from pint import UnitRegistry
    ureg = UnitRegistry()
    Q = ureg.Quantity
except (ImportError, ModuleNotFoundError):
    pass # raise NameError ("Pint module not installed. Please install.")


##############################################################################
###

BASIS_ORDER = {
    "Zero Order": 0,
    "First Order": 1,
    "Second Order": 2,
    "Mixed Order": -1
}

# UNITS
# LENGTH_UNIT         --- HFSS UNITS
# #Assumed default input units for ansys hfss
LENGTH_UNIT = 'meter'
# LENGTH_UNIT_ASSUMED --- USER UNITS
# if a user inputs a blank number with no units in `parse_fix`,
# we can assume the following using
LENGTH_UNIT_ASSUMED = 'mm'


def simplify_arith_expr(expr):
    try:
        out = repr(sympy_parser.parse_expr(str(expr)))
        return out
    except:
        print("Couldn't parse", expr)
        raise


def increment_name(base, existing):
    if not base in existing:
        return base
    n = 1

    def make_name():
        return base + str(n)

    while make_name() in existing:
        n += 1
    return make_name()


def extract_value_unit(expr, units):
    """
    :type expr: str
    :type units: str
    :return: float
    """
    try:
        return Q(expr).to(units).magnitude
    except Exception:
        try:
            return float(expr)
        except Exception:
            return expr


def extract_value_dim(expr):
    """
    type expr: str
    """
    return str(Q(expr).dimensionality)


def parse_entry(entry, convert_to_unit=LENGTH_UNIT):
    '''
    Should take a list of tuple of list... of int, float or str...
    For iterables, returns lists
    '''
    if not isinstance(entry, list) and not isinstance(entry, tuple):
        return extract_value_unit(entry, convert_to_unit)
    else:
        entries = entry
        _entry = []
        for entry in entries:
            _entry.append(parse_entry(entry, convert_to_unit=convert_to_unit))
        return _entry


def fix_units(x, unit_assumed=None):
    '''
    Convert all numbers to string and append the assumed units if needed.
    For an iterable, returns a list
    '''
    unit_assumed = LENGTH_UNIT_ASSUMED if unit_assumed is None else unit_assumed
    if isinstance(x, str):
        # Check if there are already units defined, assume of form 2.46mm  or 2.0 or 4.
        if x[-1].isdigit() or x[-1] == '.':  # number
            return x + unit_assumed
        else:  # units are already applied
            return x

    elif isinstance(x, Number):
        return fix_units(str(x) + unit_assumed, unit_assumed=unit_assumed)

    elif isinstance(x, Iterable):  # hasattr(x, '__iter__'):
        return [fix_units(y, unit_assumed=unit_assumed) for y in x]
    else:
        return x


def parse_units(x):
    '''
    Convert number, string, and lists/arrays/tuples to numbers scaled
    in HFSS units.

    Converts to                  LENGTH_UNIT = meters  [HFSS UNITS]
    Assumes input units  LENGTH_UNIT_ASSUMED = mm      [USER UNITS]

    [USER UNITS] ----> [HFSS UNITS]
    '''
    return parse_entry(fix_units(x))


def unparse_units(x):
    '''
        Undo effect of parse_unit.

        Converts to     LENGTH_UNIT_ASSUMED = mm     [USER UNITS]
        Assumes input units     LENGTH_UNIT = meters [HFSS UNITS]

        [HFSS UNITS] ----> [USER UNITS]
    '''
    return parse_entry(fix_units(x, unit_assumed=LENGTH_UNIT),
                       LENGTH_UNIT_ASSUMED)


def parse_units_user(x):
    '''
        Convert from user assumed units to user assumed units
        [USER UNITS] ----> [USER UNITS]
    '''
    return parse_entry(fix_units(x, LENGTH_UNIT_ASSUMED), LENGTH_UNIT_ASSUMED)


class VariableString(str):
    def __add__(self, other):
        return var("(%s) + (%s)" % (self, other))

    def __radd__(self, other):
        return var("(%s) + (%s)" % (other, self))

    def __sub__(self, other):
        return var("(%s) - (%s)" % (self, other))

    def __rsub__(self, other):
        return var("(%s) - (%s)" % (other, self))

    def __mul__(self, other):
        return var("(%s) * (%s)" % (self, other))

    def __rmul__(self, other):
        return var("(%s) * (%s)" % (other, self))

    def __div__(self, other):
        return var("(%s) / (%s)" % (self, other))

    def __rdiv__(self, other):
        return var("(%s) / (%s)" % (other, self))

    def __truediv__(self, other):
        return var("(%s) / (%s)" % (self, other))

    def __rtruediv__(self, other):
        return var("(%s) / (%s)" % (other, self))

    def __pow__(self, other):
        return var("(%s) ^ (%s)" % (self, other))

    def __rpow__(self, other):
        return var("(%s) ^ (%s)" % (other, self))

    def __neg__(self):
        return var("-(%s)" % self)

    def __abs__(self):
        return var("abs(%s)" % self)


def var(x):
    if isinstance(x, str):
        return VariableString(simplify_arith_expr(x))
    return x


_release_fns = []


def _add_release_fn(fn):
    global _release_fns
    _release_fns.append(fn)
    atexit.register(fn)
    signal.signal(signal.SIGTERM, fn)
    signal.signal(signal.SIGABRT, fn)


def release():
    '''
    Release COM connection to Ansys.
    '''
    global _release_fns
    for fn in _release_fns:
        fn()
    time.sleep(0.1)

    # Note that _GetInterfaceCount is a member
    refcount = pythoncom._GetInterfaceCount()  # pylint: disable=no-member

    if refcount > 0:
        print("Warning! %d COM references still alive" % (refcount))
        print("Ansys will likely refuse to shut down")


class COMWrapper(object):
    def __init__(self):
        _add_release_fn(self.release)

    def release(self):
        for k, v in self.__dict__.items():
            if isinstance(v, CDispatch):
                setattr(self, k, None)


class HfssPropertyObject(COMWrapper):
    prop_holder = None
    prop_tab = None
    prop_server = None


def make_str_prop(name, prop_tab=None, prop_server=None):
    return make_prop(name, prop_tab=prop_tab, prop_server=prop_server)


def make_int_prop(name, prop_tab=None, prop_server=None):
    return make_prop(name,
                     prop_tab=prop_tab,
                     prop_server=prop_server,
                     prop_args=["MustBeInt:=", True])


def make_float_prop(name, prop_tab=None, prop_server=None):
    return make_prop(name,
                     prop_tab=prop_tab,
                     prop_server=prop_server,
                     prop_args=["MustBeInt:=", False])


def make_prop(name, prop_tab=None, prop_server=None, prop_args=None):
    def set_prop(self,
                 value,
                 prop_tab=prop_tab,
                 prop_server=prop_server,
                 prop_args=prop_args):
        prop_tab = self.prop_tab if prop_tab is None else prop_tab
        prop_server = self.prop_server if prop_server is None else prop_server
        if isinstance(prop_tab, types.FunctionType):
            prop_tab = prop_tab(self)
        if isinstance(prop_server, types.FunctionType):
            prop_server = prop_server(self)
        if prop_args is None:
            prop_args = []
        self.prop_holder.ChangeProperty([
            "NAME:AllTabs",
            [
                "NAME:" + prop_tab, ["NAME:PropServers", prop_server],
                [
                    "NAME:ChangedProps",
                    ["NAME:" + name, "Value:=", value] + prop_args
                ]
            ]
        ])

    def get_prop(self, prop_tab=prop_tab, prop_server=prop_server):
        prop_tab = self.prop_tab if prop_tab is None else prop_tab
        prop_server = self.prop_server if prop_server is None else prop_server
        if isinstance(prop_tab, types.FunctionType):
            prop_tab = prop_tab(self)
        if isinstance(prop_server, types.FunctionType):
            prop_server = prop_server(self)
        return self.prop_holder.GetPropertyValue(prop_tab, prop_server, name)

    return property(get_prop, set_prop)


def set_property(prop_holder,
                 prop_tab,
                 prop_server,
                 name,
                 value,
                 prop_args=None):
    '''
    More general non obj oriented, functional version
    prop_args = [] by default
    '''
    if not isinstance(prop_server, list):
        prop_server = [prop_server]
    return prop_holder.ChangeProperty([
        "NAME:AllTabs",
        [
            "NAME:" + prop_tab, ["NAME:PropServers", *prop_server],
            [
                "NAME:ChangedProps",
                ["NAME:" + name, "Value:=", value] + (prop_args or [])
            ]
        ]
    ])


class HfssApp(COMWrapper):
    def __init__(self, ProgID='AnsoftHfss.HfssScriptInterface'):
        '''
         Connect to IDispatch-based COM object.
             Parameter is the ProgID or CLSID of the COM object.
             This is found in the regkey.

         Version changes for Ansys HFSS for the main object
             v2016 - 'Ansoft.ElectronicsDesktop'
             v2017 and subsequent - 'AnsoftHfss.HfssScriptInterface'

        '''
        super(HfssApp, self).__init__()
        self._app = Dispatch(ProgID)

    def get_app_desktop(self):
        return HfssDesktop(self, self._app.GetAppDesktop())
        # in v2016, there is also getApp - which can be called with HFSS


class HfssDesktop(COMWrapper):
    def __init__(self, app, desktop):
        """
        :type app: HfssApp
        :type desktop: Dispatch
        """
        super(HfssDesktop, self).__init__()
        self.parent = app
        self._desktop = desktop

        # ansys version, needed to check for command changes,
        # since some commands have changed over the years
        self.version = self.get_version()

    def close_all_windows(self):
        self._desktop.CloseAllWindows()

    def project_count(self):
        count = len(self._desktop.GetProjects())
        return count

    def get_active_project(self):
        return HfssProject(self, self._desktop.GetActiveProject())

    def get_projects(self):
        return [HfssProject(self, p) for p in self._desktop.GetProjects()]

    def get_project_names(self):
        return self._desktop.GetProjectList()

    def get_messages(self, project_name="", design_name="", level=0):
        """Use:  Collects the messages from a specified project and design.
        Syntax:              GetMessages <ProjectName>, <DesignName>, <SeverityName>
        Return Value:    A simple array of strings.

        Parameters:
        <ProjectName>
            Type:<string>
            Name of the project for which to collect messages.
            An incorrect project name results in no messages (design is ignored)
            An empty project name results in all messages (design is ignored)

        <DesignName>
            Type: <string>
            Name of the design in the named project for which to collect messages
            An incorrect design name results in no messages for the named project
            An empty design name results in all messages for the named project

        <SeverityName>
            Type: <integer>
            Severity is 0-3, and is tied in to info/warning/error/fatal types as follows:
                0 is info and above
                1 is warning and above
                2 is error and fatal
                3 is fatal only (rarely used)
        """
        return self._desktop.GetMessages(project_name, design_name, level)

    def get_version(self):
        return self._desktop.GetVersion()

    def new_project(self):
        return HfssProject(self, self._desktop.NewProject())

    def open_project(self, path):
        ''' returns error if already open '''
        return HfssProject(self, self._desktop.OpenProject(path))

    def set_active_project(self, name):
        self._desktop.SetActiveProject(name)

    @property
    def project_directory(self):
        return self._desktop.GetProjectDirectory()

    @project_directory.setter
    def project_directory(self, path):
        self._desktop.SetProjectDirectory(path)

    @property
    def library_directory(self):
        return self._desktop.GetLibraryDirectory()

    @library_directory.setter
    def library_directory(self, path):
        self._desktop.SetLibraryDirectory(path)

    @property
    def temp_directory(self):
        return self._desktop.GetTempDirectory()

    @temp_directory.setter
    def temp_directory(self, path):
        self._desktop.SetTempDirectory(path)


class HfssProject(COMWrapper):
    def __init__(self, desktop, project):
        """
        :type desktop: HfssDesktop
        :type project: Dispatch
        """
        super(HfssProject, self).__init__()
        self.parent = desktop
        self._project = project
        #self.name = project.GetName()
        self._ansys_version = self.parent.version

    def close(self):
        self._project.Close()

    def make_active(self):
        self.parent.set_active_project(self.name)

    def get_designs(self):
        return [HfssDesign(self, d) for d in self._project.GetDesigns()]

    def get_design_names(self):
        return [d.GetName() for d in self._project.GetDesigns()]

    def save(self, path=None):
        if path is None:
            self._project.Save()
        else:
            self._project.SaveAs(path, True)

    def simulate_all(self):
        self._project.SimulateAll()

    def import_dataset(self, path):
        self._project.ImportDataset(path)

    def rename_design(self, design, rename):
        if design in self.get_designs():
            design.rename_design(design.name, rename)
        else:
            raise ValueError('%s design does not exist' % design.name)

    def duplicate_design(self, target, source):
        src_design = self.get_design(source)
        return src_design.duplicate(name=target)

    def get_variable_names(self):
        return [VariableString(s) for s in self._project.GetVariables()]

    def get_variables(self):
        """ Returns the project variables only, which start with $. These are global variables. """
        return {
            VariableString(s): self.get_variable_value(s)
            for s in self._project.GetVariables()
        }

    def get_variable_value(self, name):
        return self._project.GetVariableValue(name)

    def create_variable(self, name, value):
        self._project.ChangeProperty([
            "NAME:AllTabs",
            [
                "NAME:ProjectVariableTab",
                ["NAME:PropServers", "ProjectVariables"],
                [
                    "Name:NewProps",
                    [
                        "NAME:" + name, "PropType:=", "VariableProp",
                        "UserDef:=", True, "Value:=", value
                    ]
                ]
            ]
        ])

    def set_variable(self, name, value):
        if name not in self._project.GetVariables():
            self.create_variable(name, value)
        else:
            self._project.SetVariableValue(name, value)
        return VariableString(name)

    def get_path(self):
        if self._project:
            return self._project.GetPath()
        else:
            raise Exception('''Error: HFSS Project does not have a path.
        Either there is no HFSS project open, or it is not saved.''')

    def new_design(self, design_name, solution_type, design_type="HFSS"):
        design_name_int = increment_name(
            design_name, [d.GetName() for d in self._project.GetDesigns()])
        return HfssDesign(
            self,
            self._project.InsertDesign(design_type, design_name_int,
                                       solution_type, ""))

    def get_design(self, name):
        return HfssDesign(self, self._project.GetDesign(name))

    def get_active_design(self):
        d = self._project.GetActiveDesign()
        if d is None:
            raise EnvironmentError("No Design Active")
        return HfssDesign(self, d)

    def new_dm_design(self, name: str):
        """Create a new driven model design

        Args:
            name (str): Name of driven modal design
        """
        return self.new_design(name, "DrivenModal")

    def new_em_design(self, name: str):
        """Create a new eigenmode design

        Args:
            name (str): Name of eigenmode design
        """
        return self.new_design(name, "Eigenmode")

    def new_q3d_design(self, name: str):
        """Create a new Q3D design.
        Args:
            name (str): Name of Q3D design
        """
        return self.new_design(name, "Q3D", "Q3D Extractor")

    @property  # v2016
    def name(self):
        return self._project.GetName()


class HfssDesign(COMWrapper):
    def __init__(self, project, design):
        super(HfssDesign, self).__init__()
        self.parent = project
        self._design = design
        self.name = design.GetName()
        self._ansys_version = self.parent._ansys_version

        try:
            # This function does not exist if the design is not HFSS
            self.solution_type = design.GetSolutionType()
        except Exception as e:
            logger.debug(
                f'Exception occurred at design.GetSolutionType() {e}. Assuming Q3D design'
            )
            self.solution_type = 'Q3D'

        if design is None:
            return
        self._setup_module = design.GetModule("AnalysisSetup")
        self._solutions = design.GetModule("Solutions")
        self._fields_calc = design.GetModule("FieldsReporter")
        self._output = design.GetModule("OutputVariable")
        self._boundaries = design.GetModule("BoundarySetup")
        self._reporter = design.GetModule("ReportSetup")
        self._modeler = design.SetActiveEditor("3D Modeler")
        self._optimetrics = design.GetModule("Optimetrics")
        self._mesh = design.GetModule("MeshSetup")
        self.modeler = HfssModeler(self, self._modeler, self._boundaries,
                                   self._mesh)
        self.optimetrics = Optimetrics(self)

    def add_message(self, message: str, severity: int = 0):
        """
        Add a message to HFSS log with severity and context to message window.

        Keyword Args:
            severity (int) : 0 = Informational, 1 = Warning, 2 = Error, 3 = Fatal..
        """
        project = self.parent
        desktop = project.parent
        oDesktop = desktop._desktop
        oDesktop.AddMessage(project.name, self.name, severity, message)

    def save_screenshot(self, path: str = None, show: bool = True):
        if not path:
            path = Path().absolute() / 'ansys.png'  # TODO find better
        self._modeler.ExportModelImageToFile(
            str(path),
            0,
            0,  # can be 0 For the default, use 0, 0. For higher resolution, set desired <width> and <height>, for example for 8k export as: 7680, 4320.
            [
                "NAME:SaveImageParams", "ShowAxis:=", "True", "ShowGrid:=",
                "True", "ShowRuler:=", "True", "ShowRegion:=", "Default",
                "Selections:=", "", "Orientation:=", ""
            ])

        if show:
            from IPython.display import display, Image
            display(Image(str(path)))

        return path

    def rename_design(self, name):
        old_name = self._design.GetName()
        self._design.RenameDesignInstance(old_name, name)

    def copy_to_project(self, project):
        project.make_active()
        project._project.CopyDesign(self.name)
        project._project.Paste()
        return project.get_active_design()

    def duplicate(self, name=None):
        dup = self.copy_to_project(self.parent)
        if name is not None:
            dup.rename_design(name)
        return dup

    def get_setup_names(self):
        return self._setup_module.GetSetups()

    def get_setup(self, name=None):
        """
        :rtype: HfssSetup
        """
        setups = self.get_setup_names()
        if not setups:
            raise EnvironmentError(" *** No Setups Present ***")
        if name is None:
            name = setups[0]
        elif name not in setups:
            raise EnvironmentError("Setup {} not found: {}".format(
                name, setups))

        if self.solution_type == "Eigenmode":
            return HfssEMSetup(self, name)
        elif self.solution_type == "DrivenModal":
            return HfssDMSetup(self, name)
        elif self.solution_type == "DrivenTerminal":
            return HfssDTSetup(self, name)
        elif self.solution_type == "Q3D":
            return AnsysQ3DSetup(self, name)

    def create_q3d_setup(self,
                         freq_ghz=5.,
                         name="Setup",
                         save_fields=False,
                         enabled=True,
                         max_passes=15,
                         min_passes=2,
                         min_converged_passes=2,
                         percent_error=0.5,
                         percent_refinement=30,
                         auto_increase_solution_order=True,
                         solution_order="High",
                         solver_type='Iterative'):
        name = increment_name(name, self.get_setup_names())
        self._setup_module.InsertSetup("Matrix", [
            f"NAME:{name}", "AdaptiveFreq:=", f"{freq_ghz}GHz", "SaveFields:=",
            save_fields, "Enabled:=", enabled,
            [
                "NAME:Cap", "MaxPass:=", max_passes, "MinPass:=", min_passes,
                "MinConvPass:=", min_converged_passes, "PerError:=",
                percent_error, "PerRefine:=", percent_refinement,
                "AutoIncreaseSolutionOrder:=", auto_increase_solution_order,
                "SolutionOrder:=", solution_order, "Solver Type:=", solver_type
            ]
        ])
        return AnsysQ3DSetup(self, name)

    def create_dm_setup(self,
                        freq_ghz=1,
                        name="Setup",
                        max_delta_s=0.1,
                        max_passes=10,
                        min_passes=1,
                        min_converged=1,
                        pct_refinement=30,
                        basis_order=-1):

        name = increment_name(name, self.get_setup_names())
        self._setup_module.InsertSetup("HfssDriven", [
            "NAME:" + name, "Frequency:=",
            str(freq_ghz) + "GHz", "MaxDeltaS:=", max_delta_s,
            "MaximumPasses:=", max_passes, "MinimumPasses:=", min_passes,
            "MinimumConvergedPasses:=", min_converged, "PercentRefinement:=",
            pct_refinement, "IsEnabled:=", True, "BasisOrder:=", basis_order
        ])
        return HfssDMSetup(self, name)

    def create_dt_setup(self,
                        freq_ghz=1,
                        name="Setup",
                        max_delta_s=0.1,
                        max_passes=10,
                        min_passes=1,
                        min_converged=1,
                        pct_refinement=30,
                        basis_order=-1):

        name = increment_name(name, self.get_setup_names())
        self._setup_module.InsertSetup("HfssDriven", [
            "NAME:" + name, "Frequency:=",
            str(freq_ghz) + "GHz", "MaxDeltaS:=", max_delta_s,
            "MaximumPasses:=", max_passes, "MinimumPasses:=", min_passes,
            "MinimumConvergedPasses:=", min_converged, "PercentRefinement:=",
            pct_refinement, "IsEnabled:=", True, "BasisOrder:=", basis_order
        ])
        return HfssDTSetup(self, name)

    def create_em_setup(self,
                        name="Setup",
                        min_freq_ghz=1,
                        n_modes=1,
                        max_delta_f=0.1,
                        max_passes=10,
                        min_passes=1,
                        min_converged=1,
                        pct_refinement=30,
                        basis_order=-1):

        name = increment_name(name, self.get_setup_names())
        self._setup_module.InsertSetup("HfssEigen", [
            "NAME:" + name, "MinimumFrequency:=",
            str(min_freq_ghz) + "GHz", "NumModes:=", n_modes, "MaxDeltaFreq:=",
            max_delta_f, "ConvergeOnRealFreq:=", True, "MaximumPasses:=",
            max_passes, "MinimumPasses:=", min_passes,
            "MinimumConvergedPasses:=", min_converged, "PercentRefinement:=",
            pct_refinement, "IsEnabled:=", True, "BasisOrder:=", basis_order
        ])
        return HfssEMSetup(self, name)

    def delete_setup(self, name):
        if name in self.get_setup_names():
            self._setup_module.DeleteSetups(name)

    def delete_full_variation(self,
                              DesignVariationKey="All",
                              del_linked_data=False):
        """
        DeleteFullVariation
        Use:                   Use to selectively make deletions or delete all solution data.
        Command:         HFSS>Results>Clean Up Solutions...
        Syntax:              DeleteFullVariation Array(<parameters>), boolean
        Parameters:      All | <DataSpecifierArray>
                        If, All, all data of existing variations is deleted.
                        Array(<DesignVariationKey>, )
                        <DesignVariationKey>
                            Type: <string>
                            Design variation string.
                        <Boolean>
                        Type: boolean
                        Whether to also delete linked data.
        """
        self._design.DeleteFullVariation("All", False)

    def get_nominal_variation(self):
        """
        Use: Gets the nominal variation string
        Return Value: Returns a string representing the nominal variation
        Returns string such as "Height='0.06mm' Lj='13.5nH'"
        """
        return self._design.GetNominalVariation()

    def create_variable(self, name, value, postprocessing=False):
        if postprocessing == True:
            variableprop = "PostProcessingVariableProp"
        else:
            variableprop = "VariableProp"

        self._design.ChangeProperty([
            "NAME:AllTabs",
            [
                "NAME:LocalVariableTab",
                ["NAME:PropServers", "LocalVariables"],
                [
                    "Name:NewProps",
                    [
                        "NAME:" + name, "PropType:=", variableprop,
                        "UserDef:=", True, "Value:=", value
                    ]
                ]
            ]
        ])

    def _variation_string_to_variable_list(self,
                                           variation_string: str,
                                           for_prop_server=True):
        """Example:
            Takes
                "Cj='2fF' Lj='13.5nH'"
            for for_prop_server=True into
                [['NAME:Cj', 'Value:=', '2fF'], ['NAME:Lj', 'Value:=', '13.5nH']]
            or for for_prop_server=False into
                [['Cj', '2fF'], ['Lj', '13.5nH']]
        """
        s = variation_string
        s = s.split(' ')
        s = [s1.strip().strip("''").split("='") for s1 in s]

        if for_prop_server:
            local, project = [], []

            for arr in s:
                to_add = [f'NAME:{arr[0]}', "Value:=", arr[1]]
                if arr[0][0] == '$':
                    project += [to_add]  # global variable
                else:
                    local += [to_add]  # local variable

            return local, project

        else:
            return s

    def set_variables(self, variation_string: str):
        """
        Set all variables to match a solved variation string.

        Args:
            variation_string (str) :  Variation string such as
                "Cj='2fF' Lj='13.5nH'"
        """
        assert isinstance(variation_string, str)

        content = ["NAME:ChangedProps"]
        local, project = self._variation_string_to_variable_list(
            variation_string)
        #print('\nlocal=', local, '\nproject=', project)

        if len(project) > 0:
            self._design.ChangeProperty([
                "NAME:AllTabs",
                [
                    "NAME:ProjectVariableTab",
                    ["NAME:PropServers", "ProjectVariables"], content + project
                ]
            ])

        if len(local) > 0:
            self._design.ChangeProperty([
                "NAME:AllTabs",
                [
                    "NAME:LocalVariableTab",
                    ["NAME:PropServers", "LocalVariables"], content + local
                ]
            ])

    def set_variable(self, name: str, value: str, postprocessing=False):
        """Warning: THis is case sensitive,

        Arguments:
            name {str} -- Name of variable to set, such as 'Lj_1'.
                          This is not the same as as 'LJ_1'.
                          You must use the same casing.
            value {str} -- Value, such as '10nH'

        Keyword Arguments:
            postprocessing {bool} -- Postprocessing variable only or not.
                          (default: {False})

        Returns:
            VariableString
        """
        # TODO: check if variable does not exist and quit if it doesn't?
        if name not in self.get_variable_names():
            self.create_variable(name, value, postprocessing=postprocessing)
        else:
            self._design.SetVariableValue(name, value)

        return VariableString(name)

    def get_variable_value(self, name):
        """ Can only access the design variables, i.e., the local ones
            Cannot access the project (global) variables, which start with $. """
        return self._design.GetVariableValue(name)

    def get_variable_names(self):
        """ Returns the local design variables.
            Does not return the project (global) variables, which start with $. """
        return [
            VariableString(s) for s in self._design.GetVariables() +
            self._design.GetPostProcessingVariables()
        ]

    def get_variables(self):
        """ Returns dictionary of local design variables and their values.
            Does not return the project (global) variables and their values,
            whose names start with $. """
        local_variables = self._design.GetVariables(
        ) + self._design.GetPostProcessingVariables()
        return {lv: self.get_variable_value(lv) for lv in local_variables}

    def copy_design_variables(self, source_design):
        ''' does not check that variables are all present '''

        # don't care about values
        source_variables = source_design.get_variables()

        for name, value in source_variables.items():
            self.set_variable(name, value)

    def get_excitations(self):
        self._boundaries.GetExcitations()

    def _evaluate_variable_expression(self, expr, units):
        """
        :type expr: str
        :type units: str
        :return: float
        """
        try:
            sexp = sympy_parser.parse_expr(expr)
        except SyntaxError:
            return Q(expr).to(units).magnitude

        sub_exprs = {
            fs: self.get_variable_value(fs.name)
            for fs in sexp.free_symbols
        }

        return float(
            sexp.subs({
                fs: self._evaluate_variable_expression(e, units)
                for fs, e in sub_exprs.items()
            }))

    def eval_expr(self, expr, units="mm"):
        return str(self._evaluate_variable_expression(expr, units)) + units

    def Clear_Field_Clac_Stack(self):
        self._fields_calc.CalcStack("Clear")

    def clean_up_solutions(self):
        self._design.DeleteFullVariation('All', True)  # Delete existing solutions


class HfssSetup(HfssPropertyObject):
    prop_tab = "HfssTab"
    passes = make_int_prop("Passes")  # see EditSetup
    n_modes = make_int_prop("Modes")
    pct_refinement = make_float_prop("Percent Refinement")
    delta_f = make_float_prop("Delta F")
    min_freq = make_float_prop("Min Freq")
    basis_order = make_str_prop("Basis Order")

    def __init__(self, design, setup: str):
        """
        :type design: HfssDesign
        :type setup: Dispatch

        :COM Scripting Help: "Analysis Setup Module Script Commands"

        Get properties:
            setup.parent._design.GetProperties("HfssTab",'AnalysisSetup:Setup1')
        """
        super(HfssSetup, self).__init__()
        self.parent = design
        self.prop_holder = design._design
        self._setup_module = design._setup_module
        self._reporter = design._reporter
        self._solutions = design._solutions
        self.name = setup
        self.solution_name = setup + " : LastAdaptive"
        #self.solution_name_pass = setup + " : AdaptivePass"
        self.prop_server = "AnalysisSetup:" + setup
        self.expression_cache_items = []
        self._ansys_version = self.parent._ansys_version

    def analyze(self, name=None):
        '''
        Use:             Solves a single solution setup and all of its frequency sweeps.
        Command:         Right-click a solution setup in the project tree, and then click Analyze
                         on the shortcut menu.
        Syntax:          Analyze(<SetupName>)
        Parameters:      <setupName>
        Return Value:    None
        -----------------------------------------------------

        Will block the until the analysis is completely done.
        Will raise a com_error if analysis is aborted in HFSS.
        '''
        if name is None:
            name = self.name
        logger.info(f'Analyzing setup {name}')
        return self.parent._design.Analyze(name)

    def solve(self, name=None):
        '''
        Use:             Performs a blocking simulation.
                         The next script command will not be executed
                         until the simulation is complete.

        Command:         HFSS>Analyze
        Syntax:          Solve <SetupNameArray>
        Return Value:   Type: <int>
                        -1: simulation error
                        0: normal completion
        Parameters:      <SetupNameArray>: Array(<SetupName>, <SetupName>, ...)
           <SetupName>
        Type: <string>
        Name of the solution setup to solve.
        Example:
            return_status = oDesign.Solve Array("Setup1", "Setup2")
        -----------------------------------------------------

        HFSS abort: still returns 0 , since termination by user.

        '''
        if name is None:
            name = self.name
        return self.parent._design.Solve(name)

    def insert_sweep(self,
                     start_ghz,
                     stop_ghz,
                     count=None,
                     step_ghz=None,
                     name="Sweep",
                     type="Fast",
                     save_fields=False):

        if not type in ['Fast', 'Interpolating', 'Discrete']:
            logger.error(
                "insert_sweep: Error type was not in  ['Fast', 'Interpolating', 'Discrete']"
            )

        name = increment_name(name, self.get_sweep_names())
        params = [
            "NAME:" + name,
            "IsEnabled:=",     True,
            "Type:=",          type,
            "SaveFields:=",    save_fields,
            "SaveRadFields:=", False,
            # "GenerateFieldsForAllFreqs:="
            "ExtrapToDC:=",    False,
        ]

        # not sure when exactly this changed between 2016 and 2019
        if self._ansys_version >= '2019':
            if count:
                params.extend([
                    "RangeType:=", 'LinearCount', "RangeStart:=",
                    f"{start_ghz:f}GHz", "RangeEnd:=", f"{stop_ghz:f}GHz",
                    "RangeCount:=", count
                ])
            if step_ghz:
                params.extend([
                    "RangeType:=", 'LinearStep', "RangeStart:=",
                    f"{start_ghz:f}GHz", "RangeEnd:=", f"{stop_ghz:f}GHz",
                    "RangeStep:=", step_ghz
                ])

            if (count and step_ghz) or ((not count) and (not step_ghz)):
                logger.error(
                    'ERROR: you should provide either step_ghz or count \
                    when inserting an HFSS driven model freq sweep. \
                    YOu either provided both or neither! See insert_sweep.')
        else:
            params.extend([
                "StartValue:=",
                "%fGHz" % start_ghz, "StopValue:=",
                "%fGHz" % stop_ghz
            ])
            if step_ghz is not None:
                params.extend([
                    "SetupType:=", "LinearSetup", "StepSize:=",
                    "%fGHz" % step_ghz
                ])
            else:
                params.extend(["SetupType:=", "LinearCount", "Count:=", count])

        self._setup_module.InsertFrequencySweep(self.name, params)

        return HfssFrequencySweep(self, name)

    def delete_sweep(self, name):
        self._setup_module.DeleteSweep(self.name, name)


#    def add_fields_convergence_expr(self, expr, pct_delta, phase=0):
#        """note: because of hfss idiocy, you must call "commit_convergence_exprs"
#         after adding all exprs"""
#        assert isinstance(expr, NamedCalcObject)
#        self.expression_cache_items.append(
#            ["NAME:CacheItem",
#             "Title:=", expr.name+"_conv",
#             "Expression:=", expr.name,
#             "Intrinsics:=", "Phase='{}deg'".format(phase),
#             "IsConvergence:=", True,
#             "UseRelativeConvergence:=", 1,
#             "MaxConvergenceDelta:=", pct_delta,
#             "MaxConvergeValue:=", "0.05",
#             "ReportType:=", "Fields",
#             ["NAME:ExpressionContext"]])

#    def commit_convergence_exprs(self):
#        """note: this will eliminate any convergence expressions not added
#           through this interface"""
#        args = [
#            "NAME:"+self.name,
#            ["NAME:ExpressionCache", self.expression_cache_items]
#        ]
#        self._setup_module.EditSetup(self.name, args)

    def get_sweep_names(self):
        return self._setup_module.GetSweeps(self.name)

    def get_sweep(self, name=None):
        sweeps = self.get_sweep_names()
        if not sweeps:
            raise EnvironmentError("No Sweeps Present")
        if name is None:
            name = sweeps[0]
        elif name not in sweeps:
            raise EnvironmentError("Sweep {} not found in {}".format(
                name, sweeps))
        return HfssFrequencySweep(self, name)

    def add_fields_convergence_expr(self, expr, pct_delta, phase=0):
        """note: because of hfss idiocy, you must call "commit_convergence_exprs"
        after adding all exprs"""
        assert isinstance(expr, NamedCalcObject)
        self.expression_cache_items.append([
            "NAME:CacheItem", "Title:=", expr.name + "_conv", "Expression:=",
            expr.name, "Intrinsics:=", "Phase='{}deg'".format(phase),
            "IsConvergence:=", True, "UseRelativeConvergence:=", 1,
            "MaxConvergenceDelta:=", pct_delta, "MaxConvergeValue:=", "0.05",
            "ReportType:=", "Fields", ["NAME:ExpressionContext"]
        ])

    def commit_convergence_exprs(self):
        """note: this will eliminate any convergence expressions not added through this interface"""
        args = [
            "NAME:" + self.name,
            ["NAME:ExpressionCache", self.expression_cache_items]
        ]
        self._setup_module.EditSetup(self.name, args)

    def get_convergence(self, variation="", pre_fn_args=[], overwrite=True):
        '''
        Returns converge as a dataframe
            Variation should be in the form
            variation = "scale_factor='1.2001'" ...
        '''
        # TODO: (Daniel) I think this data should be store in a more comfortable datatype (dictionary maybe?)
        # Write file
        temp = tempfile.NamedTemporaryFile()
        temp.close()
        temp = temp.name + '.conv'
        self.parent._design.ExportConvergence(self.name, variation,
                                              *pre_fn_args, temp, overwrite)

        # Read File
        temp = Path(temp)
        if not temp.is_file():
            logger.error(
                f'''ERROR!  Error in trying to read temporary convergence file.
                        `get_convergence` did not seem to have the file written {str(temp)}.
                        Perhaps there was no convergence?  Check to see if there is a CONV available for this current variation. If the nominal design is not solved, it will not have a CONV., but will show up as a variation
                        Check for error messages in HFSS.
                        Retuning None''')
            return None, ''
        text = temp.read_text()

        # Parse file
        text2 = text.split(r'==================')
        if len(text) >= 3:
            df = pd.read_csv(io.StringIO(text2[3].strip()),
                             sep='|',
                             skipinitialspace=True,
                             index_col=0).drop('Unnamed: 3', axis=1)
        else:
            logger.error(f'ERROR IN reading in {temp}:\n{text}')
            df = None

        return df, text

    def get_mesh_stats(self, variation=""):
        '''  variation should be in the form
             variation = "scale_factor='1.2001'" ...
        '''
        temp = tempfile.NamedTemporaryFile()
        temp.close()
        # print(temp.name0
        # seems broken in 2016 because of extra text added to the top of the file
        self.parent._design.ExportMeshStats(self.name, variation,
                                            temp.name + '.mesh', True)
        try:
            df = pd.read_csv(temp.name + '.mesh',
                             delimiter='|',
                             skipinitialspace=True,
                             skiprows=7,
                             skipfooter=1,
                             skip_blank_lines=True,
                             engine='python')
            df = df.drop('Unnamed: 9', axis=1)
        except Exception as e:
            print("ERROR in MESH reading operation.")
            print(e)
            print(
                'ERROR!  Error in trying to read temporary MESH file ' +
                temp.name +
                '\n. Check to see if there is a mesh available for this current variation.\
                   If the nominal design is not solved, it will not have a mesh., \
                   but will show up as a variation.')
            df = None
        return df

    def get_profile(self, variation=""):
        fn = tempfile.mktemp()
        self.parent._design.ExportProfile(self.name, variation, fn, False)
        df = pd.read_csv(fn,
                         delimiter='\t',
                         skipinitialspace=True,
                         skiprows=6,
                         skipfooter=1,
                         skip_blank_lines=True,
                         engine='python')
        # just broken down by new lines
        return df

    def get_fields(self):
        return HfssFieldsCalc(self)


class HfssDMSetup(HfssSetup):
    """
    Driven modal setup
    """
    solution_freq = make_float_prop("Solution Freq")
    delta_s = make_float_prop("Delta S")
    solver_type = make_str_prop("Solver Type")

    def setup_link(self, linked_setup):
        '''
            type: linked_setup <HfssSetup>
        '''
        args = [
            "NAME:" + self.name,
            [
                "NAME:MeshLink",
                "Project:=",
                "This Project*",
                "Design:=",
                linked_setup.parent.name,
                "Soln:=",
                linked_setup.solution_name,
                self._map_variables_by_name(),
                "ForceSourceToSolve:=",
                True,
                "PathRelativeTo:=",
                "TargetProject",
            ],
        ]
        self._setup_module.EditSetup(self.name, args)

    def _map_variables_by_name(self):
        ''' does not check that variables are all present '''
        # don't care about values
        project_variables = self.parent.parent.get_variable_names()
        design_variables = self.parent.get_variable_names()

        # build array
        args = [
            "NAME:Params",
        ]
        for name in project_variables:
            args.extend([str(name) + ":=", str(name)])
        for name in design_variables:
            args.extend([str(name) + ":=", str(name)])
        return args

    def get_solutions(self):
        return HfssDMDesignSolutions(self, self.parent._solutions)

class HfssDTSetup(HfssDMSetup):

    def get_solutions(self):
        return HfssDTDesignSolutions(self, self.parent._solutions)


class HfssEMSetup(HfssSetup):
    """
    Eigenmode setup
    """
    min_freq = make_float_prop("Min Freq")
    n_modes = make_int_prop("Modes")
    delta_f = make_float_prop("Delta F")

    def get_solutions(self):
        return HfssEMDesignSolutions(self, self.parent._solutions)


class AnsysQ3DSetup(HfssSetup):
    """
    Q3D setup
    """
    prop_tab = "CG"
    max_pass = make_int_prop("Max. Number of Passes")
    min_pass = make_int_prop("Min. Number of Passes")
    pct_error = make_int_prop("Percent Error")
    frequency = make_str_prop("Adaptive Freq", 'General')  # e.g., '5GHz'
    n_modes = 0  # for compatibility with eigenmode

    def get_frequency_Hz(self):
        return int(ureg(self.frequency).to('Hz').magnitude)

    def get_solutions(self):
        return HfssQ3DDesignSolutions(self, self.parent._solutions)

    def get_convergence(self, variation=""):
        '''
        Returns df
                    # Triangle   Delta %
            Pass
            1            164       NaN
        '''
        return super().get_convergence(variation, pre_fn_args=['CG'])

    def get_matrix(
            self,
            variation='',
            pass_number=0,
            frequency=None,
            MatrixType='Maxwell',
            solution_kind='LastAdaptive',  # AdaptivePass
            ACPlusDCResistance=False,
            soln_type="C"):
        '''
        Arguments:
        -----------
            variation: an empty string returns nominal variation.
                        Otherwise need the list
            frequency: in Hz
            soln_type = "C", "AC RL" and "DC RL"
            solution_kind = 'LastAdaptive' # AdaptivePass
        Internals:
        -----------
            Uses self.solution_name  = Setup1 : LastAdaptive

        Returns:
        ---------------------
            df_cmat, user_units, (df_cond, units_cond), design_variation
        '''
        if frequency is None:
            frequency = self.get_frequency_Hz()

        temp = tempfile.NamedTemporaryFile()
        temp.close()
        path = temp.name + '.txt'
        # <FileName>, <SolnType>, <DesignVariationKey>, <Solution>, <Matrix>, <ResUnit>,
        # <IndUnit>, <CapUnit>, <CondUnit>, <Frequency>, <MatrixType>, <PassNumber>,
        # <ACPlusDCResistance>
        logger.info(f'Exporting matrix data to ({path}, {soln_type}, {variation}, '
                                             f'{self.name}:{solution_kind}, '
                                             '"Original", "ohm", "nH", "fF", '
                                             f'"mSie", {frequency}, {MatrixType}, '
                                             f'{pass_number}, {ACPlusDCResistance}')
        self.parent._design.ExportMatrixData(path, soln_type, variation,
                                             f'{self.name}:{solution_kind}',
                                             "Original", "ohm", "nH", "fF",
                                             "mSie", frequency, MatrixType,
                                             pass_number, ACPlusDCResistance)

        df_cmat, user_units, (df_cond, units_cond), design_variation = \
            self.load_q3d_matrix(path)
        return df_cmat, user_units, (df_cond, units_cond), design_variation

    @staticmethod
    def _readin_Q3D_matrix(path: str):
        """
        Read in the txt file created from q3d export
        and output the capacitance matrix

        When exporting pick "save as type: data table"

        See Zlatko

        RETURNS: Dataframe

        Example file:
        ```
        DesignVariation:$BBoxL='650um' $boxH='750um' $boxL='2mm' $QubitGap='30um' \
                        $QubitH='90um' \$QubitL='450um' Lj_1='13nH'
        Setup1:LastAdaptive
        Problem Type:C
        C Units:farad, G Units:mSie
        Reduce Matrix:Original
        Frequency: 5.5E+09 Hz

        Capacitance Matrix
            ground_plane	Q1_bus_Q0_connector_pad	Q1_bus_Q2_connector_pad	Q1_pad_bot	Q1_pad_top1	Q1_readout_connector_pad
        ground_plane	2.8829E-13	-3.254E-14	-3.1978E-14	-4.0063E-14	-4.3842E-14	-3.0053E-14
        Q1_bus_Q0_connector_pad	-3.254E-14	4.7257E-14	-2.2765E-16	-1.269E-14	-1.3351E-15	-1.451E-16
        Q1_bus_Q2_connector_pad	-3.1978E-14	-2.2765E-16	4.5327E-14	-1.218E-15	-1.1552E-14	-5.0414E-17
        Q1_pad_bot	-4.0063E-14	-1.269E-14	-1.218E-15	9.5831E-14	-3.2415E-14	-8.3665E-15
        Q1_pad_top1	-4.3842E-14	-1.3351E-15	-1.1552E-14	-3.2415E-14	9.132E-14	-1.0199E-15
        Q1_readout_connector_pad	-3.0053E-14	-1.451E-16	-5.0414E-17	-8.3665E-15	-1.0199E-15	3.9884E-14

        Conductance Matrix
            ground_plane	Q1_bus_Q0_connector_pad	Q1_bus_Q2_connector_pad	Q1_pad_bot	Q1_pad_top1	Q1_readout_connector_pad
        ground_plane	0	0	0	0	0	0
        Q1_bus_Q0_connector_pad	0	0	0	0	0	0
        Q1_bus_Q2_connector_pad	0	0	0	0	0	0
        Q1_pad_bot	0	0	0	0	0	0
        Q1_pad_top1	0	0	0	0	0	0
        Q1_readout_connector_pad	0	0	0	0	0	0
        ```
        """

        text = Path(path).read_text()

        s1 = text.split('Capacitance Matrix')
        assert len(s1) == 2, "Could not split text to `Capacitance Matrix`"

        s2 = s1[1].split('Conductance Matrix')

        df_cmat = pd.read_csv(io.StringIO(s2[0].strip()),
                              delim_whitespace=True,
                              skipinitialspace=True,
                              index_col=0)
        units = re.findall(r'C Units:(.*?),', text)[0]

        if len(s2) > 1:
            df_cond = pd.read_csv(io.StringIO(s2[1].strip()),
                                  delim_whitespace=True,
                                  skipinitialspace=True,
                                  index_col=0)
            units_cond = re.findall(r'G Units:(.*?)\n', text)[0]
        else:
            df_cond = None

        var = re.findall(r'DesignVariation:(.*?)\n',
                         text)  # this changed circa v2020
        if len(var) < 1:  # didnt find
            var = re.findall(r'Design Variation:(.*?)\n', text)
            if len(var) < 1:  # didnt find
                # May not be present if there are no design variations to begin
                # with and no variables in the design.
                pass  #logger.error(f'Failed to parse Q3D matrix Design Variation:\nFile:{path}\nText:{text}')

                var = ['']
        design_variation = var[0]

        return df_cmat, units, design_variation, df_cond, units_cond

    @staticmethod
    def load_q3d_matrix(path, user_units='fF'):
        """Load Q3D capacitance file exported as Maxwell matrix.
        Exports also conductance conductance.
        Units are read in automatically and converted to user units.

        Arguments:
            path {[str or Path]} -- [path to file text with matrix]

        Returns:
            df_cmat, user_units, (df_cond, units_cond), design_variation

            dataframes: df_cmat, df_cond
        """
        df_cmat, Cunits, design_variation, df_cond, units_cond = AnsysQ3DSetup._readin_Q3D_matrix(
            path)

        # Unit convert
        q = ureg.parse_expression(Cunits).to(user_units)
        df_cmat = df_cmat * q.magnitude  # scale to user units

        #print("Imported capacitance matrix with UNITS: [%s] now converted to USER UNITS:[%s] from file:\n\t%s"%(Cunits, user_units, path))

        return df_cmat, user_units, (df_cond, units_cond), design_variation


class HfssDesignSolutions(COMWrapper):
    def __init__(self, setup, solutions):
        '''
        :type setup: HfssSetup
        '''
        super(HfssDesignSolutions, self).__init__()
        self.parent = setup
        self._solutions = solutions
        self._ansys_version = self.parent._ansys_version

    def get_valid_solution_list(self):
        '''
         Gets all available solution names that exist in a design.
         Return example:
            ('Setup1 : AdaptivePass', 'Setup1 : LastAdaptive')
        '''
        return self._solutions.GetValidISolutionList()

    def list_variations(self, setup_name: str = None):
        """
        Get a list of solved variations.

        Args:
            setup_name(str) : Example name ("Setup1 : LastAdaptive") Defaults to None.

        Returns:
             An array of strings corresponding to solved variations.

             .. code-block:: python

                ("Cj='2fF' Lj='12nH'",
                "Cj='2fF' Lj='12.5nH'",
                "Cj='2fF' Lj='13nH'",
                "Cj='2fF' Lj='13.5nH'",
                "Cj='2fF' Lj='14nH'")
        """
        if setup_name is None:
            setup_name = str(self.parent.solution_name)
        return self._solutions.ListVariations(setup_name)


class HfssEMDesignSolutions(HfssDesignSolutions):
    def eigenmodes(self, lv=""):
        '''
        Returns the eigenmode data of freq and kappa/2p
        '''
        fn = tempfile.mktemp()
        #print(self.parent.solution_name, lv, fn)
        self._solutions.ExportEigenmodes(self.parent.solution_name, lv, fn)
        data = np.genfromtxt(fn, dtype='str')
        # Update to Py 3:
        # np.loadtxt and np.genfromtxt operate in byte mode, which is the default string type in Python 2.
        # But Python 3 uses unicode, and marks bytestrings with this b.
        # getting around the very annoying fact that
        if np.size(np.shape(data)) == 1:
            # in Python a 1D array does not have shape (N,1)
            data = np.array([data])
        else:  # but rather (N,) ....
            pass
        if np.size(data[0, :]) == 6:  # checking if values for Q were saved
            # eigvalue=(omega-i*kappa/2)/2pi
            kappa_over_2pis = [2 * float(ii) for ii in data[:, 3]]
            # so kappa/2pi = 2*Im(eigvalue)
        else:
            kappa_over_2pis = None

        # print(data[:,1])
        freqs = [float(ii) for ii in data[:, 1]]
        return freqs, kappa_over_2pis

    """
    Export eigenmodes vs pass number
    Did not figure out how to set pass number in a hurry.


    import tempfile
    self = epr_hfss.solutions

    '''
    HFSS: Exports a tab delimited table of Eigenmodes in HFSS. Not in HFSS-IE.
    <setupName> <solutionName> <DesignVariationKey>
    <filename>
    Return Value:    None

    Parameters:
        <SolutionName>
            Type: <string>
            Name of the solutions within the solution setup.
        <DesignVariationKey>
            Type: <string>
            Design variation string.
    '''
    setup = self.parent
    fn = tempfile.mktemp()
    variation_list=''
    soln_name = f'{setup.name} : AdaptivePas'
    available_solns = self._solutions.GetValidISolutionList()
    if not(soln_name in available_solns):
        logger.error(f'ERROR Tried to export freq vs pass number, but solution  `{soln_name}` was not in available `{available_solns}`. Returning []')
        #return []
    self._solutions.ExportEigenmodes(soln_name, ['Pass:=5'], fn) # ['Pass:=5'] fails  can do with ''
    """

    def set_mode(self, n, phase=0, FieldType='EigenStoredEnergy'):
        '''
        Indicates which source excitations should be used for fields post processing.
        HFSS>Fields>Edit Sources

        Mode count starts at 1

        Amplitude is set to 1

        No error is thrown if a number exceeding number of modes is set

            FieldType -- EigenStoredEnergy or EigenPeakElecticField
        '''
        n_modes = int(self.parent.n_modes)

        if n < 1:
            err = f'ERROR: You tried to set a mode < 1. {n}/{n_modes}'
            logger.error(err)
            raise Exception(err)

        if n > n_modes:
            err = f'ERROR: You tried to set a mode > number of modes {n}/{n_modes}'
            logger.error(err)
            raise Exception(err)

        if self._ansys_version >= '2019':
            # THIS WORKS FOR v2019R2
            self._solutions.EditSources(
                [["FieldType:=", "EigenPeakElectricField"],
                 [
                     "Name:=", "Modes", "Magnitudes:=",
                     ["1" if i + 1 == n else "0" for i in range(n_modes)],
                     "Phases:=",
                     [
                         str(phase) if i + 1 == n else "0"
                         for i in range(n_modes)
                     ]
                 ]])
        else:
            # The syntax has changed for AEDT 18.2.
            # see https://ansyshelp.ansys.com/account/secured?returnurl=/Views/Secured/Electronics/v195//Subsystems/HFSS/Subsystems/HFSS%20Scripting/HFSS%20Scripting.htm

            self._solutions.EditSources(
                "EigenStoredEnergy", ["NAME:SourceNames", "EigenMode"],
                ["NAME:Modes", n_modes], ["NAME:Magnitudes"] +
                [1 if i + 1 == n else 0
                 for i in range(n_modes)], ["NAME:Phases"] +
                [phase if i + 1 == n else 0 for i in range(n_modes)],
                ["NAME:Terminated"], ["NAME:Impedances"])

    def has_fields(self, variation_string=None):
        '''
        Determine if fields exist for a particular solution.

        variation_string : str | None
            This must the string that describes the variation in hFSS, not 0 or 1, but
            the string of variables, such as
                "Cj='2fF' Lj='12.75nH'"
            If None, gets the nominal variation
        '''
        if variation_string is None:
            variation_string = self.parent.parent.get_nominal_variation()

        return bool(
            self._solutions.HasFields(self.parent.solution_name,
                                      variation_string))

    def create_report(self,
                      plot_name,
                      xcomp,
                      ycomp,
                      params,
                      pass_name='LastAdaptive'):
        '''
        pass_name: AdaptivePass, LastAdaptive

        Example
        ------------------------------------------------------
        Example plot for a single variation all pass converge of mode freq
        .. code-block python
            ycomp = [f"re(Mode({i}))" for i in range(1,1+epr_hfss.n_modes)]
            params = ["Pass:=", ["All"]]+variation
            setup.create_report("Freq. vs. pass", "Pass", ycomp, params, pass_name='AdaptivePass')
        '''
        assert isinstance(ycomp, list)
        assert isinstance(params, list)

        setup = self.parent
        reporter = setup._reporter
        return reporter.CreateReport(
            plot_name, "Eigenmode Parameters", "Rectangular Plot",
            f"{setup.name} : {pass_name}", [], params,
            ["X Component:=", xcomp, "Y Component:=", ycomp], [])


class HfssDMDesignSolutions(HfssDesignSolutions):
    pass

class HfssDTDesignSolutions(HfssDesignSolutions):
    pass

class HfssQ3DDesignSolutions(HfssDesignSolutions):
    pass


class HfssFrequencySweep(COMWrapper):
    prop_tab = "HfssTab"
    start_freq = make_float_prop("Start")
    stop_freq = make_float_prop("Stop")
    step_size = make_float_prop("Step Size")
    count = make_float_prop("Count")
    sweep_type = make_str_prop("Type")

    def __init__(self, setup, name):
        """
        :type setup: HfssSetup
        :type name: str
        """
        super(HfssFrequencySweep, self).__init__()
        self.parent = setup
        self.name = name
        self.solution_name = self.parent.name + " : " + name
        self.prop_holder = self.parent.prop_holder
        self.prop_server = self.parent.prop_server + ":" + name
        self._ansys_version = self.parent._ansys_version

    def analyze_sweep(self):
        self.parent.analyze(self.solution_name)

    def get_network_data(self, formats):
        if isinstance(formats, str):
            formats = formats.split(",")
        formats = [f.upper() for f in formats]

        fmts_lists = {'S': [], 'Y': [], 'Z': []}

        for f in formats:
            fmts_lists[f[0]].append((int(f[1]), int(f[2])))
        ret = [None] * len(formats)
        freq = None

        for data_type, list in fmts_lists.items():
            if list:
                fn = tempfile.mktemp()
                self.parent._solutions.ExportNetworkData(
                    [], self.parent.name + " : " + self.name, 2, fn, ["all"],
                    False, 0, data_type, -1, 1, 15)
                with open(fn) as f:
                    f.readline()
                    colnames = f.readline().split()
                array = np.loadtxt(fn, skiprows=2)
                # WARNING for python 3 probably need to use genfromtxt
                if freq is None:
                    freq = array[:, 0]
                # TODO: If Ansys version is 2019, use 'Real' and 'Imag'
                # in place of 'Re' and 'Im
                for i, j in list:
                    real_idx = colnames.index("%s[%d,%d]_Re" %
                                            (data_type, i, j))
                    imag_idx = colnames.index("%s[%d,%d]_Im" %
                                            (data_type, i, j))
                    c_arr = array[:, real_idx] + 1j * array[:, imag_idx]
                    ret[formats.index("%s%d%d" % (data_type, i, j))] = c_arr

        return freq, ret

    def create_report(self, name, expr):
        existing = self.parent._reporter.GetAllReportNames()
        name = increment_name(name, existing)
        var_names = self.parent.parent.get_variable_names()
        var_args = sum([["%s:=" % v_name, ["Nominal"]]
                        for v_name in var_names], [])
        self.parent._reporter.CreateReport(
            name, "Modal Solution Data", "Rectangular Plot",
            self.solution_name, ["Domain:=", "Sweep"],
            ["Freq:=", ["All"]] + var_args,
            ["X Component:=", "Freq", "Y Component:=", [expr]], [])
        return HfssReport(self.parent.parent, name)

    def get_report_arrays(self, expr):
        r = self.create_report("Temp", expr)
        return r.get_arrays()


class HfssReport(COMWrapper):
    def __init__(self, design, name):
        """
        :type design: HfssDesign
        :type name: str
        """
        super(HfssReport, self).__init__()
        self.parent_design = design
        self.name = name

    def export_to_file(self, filename):
        filepath = os.path.abspath(filename)
        self.parent_design._reporter.ExportToFile(self.name, filepath)

    def get_arrays(self):
        fn = tempfile.mktemp(suffix=".csv")
        self.export_to_file(fn)
        return np.loadtxt(fn, skiprows=1, delimiter=',').transpose()
        # warning for python 3 probably need to use genfromtxt


class Optimetrics(COMWrapper):
    """
    Optimetrics script commands executed by the "Optimetrics" module.

    Example use:
    .. code-block python
            opti = Optimetrics(pinfo.design)
            names = opti.get_setup_names()
            print('Names of optimetrics: ', names)
            opti.solve_setup(names[0])

    Note that running optimetrics requires the license for Optimetrics by Ansys.
    """
    def __init__(self, design):
        super(Optimetrics, self).__init__()

        self.design = design  # parent
        self._optimetrics = self.design._optimetrics  # <COMObject GetModule>
        self.setup_names = None

    def get_setup_names(self):
        """
        Return list of Optimetrics setup names
        """
        self.setup_names = list(self._optimetrics.GetSetupNames())
        return self.setup_names.copy()

    def solve_setup(self, setup_name: str):
        """
        Solves the specified Optimetrics setup.
        Corresponds to:  Right-click the setup in the project tree, and then click
        Analyze on the shortcut menu.

        setup_name (str) : name of setup, should be in get_setup_names

        Blocks execution until ready to use.

        Note that this requires the license for Optimetrics by Ansys.
        """
        return self._optimetrics.SolveSetup(setup_name)

    def create_setup(self,
                     variable,
                     swp_params,
                     name="ParametricSetup1",
                     swp_type='linear_step',
                     setup_name=None,
                     save_fields=True,
                     copy_mesh=True,
                     solve_with_copied_mesh_only=True,
                     setup_type='parametric'):
        """
        Inserts a new parametric setup of one variable. Either with sweep 
        definition or from file.

        Corresponds to ui access:
        Right-click the Optimetrics folder in the project tree, and then click 
        Add> Parametric on the shortcut menu.

        Ansys provides six sweep definitions types specified using the swp_type
        variable.

        Sweep type definitions:
        - 'single_value'          	
            Specify a single value for the sweep definition.
        - 'linear_step'
            Specify a linear range of values with a constant step size.
        - 'linear_count'
            Specify a linear range of values and the number, or count of points
            within this range.
        - 'decade_count'
            Specify a logarithmic (base 10) series of values, and the number of
            values to calculate in each decade.
        - 'octave_count'
            Specify a logarithmic (base 2) series of values, and the number of 
            values to calculate in each octave.
        - 'exponential_count'
            Specify an exponential (base e) series of values, and the number of 
            values to calculate.

        For swp_type='single_value' swp_params is the single value.

        For  swp_type='linear_step' swp_params is start, stop, step:
            swp_params = ("12.8nH", "13.6nH", "0.2nH")
        
        All other types swp_params is start, stop, count:
            swp_params = ("12.8nH", "13.6nH", 4)
            The definition of count varies amongst the available types. 

        For Decade count and Octave count, the Count value specifies the number
        of points to calculate in every decade or octave. For Exponential count,
        the Count value is the total number of points. The total number of 
        points includes the start and stop values.

        For parametric from file, setup_type='parametric_file', pass in a file
        name and path to swp_params like "C:\\test.csv" or "C:\\test.txt" for 
        example.

        Example csv formatting:
        *,Lj_qubit
        1,12.2nH
        2,9.7nH
        3,10.2nH

        See Ansys documentation for additional formatting instructions. 
        """
        setup_name = setup_name or self.design.get_setup_names()[0]
        print(
            f"Inserting optimetrics setup `{name}` for simulation setup: `{setup_name}`"
        )

        if setup_type == 'parametric':
            valid_swp_types = ['single_value', 'linear_step', 'linear_count', 
            'decade_count', 'octave_count', 'exponential_count']

            if swp_type not in valid_swp_types:
                raise NotImplementedError()
            else:
                if swp_type == 'single_value':
                    # Single takes string of single variable no swp_type_name
                    swp_str = f"{swp_params}"

                else:
                    # correct number of inputs
                    assert len(swp_params) == 3, "Incorrect number of sweep parameters."

                    # Not checking for compatible unit types
                    if swp_type == 'linear_step':
                        swp_type_name = "LIN"
                    else:
                        # counts needs to be an integer number
                        assert isinstance(swp_params[2], int), "Count must be integer."

                        if swp_type == 'linear_count':
                            swp_type_name = "LINC"
                        elif swp_type == 'decade_count':
                            swp_type_name = "DEC"
                        elif swp_type == 'octave_count':
                            swp_type_name = "OCT"
                        elif swp_type == 'exponential_count':
                            swp_type_name = "ESTP"

                    # prepare the string to pass to Ansys
                    swp_str = f"{swp_type_name} {swp_params[0]} {swp_params[1]} {swp_params[2]}"

            # talk with Ansys
            self._optimetrics.InsertSetup("OptiParametric", [
                f"NAME:{name}", "IsEnabled:=", True,
                [
                    "NAME:ProdOptiSetupDataV2",
                    "SaveFields:=",
                    save_fields,
                    "CopyMesh:=",
                    copy_mesh,
                    "SolveWithCopiedMeshOnly:=",
                    solve_with_copied_mesh_only,
                ], ["NAME:StartingPoint"], "Sim. Setups:=", [setup_name],
                [
                    "NAME:Sweeps",
                    [
                        "NAME:SweepDefinition", "Variable:=", variable, "Data:=",
                        swp_str, "OffsetF1:=", False, "Synchronize:=", 0
                    ]
                ], ["NAME:Sweep Operations"], ["NAME:Goals"]
            ])
        elif setup_type == 'parametric_file':
            # Uses the file name as the swp_params 
            filename = swp_params

            self._optimetrics.ImportSetup("OptiParametric",
                [
                f"NAME:{name}",
                filename,
                ])
            self._optimetrics.EditSetup(f"{name}",
                [
                    f"NAME:{name}",
            		[
            			"NAME:ProdOptiSetupDataV2",
            			"SaveFields:="		, save_fields,
            			"CopyMesh:="		, copy_mesh,
            			"SolveWithCopiedMeshOnly:=", solve_with_copied_mesh_only,
            		],
            ])
        else:
            raise NotImplementedError()


class HfssModeler(COMWrapper):
    def __init__(self, design, modeler, boundaries, mesh):
        """
        :type design: HfssDesign
        """
        super(HfssModeler, self).__init__()
        self.parent = design
        self._modeler = modeler
        self._boundaries = boundaries
        self._mesh = mesh  # Mesh module

    def set_units(self, units, rescale=True):
        self._modeler.SetModelUnits(
            ["NAME:Units Parameter", "Units:=", units, "Rescale:=", rescale])

    def get_units(self):
        """Get the model units.
            Return Value:    A string contains current model units. """
        return str(self._modeler.GetModelUnits())

    def get_all_properties(self, obj_name, PropTab='Geometry3DAttributeTab'):
        '''
            Get all properties for modeler PropTab, PropServer
        '''
        PropServer = obj_name
        properties = {}
        for key in self._modeler.GetProperties(PropTab, PropServer):
            properties[key] = self._modeler.GetPropertyValue(
                PropTab, PropServer, key)
        return properties

    def _attributes_array(
            self,
            name=None,
            nonmodel=False,
            wireframe=False,
            color=None,
            transparency=0.9,
            material=None,  # str
            solve_inside=None,  # bool
            coordinate_system="Global"):
        arr = ["NAME:Attributes", "PartCoordinateSystem:=", coordinate_system]
        if name is not None:
            arr.extend(["Name:=", name])

        if nonmodel or wireframe:
            flags = 'NonModel' if nonmodel else ''  # can be done smarter
            if wireframe:
                flags += '#' if len(flags) > 0 else ''
                flags += 'Wireframe'
            arr.extend(["Flags:=", flags])

        if color is not None:
            arr.extend(["Color:=", "(%d %d %d)" % color])
        if transparency is not None:
            arr.extend(["Transparency:=", transparency])
        if material is not None:
            arr.extend(["MaterialName:=", material])
        if solve_inside is not None:
            arr.extend(["SolveInside:=", solve_inside])

        return arr

    def _selections_array(self, *names):
        return ["NAME:Selections", "Selections:=", ",".join(names)]

    def mesh_length(self,
                    name_mesh,
                    objects: list,
                    MaxLength='0.1mm',
                    **kwargs):
        '''
        "RefineInside:="	, False,
        "Enabled:="		, True,
        "RestrictElem:="	, False,
        "NumMaxElem:="		, "1000",
        "RestrictLength:="	, True,
        "MaxLength:="		, "0.1mm"

        Example use:
        modeler.assign_mesh_length('mesh2', ["Q1_mesh"], MaxLength=0.1)
        '''
        assert isinstance(objects, list)

        arr = [
            f"NAME:{name_mesh}", "Objects:=", objects, 'MaxLength:=', MaxLength
        ]
        ops = [
            'RefineInside', 'Enabled', 'RestrictElem', 'NumMaxElem',
            'RestrictLength'
        ]
        for key, val in kwargs.items():
            if key in ops:
                arr += [key + ':=', str(val)]
            else:
                logger.error('KEY `{key}` NOT IN ops!')

        self._mesh.AssignLengthOp(arr)

    def mesh_reassign(self, name_mesh, objects: list):
        assert isinstance(objects, list)
        self._mesh.ReassignOp(name_mesh, ["Objects:=", objects])

    def mesh_get_names(self, kind="Length Based"):
        ''' "Length Based", "Skin Depth Based", ...'''
        return list(self._mesh.GetOperationNames(kind))

    def mesh_get_all_props(self, mesh_name):
        # TODO: make mesh tis own  class with properties
        prop_tab = 'MeshSetupTab'
        prop_server = f'MeshSetup:{mesh_name}'
        prop_names = self.parent._design.GetProperties('MeshSetupTab',
                                                       prop_server)
        dic = {}
        for name in prop_names:
            dic[name] = self._modeler.GetPropertyValue(prop_tab, prop_server,
                                                       name)
        return dic

    def draw_box_corner(self, pos, size, **kwargs):
        name = self._modeler.CreateBox([
            "NAME:BoxParameters", "XPosition:=",
            str(pos[0]), "YPosition:=",
            str(pos[1]), "ZPosition:=",
            str(pos[2]), "XSize:=",
            str(size[0]), "YSize:=",
            str(size[1]), "ZSize:=",
            str(size[2])
        ], self._attributes_array(**kwargs))
        return Box(name, self, pos, size)

    def draw_box_center(self, pos, size, **kwargs):
        """
        Creates a 3-D box centered at pos [x0, y0, z0], with width
        size [xwidth, ywidth, zwidth] along each respective direction.

        Args:
            pos (list): Coordinates of center of box, [x0, y0, z0]
            size (list): Width of box along each direction, [xwidth, ywidth, zwidth]
        """
        corner_pos = [var(p) - var(s) / 2 for p, s in zip(pos, size)]
        return self.draw_box_corner(corner_pos, size, **kwargs)

    def draw_polyline(self, points, closed=True, **kwargs):
        """
        Draws a closed or open polyline.
        If closed = True, then will make into a sheet.
        points : need to be in the correct units

        For optional arguments, see _attributes_array; these include:
        ```
            nonmodel=False,
            wireframe=False,
            color=None,
            transparency=0.9,
            material=None,  # str
            solve_inside=None,  # bool
            coordinate_system="Global"
        ```
        """
        pointsStr = ["NAME:PolylinePoints"]
        indexsStr = ["NAME:PolylineSegments"]
        for ii, point in enumerate(points):
            pointsStr.append([
                "NAME:PLPoint", "X:=",
                str(point[0]), "Y:=",
                str(point[1]), "Z:=",
                str(point[2])
            ])
            indexsStr.append([
                "NAME:PLSegment", "SegmentType:=", "Line", "StartIndex:=", ii,
                "NoOfPoints:=", 2
            ])
        if closed:
            pointsStr.append([
                "NAME:PLPoint", "X:=",
                str(points[0][0]), "Y:=",
                str(points[0][1]), "Z:=",
                str(points[0][2])
            ])
            params_closed = [
                "IsPolylineCovered:=", True, "IsPolylineClosed:=", True
            ]
        else:
            indexsStr = indexsStr[:-1]
            params_closed = [
                "IsPolylineCovered:=", True, "IsPolylineClosed:=", False
            ]

        name = self._modeler.CreatePolyline(
            ["NAME:PolylineParameters", *params_closed, pointsStr, indexsStr],
            self._attributes_array(**kwargs))

        if closed:
            return Polyline(name, self, points)
        else:
            return OpenPolyline(name, self, points)

    def draw_rect_corner(self, pos, x_size=0, y_size=0, z_size=0, **kwargs):
        size = [x_size, y_size, z_size]
        assert 0 in size
        axis = "XYZ"[size.index(0)]
        w_idx, h_idx = {'X': (1, 2), 'Y': (2, 0), 'Z': (0, 1)}[axis]

        name = self._modeler.CreateRectangle([
            "NAME:RectangleParameters", "XStart:=",
            str(pos[0]), "YStart:=",
            str(pos[1]), "ZStart:=",
            str(pos[2]), "Width:=",
            str(size[w_idx]), "Height:=",
            str(size[h_idx]), "WhichAxis:=", axis
        ], self._attributes_array(**kwargs))
        return Rect(name, self, pos, size)

    def draw_rect_center(self, pos, x_size=0, y_size=0, z_size=0, **kwargs):
        """
        Creates a rectangle centered at pos [x0, y0, z0].
        It is assumed that the rectangle lies parallel to the xy, yz, or xz plane.
        User inputs 2 of 3 of the following: x_size, y_size, and z_size
        depending on how the rectangle is oriented.

        Args:
            pos (list): Coordinates of rectangle center, [x0, y0, z0]
            x_size (int, optional): Width along the x direction. Defaults to 0.
            y_size (int, optional):  Width along the y direction. Defaults to 0.
            z_size (int, optional):  Width along the z direction]. Defaults to 0.
        """
        corner_pos = [
            var(p) - var(s) / 2. for p, s in zip(pos, [x_size, y_size, z_size])
        ]
        return self.draw_rect_corner(corner_pos, x_size, y_size, z_size,
                                     **kwargs)

    def draw_cylinder(self, pos, radius, height, axis, **kwargs):
        assert axis in "XYZ"
        return self._modeler.CreateCylinder([
            "NAME:CylinderParameters", "XCenter:=", pos[0], "YCenter:=",
            pos[1], "ZCenter:=", pos[2], "Radius:=", radius, "Height:=",
            height, "WhichAxis:=", axis, "NumSides:=", 0
        ], self._attributes_array(**kwargs))

    def draw_cylinder_center(self, pos, radius, height, axis, **kwargs):
        axis_idx = ["X", "Y", "Z"].index(axis)
        edge_pos = copy(pos)
        edge_pos[axis_idx] = var(pos[axis_idx]) - var(height) / 2
        return self.draw_cylinder(edge_pos, radius, height, axis, **kwargs)

    def draw_wirebond(self,
                      pos,
                      ori,
                      width,
                      height='0.1mm',
                      z=0,
                      wire_diameter="0.02mm",
                      NumSides=6,
                      **kwargs):
        '''
            Args:
                pos: 2D position vector  (specify center point)
                ori: should be normed
                z: z position

            # TODO create Wirebond class
            position is the origin of one point
            ori is the orientation vector, which gets normalized
        '''
        p = np.array(pos)
        o = np.array(ori)
        pad1 = p - o * width / 2.
        name = self._modeler.CreateBondwire([
            "NAME:BondwireParameters", "WireType:=", "Low", "WireDiameter:=",
            wire_diameter, "NumSides:=", NumSides, "XPadPos:=", pad1[0],
            "YPadPos:=", pad1[1], "ZPadPos:=", z, "XDir:=", ori[0], "YDir:=",
            ori[1], "ZDir:=", 0, "Distance:=", width, "h1:=", height, "h2:=",
            "0mm", "alpha:=", "80deg", "beta:=", "80deg", "WhichAxis:=", "Z"
        ], self._attributes_array(**kwargs))

        return name

    def draw_region(self,
                    Padding,
                    PaddingType="Percentage Offset",
                    name='Region',
                    material="\"vacuum\""):
        """
            PaddingType : 'Absolute Offset', "Percentage Offset"
        """
        # TODO: Add option to modify these
        RegionAttributes = [
            "NAME:Attributes", "Name:=", name, "Flags:=", "Wireframe#",
            "Color:=", "(255 0 0)", "Transparency:=", 1,
            "PartCoordinateSystem:=", "Global", "UDMId:=", "",
            "IsAlwaysHiden:=", False, "MaterialValue:=", material,
            "SolveInside:=", True
        ]

        self._modeler.CreateRegion([
            "NAME:RegionParameters", "+XPaddingType:=", PaddingType,
            "+XPadding:=", Padding[0][0], "-XPaddingType:=", PaddingType,
            "-XPadding:=", Padding[0][1], "+YPaddingType:=", PaddingType,
            "+YPadding:=", Padding[1][0], "-YPaddingType:=", PaddingType,
            "-YPadding:=", Padding[1][1], "+ZPaddingType:=", PaddingType,
            "+ZPadding:=", Padding[2][0], "-ZPaddingType:=", PaddingType,
            "-ZPadding:=", Padding[2][1]
        ], RegionAttributes)

    def unite(self, names, keep_originals=False):
        self._modeler.Unite(
            self._selections_array(*names),
            ["NAME:UniteParameters", "KeepOriginals:=", keep_originals])
        return names[0]

    def intersect(self, names, keep_originals=False):
        self._modeler.Intersect(
            self._selections_array(*names),
            ["NAME:IntersectParameters", "KeepOriginals:=", keep_originals])
        return names[0]

    def translate(self, name, vector):
        self._modeler.Move(self._selections_array(name), [
            "NAME:TranslateParameters", "TranslateVectorX:=", vector[0],
            "TranslateVectorY:=", vector[1], "TranslateVectorZ:=", vector[2]
        ])

    def get_boundary_assignment(self, boundary_name: str):
        # Gets a list of face IDs associated with the given boundary or excitation assignment.
        objects = self._boundaries.GetBoundaryAssignment(boundary_name)
        # Gets an object name corresponding to the input face id. Returns the name of the corresponding object name.
        objects = [self._modeler.GetObjectNameByFaceID(k) for k in objects]
        return objects

    def append_PerfE_assignment(self, boundary_name: str, object_names: list):
        '''
            This will create a new boundary if need, and will
            otherwise append given names to an existing boundary
        '''
        # enforce
        boundary_name = str(boundary_name)
        if isinstance(object_names, str):
            object_names = [object_names]
        object_names = list(object_names)  # enforce list

        # do actual work
        if boundary_name not in self._boundaries.GetBoundaries(
        ):  # GetBoundariesOfType("Perfect E")
            # need to make a new boundary
            self.assign_perfect_E(object_names, name=boundary_name)
        else:
            # need to append
            objects = list(self.get_boundary_assignment(boundary_name))
            self._boundaries.ReassignBoundary([
                "NAME:" + boundary_name, "Objects:=",
                list(set(objects + object_names))
            ])

    def append_mesh(self, mesh_name: str, object_names: list, old_objs: list,
                    **kwargs):
        '''
        This will create a new boundary if need, and will
        otherwise append given names to an existing boundary
        old_obj = circ._mesh_assign
        '''
        mesh_name = str(mesh_name)
        if isinstance(object_names, str):
            object_names = [object_names]
        object_names = list(object_names)  # enforce list

        if mesh_name not in self.mesh_get_names(
        ):  # need to make a new boundary
            objs = object_names
            self.mesh_length(mesh_name, object_names, **kwargs)
        else:  # need to append
            objs = list(set(old_objs + object_names))
            self.mesh_reassign(mesh_name, objs)

        return objs

    def assign_perfect_E(self, obj: List[str], name: str = 'PerfE'):
        '''
        Assign a boundary condition to a list of objects.

        Arg:
            objs (List[str]): Takes a name of an object or a list of object names.
            name(str): If `name` is not specified `PerfE` is appended to object name for the name.
        '''
        if not isinstance(obj, list):
            obj = [obj]
            if name == 'PerfE':
                name = str(obj) + '_' + name
        name = increment_name(name, self._boundaries.GetBoundaries())
        self._boundaries.AssignPerfectE(
            ["NAME:" + name, "Objects:=", obj, "InfGroundPlane:=", False])

    def _make_lumped_rlc(self, r, l, c, start, end, obj_arr, name="LumpRLC"):
        name = increment_name(name, self._boundaries.GetBoundaries())
        params = ["NAME:" + name]
        params += obj_arr
        params.append([
            "NAME:CurrentLine",
            # for some reason here it seems to switch to use the model units, rather than meters
            "Start:=",
            fix_units(start, unit_assumed=LENGTH_UNIT),
            "End:=",
            fix_units(end, unit_assumed=LENGTH_UNIT)
        ])
        params += [
            "UseResist:=", r != 0, "Resistance:=", r, "UseInduct:=", l != 0,
            "Inductance:=", l, "UseCap:=", c != 0, "Capacitance:=", c
        ]
        self._boundaries.AssignLumpedRLC(params)

    def _make_lumped_port(self,
                          start,
                          end,
                          obj_arr,
                          z0="50ohm",
                          name="LumpPort"):
        start = fix_units(start, unit_assumed=LENGTH_UNIT)
        end = fix_units(end, unit_assumed=LENGTH_UNIT)

        name = increment_name(name, self._boundaries.GetBoundaries())
        params = ["NAME:" + name]
        params += obj_arr
        params += [
            "RenormalizeAllTerminals:=", True, "DoDeembed:=", False,
            [
                "NAME:Modes",
                [
                    "NAME:Mode1", "ModeNum:=", 1, "UseIntLine:=", True,
                    ["NAME:IntLine", "Start:=", start, "End:=",
                     end], "CharImp:=", "Zpi", "AlignmentGroup:=", 0,
                    "RenormImp:=", "50ohm"
                ]
            ], "ShowReporterFilter:=", False, "ReporterFilter:=", [True],
            "FullResistance:=", z0, "FullReactance:=", "0ohm"
        ]

        self._boundaries.AssignLumpedPort(params)

    def get_face_ids(self, obj):
        return self._modeler.GetFaceIDs(obj)

    def get_object_name_by_face_id(self, ID: str):
        ''' Gets an object name corresponding to the input face id. '''
        return self._modeler.GetObjectNameByFaceID(ID)

    def get_vertex_ids(self, obj):
        """
            Get the vertex IDs of given an object name
            oVertexIDs = oEditor.GetVertexIDsFromObject(“Box1”)
        """
        return self._modeler.GetVertexIDsFromObject(obj)

    def eval_expr(self, expr, units="mm"):
        if not isinstance(expr, str):
            return expr
        return self.parent.eval_expr(expr, units)

    def get_objects_in_group(self, group):
        """
        Use:              Returns the objects for the specified group.
        Return Value:    The objects in the group.
        Parameters:      <groupName>  Type: <string>
        One of  <materialName>, <assignmentName>, "Non Model",
                "Solids", "Unclassi­fied", "Sheets", "Lines"
        """
        if self._modeler:
            return list(self._modeler.GetObjectsInGroup(group))
        else:
            return list()

    def set_working_coordinate_system(self, cs_name="Global"):
        """
        Use:                   Sets the working coordinate system.
        Command:         Modeler>Coordinate System>Set Working CS
        """
        self._modeler.SetWCS([
            "NAME:SetWCS Parameter",
            "Working Coordinate System:=",
            cs_name,
            "RegionDepCSOk:=",
            False  # this one is prob not needed, but comes with the record tool
        ])

    def create_relative_coorinate_system_both(self,
                                              cs_name,
                                              origin=["0um", "0um", "0um"],
                                              XAxisVec=["1um", "0um", "0um"],
                                              YAxisVec=["0um", "1um", "0um"]):
        """
        Use:     Creates a relative coordinate system. Only the    Name attribute of the <AttributesArray> parameter is supported.
        Command: Modeler>Coordinate System>Create>Relative CS->Offset
        Modeler>Coordinate System>Create>Relative CS->Rotated
        Modeler>Coordinate System>Create>Relative CS->Both

        Current coordinate system is set right after this.

        cs_name : name of coord. sys
            If the name already exists, then a new coordinate system with _1 is created.

        origin, XAxisVec, YAxisVec: 3-vectors
            You can also pass in params such as origin = [0,1,0] rather than ["0um","1um","0um"], but these will be interpreted in default units, so it is safer to be explicit. Explicit over implicit.
        """
        self._modeler.CreateRelativeCS([
            "NAME:RelativeCSParameters", "Mode:=", "Axis/Position",
            "OriginX:=", origin[0], "OriginY:=", origin[1], "OriginZ:=",
            origin[2], "XAxisXvec:=", XAxisVec[0], "XAxisYvec:=", XAxisVec[1],
            "XAxisZvec:=", XAxisVec[2], "YAxisXvec:=", YAxisVec[0],
            "YAxisYvec:=", YAxisVec[1], "YAxisZvec:=", YAxisVec[1]
        ], ["NAME:Attributes", "Name:=", cs_name])

    def subtract(self, blank_name, tool_names, keep_originals=False):
        selection_array = [
            "NAME:Selections", "Blank Parts:=", blank_name, "Tool Parts:=",
            ",".join(tool_names)
        ]
        self._modeler.Subtract(
            selection_array,
            ["NAME:UniteParameters", "KeepOriginals:=", keep_originals])
        return blank_name

    def _fillet(self, radius, vertex_index, obj):
        vertices = self._modeler.GetVertexIDsFromObject(obj)
        if isinstance(vertex_index, list):
            to_fillet = [int(vertices[v]) for v in vertex_index]
        else:
            to_fillet = [int(vertices[vertex_index])]


#        print(vertices)
#        print(radius)
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj], [
            "NAME:Parameters",
            [
                "NAME:FilletParameters", "Edges:=", [], "Vertices:=",
                to_fillet, "Radius:=", radius, "Setback:=", "0mm"
            ]
        ])

    def _fillet_edges(self, radius, edge_index, obj):
        edges = self._modeler.GetEdgeIDsFromObject(obj)
        if isinstance(edge_index, list):
            to_fillet = [int(edges[e]) for e in edge_index]
        else:
            to_fillet = [int(edges[edge_index])]

        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj], [
            "NAME:Parameters",
            [
                "NAME:FilletParameters", "Edges:=", to_fillet, "Vertices:=",
                [], "Radius:=", radius, "Setback:=", "0mm"
            ]
        ])

    def _fillets(self, radius, vertices, obj):
        self._modeler.Fillet(["NAME:Selections", "Selections:=", obj], [
            "NAME:Parameters",
            [
                "NAME:FilletParameters", "Edges:=", [], "Vertices:=", vertices,
                "Radius:=", radius, "Setback:=", "0mm"
            ]
        ])

    def _sweep_along_path(self, to_sweep, path_obj):
        """
        Adds thickness to path_obj by extending to a new dimension.
        to_sweep acts as a putty knife that determines the thickness.

        Args:
            to_sweep (polyline): Small polyline running perpendicular to path_obj
                                    whose length is the desired resulting thickness
            path_obj (polyline): Original polyline; want to broaden this
        """
        self.rename_obj(path_obj, str(path_obj) + '_path')
        new_name = self.rename_obj(to_sweep, path_obj)
        names = [path_obj, str(path_obj) + '_path']
        self._modeler.SweepAlongPath(self._selections_array(*names), [
            "NAME:PathSweepParameters", "DraftAngle:=", "0deg", "DraftType:=",
            "Round", "CheckFaceFaceIntersection:=", False, "TwistAngle:=",
            "0deg"
        ])
        return Polyline(new_name, self)

    def sweep_along_vector(self, names, vector):
        self._modeler.SweepAlongVector(self._selections_array(*names), [
            "NAME:VectorSweepParameters", "DraftAngle:=", "0deg",
            "DraftType:=", "Round", "CheckFaceFaceIntersection:=", False,
            "SweepVectorX:=", vector[0], "SweepVectorY:=", vector[1],
            "SweepVectorZ:=", vector[2]
        ])

    def rename_obj(self, obj, name):
        self._modeler.ChangeProperty([
            "NAME:AllTabs",
            [
                "NAME:Geometry3DAttributeTab", ["NAME:PropServers",
                                                str(obj)],
                ["NAME:ChangedProps", ["NAME:Name", "Value:=",
                                       str(name)]]
            ]
        ])
        return name


class ModelEntity(str, HfssPropertyObject):
    prop_tab = "Geometry3DCmdTab"
    model_command = None
    transparency = make_float_prop("Transparent",
                                   prop_tab="Geometry3DAttributeTab",
                                   prop_server=lambda self: self)
    material = make_str_prop("Material",
                             prop_tab="Geometry3DAttributeTab",
                             prop_server=lambda self: self)
    wireframe = make_float_prop("Display Wireframe",
                                prop_tab="Geometry3DAttributeTab",
                                prop_server=lambda self: self)
    coordinate_system = make_str_prop("Coordinate System")

    def __new__(self, val, *args, **kwargs):
        return str.__new__(self, val)

    def __init__(self, val, modeler):
        """
        :type val: str
        :type modeler: HfssModeler
        """
        super(ModelEntity,
              self).__init__()  # val) #Comment out keyword to match arguments
        self.modeler = modeler
        self.prop_server = self + ":" + self.model_command + ":1"


class Box(ModelEntity):
    model_command = "CreateBox"
    position = make_float_prop("Position")
    x_size = make_float_prop("XSize")
    y_size = make_float_prop("YSize")
    z_size = make_float_prop("ZSize")

    def __init__(self, name, modeler, corner, size):
        """
        :type name: str
        :type modeler: HfssModeler
        :type corner: [(VariableString, VariableString, VariableString)]
        :param size: [(VariableString, VariableString, VariableString)]
        """
        super(Box, self).__init__(name, modeler)
        self.modeler = modeler
        self.prop_holder = modeler._modeler
        self.corner = corner
        self.size = size
        self.center = [c + s / 2 for c, s in zip(corner, size)]
        faces = modeler.get_face_ids(self)
        self.z_back_face, self.z_front_face = faces[0], faces[1]
        self.y_back_face, self.y_front_face = faces[2], faces[4]
        self.x_back_face, self.x_front_face = faces[3], faces[5]


class Rect(ModelEntity):
    model_command = "CreateRectangle"

    # TODO: Add a rotated rectangle object.
    # Will need to first create rect, then apply rotate operation.

    def __init__(self, name, modeler, corner, size):
        super(Rect, self).__init__(name, modeler)
        self.prop_holder = modeler._modeler
        self.corner = corner
        self.size = size
        self.center = [c + s / 2 if s else c for c, s in zip(corner, size)]

    def make_center_line(self, axis):
        '''
        Returns `start` and `end` list of 3 coordinates
        '''
        axis_idx = ["x", "y", "z"].index(axis.lower())
        start = [c for c in self.center]
        start[axis_idx] -= self.size[axis_idx] / 2
        start = [self.modeler.eval_expr(s) for s in start]
        end = [c for c in self.center]
        end[axis_idx] += self.size[axis_idx] / 2
        end = [self.modeler.eval_expr(s) for s in end]
        return start, end

    def make_rlc_boundary(self, axis, r=0, l=0, c=0, name="LumpRLC"):
        start, end = self.make_center_line(axis)
        self.modeler._make_lumped_rlc(r,
                                      l,
                                      c,
                                      start,
                                      end, ["Objects:=", [self]],
                                      name=name)

    def make_lumped_port(self, axis, z0="50ohm", name="LumpPort"):
        start, end = self.make_center_line(axis)
        self.modeler._make_lumped_port(start,
                                       end, ["Objects:=", [self]],
                                       z0=z0,
                                       name=name)


class Polyline(ModelEntity):
    '''
    Assume closed polyline, which creates a polygon.
    '''

    model_command = "CreatePolyline"

    def __init__(self, name, modeler, points=None):
        super(Polyline, self).__init__(name, modeler)
        self.prop_holder = modeler._modeler
        if points is not None:
            self.points = points
            self.n_points = len(points)
        else:
            pass
            # TODO: points = collection of points


#        axis = find_orth_axis()

# TODO: find the plane of the polyline for now, assume Z
#    def find_orth_axis():
#        X, Y, Z = (True, True, True)
#        for point in points:
#            X =

    def unite(self, list_other):
        union = self.modeler.unite(self + list_other)
        return Polyline(union, self.modeler)

    def make_center_line(self, axis):  # Expects to act on a rectangle...
        # first : find center and size
        center = [0, 0, 0]

        for point in self.points:
            center = [
                center[0] + point[0] / self.n_points,
                center[1] + point[1] / self.n_points,
                center[2] + point[2] / self.n_points
            ]
        size = [
            2 * (center[0] - self.points[0][0]),
            2 * (center[1] - self.points[0][1]),
            2 * (center[1] - self.points[0][2])
        ]
        axis_idx = ["x", "y", "z"].index(axis.lower())
        start = [c for c in center]
        start[axis_idx] -= size[axis_idx] / 2
        start = [
            self.modeler.eval_var_str(s, unit=LENGTH_UNIT) for s in start
        ]  # TODO
        end = [c for c in center]
        end[axis_idx] += size[axis_idx] / 2
        end = [self.modeler.eval_var_str(s, unit=LENGTH_UNIT) for s in end]
        return start, end

    def make_rlc_boundary(self, axis, r=0, l=0, c=0, name="LumpRLC"):
        name = str(self) + '_' + name
        start, end = self.make_center_line(axis)
        self.modeler._make_lumped_rlc(r,
                                      l,
                                      c,
                                      start,
                                      end, ["Objects:=", [self]],
                                      name=name)

    def fillet(self, radius, vertex_index):
        self.modeler._fillet(radius, vertex_index, self)

    def vertices(self):
        return self.modeler.get_vertex_ids(self)

    def rename(self, new_name):
        '''
            Warning: The increment_name only works if the sheet has not been stracted or used as a tool elsewhere.
            These names are not checked; they require modifying get_objects_in_group.

        '''
        new_name = increment_name(
            new_name, self.modeler.get_objects_in_group(
                "Sheets"))  # this is for a closed polyline

        # check to get the actual new name in case there was a substracted object with that name
        face_ids = self.modeler.get_face_ids(str(self))
        self.modeler.rename_obj(self, new_name)  # now rename
        if len(face_ids) > 0:
            new_name = self.modeler.get_object_name_by_face_id(face_ids[0])
        return Polyline(str(new_name), self.modeler)


class OpenPolyline(ModelEntity):  # Assume closed polyline
    model_command = "CreatePolyline"
    show_direction = make_prop('Show Direction',
                               prop_tab="Geometry3DAttributeTab",
                               prop_server=lambda self: self)

    def __init__(self, name, modeler, points=None):
        super(OpenPolyline, self).__init__(name, modeler)
        self.prop_holder = modeler._modeler
        if points is not None:
            self.points = points
            self.n_points = len(points)
        else:
            pass


#        axis = find_orth_axis()

# TODO: find the plane of the polyline for now, assume Z
#    def find_orth_axis():
#        X, Y, Z = (True, True, True)
#        for point in points:
#            X =

    def vertices(self):
        return self.modeler.get_vertex_ids(self)

    def fillet(self, radius, vertex_index):
        self.modeler._fillet(radius, vertex_index, self)

    def fillets(self, radius, do_not_fillet=[]):
        '''
            do_not_fillet : Index list of vertices to not fillete
        '''
        raw_list_vertices = self.modeler.get_vertex_ids(self)
        list_vertices = []
        for vertex in raw_list_vertices[1:-1]:  # ignore the start and finish
            list_vertices.append(int(vertex))
        list_vertices = list(
            map(
                int,
                np.delete(list_vertices,
                          np.array(do_not_fillet, dtype=int) - 1)))
        #print(list_vertices, type(list_vertices[0]))
        if len(list_vertices) != 0:
            self.modeler._fillets(radius, list_vertices, self)
        else:
            pass

    def sweep_along_path(self, to_sweep):
        return self.modeler._sweep_along_path(to_sweep, self)

    def rename(self, new_name):
        '''
            Warning: The  increment_name only works if the sheet has not been stracted or used as a tool elsewher.
            These names are not checked - They require modifying get_objects_in_group
        '''
        new_name = increment_name(new_name,
                                  self.modeler.get_objects_in_group("Lines"))
        # , self.points)
        return OpenPolyline(self.modeler.rename_obj(self, new_name),
                            self.modeler)

    def copy(self, new_name):
        new_obj = OpenPolyline(self.modeler.copy(self), self.modeler)
        return new_obj.rename(new_name)


class HfssFieldsCalc(COMWrapper):
    def __init__(self, setup):
        """
        :type setup: HfssSetup
        """
        self.setup = setup
        super(HfssFieldsCalc, self).__init__()
        self.parent = setup
        self.Mag_E = NamedCalcObject("Mag_E", setup)
        self.Mag_H = NamedCalcObject("Mag_H", setup)
        self.Mag_Jsurf = NamedCalcObject("Mag_Jsurf", setup)
        self.Mag_Jvol = NamedCalcObject("Mag_Jvol", setup)
        self.Vector_E = NamedCalcObject("Vector_E", setup)
        self.Vector_H = NamedCalcObject("Vector_H", setup)
        self.Vector_Jsurf = NamedCalcObject("Vector_Jsurf", setup)
        self.Vector_Jvol = NamedCalcObject("Vector_Jvol", setup)
        self.ComplexMag_E = NamedCalcObject("ComplexMag_E", setup)
        self.ComplexMag_H = NamedCalcObject("ComplexMag_H", setup)
        self.ComplexMag_Jsurf = NamedCalcObject("ComplexMag_Jsurf", setup)
        self.ComplexMag_Jvol = NamedCalcObject("ComplexMag_Jvol", setup)
        self.P_J = NamedCalcObject("P_J", setup)

        self.named_expression = {
        }  # dictionary to hold additional named expressions

    def clear_named_expressions(self):
        self.parent.parent._fields_calc.ClearAllNamedExpr()

    def declare_named_expression(self, name):
        """"
        If a named expression has been created in the fields calculator, this
        function can be called to initialize the name to work with the fields object
        """
        self.named_expression[name] = NamedCalcObject(name, self.setup)

    def use_named_expression(self, name):
        """
        Expression can be used to access dictionary of named expressions,
        Alternately user can access dictionary directly via named_expression()
        """
        return self.named_expression[name]


class CalcObject(COMWrapper):
    def __init__(self, stack, setup):
        """
        :type stack: [(str, str)]
        :type setup: HfssSetup
        """
        super(CalcObject, self).__init__()
        self.stack = stack
        self.setup = setup
        self.calc_module = setup.parent._fields_calc

    def _bin_op(self, other, op):
        if isinstance(other, (int, float)):
            other = ConstantCalcObject(other, self.setup)

        stack = self.stack + other.stack
        stack.append(("CalcOp", op))
        return CalcObject(stack, self.setup)

    def _unary_op(self, op):
        stack = self.stack[:]
        stack.append(("CalcOp", op))
        return CalcObject(stack, self.setup)

    def __add__(self, other):
        return self._bin_op(other, "+")

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self._bin_op(other, "-")

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        return self._bin_op(other, "*")

    def __rmul__(self, other):
        return self * other

    def __div__(self, other):
        return self._bin_op(other, "/")

    def __rdiv__(self, other):
        other = ConstantCalcObject(other, self.setup)
        return other / self

    def __pow__(self, other):
        return self._bin_op(other, "Pow")

    def dot(self, other):
        return self._bin_op(other, "Dot")

    def __neg__(self):
        return self._unary_op("Neg")

    def __abs__(self):
        return self._unary_op("Abs")

    def __mag__(self):
        return self._unary_op("Mag")

    def mag(self):
        return self._unary_op("Mag")

    def smooth(self):
        return self._unary_op("Smooth")

    def conj(self):
        return self._unary_op("Conj")  # make this right

    def scalar_x(self):
        return self._unary_op("ScalarX")

    def scalar_y(self):
        return self._unary_op("ScalarY")

    def scalar_z(self):
        return self._unary_op("ScalarZ")

    def norm_2(self):

        return (self.__mag__()).__pow__(2)
        # return self._unary_op("ScalarX")**2+self._unary_op("ScalarY")**2+self._unary_op("ScalarZ")**2

    def real(self):
        return self._unary_op("Real")

    def imag(self):
        return self._unary_op("Imag")

    def complexmag(self):
        return self._unary_op("CmplxMag")

    def _integrate(self, name, type):
        stack = self.stack + [(type, name), ("CalcOp", "Integrate")]
        return CalcObject(stack, self.setup)

    def _maximum(self, name, type):
        stack = self.stack + [(type, name), ("CalcOp", "Maximum")]
        return CalcObject(stack, self.setup)

    def getQty(self, name):
        stack = self.stack + [("EnterQty", name)]
        return CalcObject(stack, self.setup)

    def integrate_line(self, name):
        return self._integrate(name, "EnterLine")

    def normal2surface(self, name):
        ''' return the part normal to surface.
            Complex Vector. '''
        stack = self.stack + [("EnterSurf", name),
                                ("CalcOp",    "Normal")]
        stack.append(("CalcOp", "Dot"))
        stack.append(("EnterSurf", name))
        stack.append(("CalcOp",    "Normal"))
        stack.append(("CalcOp", "*"))
        return CalcObject(stack, self.setup)

    def tangent2surface(self, name):
        ''' return the part tangent to surface.
            Complex Vector. '''
        stack = self.stack + [("EnterSurf", name),
                                ("CalcOp",    "Normal")]
        stack.append(("CalcOp", "Dot"))
        stack.append(("EnterSurf", name))
        stack.append(("CalcOp",    "Normal"))
        stack.append(("CalcOp", "*"))
        stack = self.stack + stack
        stack.append(("CalcOp", "-"))
        return CalcObject(stack, self.setup)

    def integrate_line_tangent(self, name):
        ''' integrate line tangent to vector expression \n
            name = of line to integrate over '''
        self.stack = self.stack + [("EnterLine", name), ("CalcOp", "Tangent"),
                                   ("CalcOp", "Dot")]
        return self.integrate_line(name)

    def line_tangent_coor(self, name, coordinate):
        ''' integrate line tangent to vector expression \n
            name = of line to integrate over '''
        if coordinate not in ['X', 'Y', 'Z']:
            raise ValueError
        self.stack = self.stack + [("EnterLine", name), ("CalcOp", "Tangent"),
                                   ("CalcOp", "Scalar" + coordinate)]
        return self.integrate_line(name)

    def integrate_surf(self, name="AllObjects"):
        return self._integrate(name, "EnterSurf")

    def integrate_vol(self, name="AllObjects"):
        return self._integrate(name, "EnterVol")

    def maximum_vol(self, name='AllObjects'):
        return self._maximum(name, 'EnterVol')

    def times_eps(self):
        stack = self.stack + [("ClcMaterial", ("Permittivity (epsi)", "mult"))]
        return CalcObject(stack, self.setup)

    def times_mu(self):
        stack = self.stack + [("ClcMaterial", ("Permeability (mu)", "mult"))]
        return CalcObject(stack, self.setup)

    def write_stack(self):
        for fn, arg in self.stack:
            if np.size(arg) > 1 and fn not in ['EnterVector']:
                getattr(self.calc_module, fn)(*arg)
            else:
                getattr(self.calc_module, fn)(arg)

    def save_as(self, name):
        """if the object already exists, try clearing your
        named expressions first with fields.clear_named_expressions"""
        self.write_stack()
        self.calc_module.AddNamedExpr(name)
        return NamedCalcObject(name, self.setup)

    def evaluate(self, phase=0, lv=None, print_debug=False):  # , n_mode=1):
        self.write_stack()
        if print_debug:
            print('---------------------')
            print('writing to stack: OK')
            print('-----------------')
        #self.calc_module.set_mode(n_mode, 0)
        setup_name = self.setup.solution_name

        if lv is not None:
            args = lv
        else:
            args = []

        args.append("Phase:=")
        args.append(str(int(phase)) + "deg")

        if isinstance(self.setup, HfssDMSetup):
            args.extend(["Freq:=", self.setup.solution_freq])

        self.calc_module.ClcEval(setup_name, args)
        return float(self.calc_module.GetTopEntryValue(setup_name, args)[0])


class NamedCalcObject(CalcObject):
    def __init__(self, name, setup):
        self.name = name
        stack = [("CopyNamedExprToStack", name)]
        super(NamedCalcObject, self).__init__(stack, setup)


class ConstantCalcObject(CalcObject):
    def __init__(self, num, setup):
        stack = [("EnterScalar", num)]
        super(ConstantCalcObject, self).__init__(stack, setup)


class ConstantVecCalcObject(CalcObject):
    def __init__(self, vec, setup):
        stack = [("EnterVector", vec)]
        super(ConstantVecCalcObject, self).__init__(stack, setup)


def get_active_project():
    ''' If you see the error:
        "The requested operation requires elevation."
        then you need to run your python as an admin.
    '''
    import ctypes
    import os
    try:
        is_admin = os.getuid() == 0
    except AttributeError:
        is_admin = ctypes.windll.shell32.IsUserAnAdmin() != 0
    if not is_admin:
        print('\033[93m WARNING: you are not running as an admin! \
            You need to run as an admin. You will probably get an error next.\
                 \033[0m')

    app = HfssApp()
    desktop = app.get_app_desktop()
    return desktop.get_active_project()


def get_active_design():
    project = get_active_project()
    return project.get_active_design()


def get_report_arrays(name: str):
    d = get_active_design()
    r = HfssReport(d, name)
    return r.get_arrays()


def load_ansys_project(proj_name: str,
                       project_path: str = None,
                       extension: str = '.aedt'):
    '''
    Utility function to load an Ansys project.

    Args:
        proj_name : None  --> get active. (make sure 2 run as admin)
        extension : `aedt` is for 2016 version and newer
    '''
    if project_path:
        # convert slashes correctly for system
        project_path = Path(project_path)

        # Checks
        assert project_path.is_dir(
        ), "ERROR! project_path is not a valid directory \N{loudly crying face}.\
            Check the path, and especially \\ characters."

        project_path /= project_path / Path(proj_name + extension)

        if (project_path).is_file():
            logger.info('\tFile path to HFSS project found.')
        else:
            raise Exception(
                "ERROR! Valid directory, but invalid project filename. \N{loudly crying face} Not found!\
                     Please check your filename.\n%s\n" % project_path)

        if (project_path / '.lock').is_file():
            logger.warning(
                '\t\tFile is locked. \N{fearful face} If connection fails, delete the .lock file.'
            )

    app = HfssApp()
    logger.info("\tOpened Ansys App")

    desktop = app.get_app_desktop()
    logger.info(f"\tOpened Ansys Desktop v{desktop.get_version()}")
    #logger.debug(f"\tOpen projects: {desktop.get_project_names()}")

    if proj_name is not None:
        if proj_name in desktop.get_project_names():
            desktop.set_active_project(proj_name)
            project = desktop.get_active_project()
        else:
            project = desktop.open_project(str(project_path))
    else:
        projects_in_app = desktop.get_projects()
        if projects_in_app:
            project = desktop.get_active_project()
        else:
            project = None

    if project:
        logger.info(
            f"\tOpened Ansys Project\n\tFolder:    {project.get_path()}\n\tProject:   {project.name}"
        )
    else:
        logger.info(f"\tAnsys Project was not found.\n\t Project is None.")

    return app, desktop, project
