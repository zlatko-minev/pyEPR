'''
 Copyright Zlatko Minev and Zaki Leghtas
 2015, 2016, 2017, 2018
'''
from __future__ import print_function    # Python 2.7 and 3 compatibility
import os
import sys
import time
import shutil
#import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
# Standard imports
from numpy        import pi
from numpy.linalg import inv
from stat         import S_ISREG, ST_CTIME, ST_MODE
from pandas       import HDFStore, Series, DataFrame
from collections  import OrderedDict
from pathlib      import Path

# pyEPR custom imports
from . import hfss
from . import logger
from . import config
from . import AttrDict
from .hfss        import ureg, CalcObject, ConstantVecCalcObject, set_property
from .toolbox     import print_NoNewLine, print_color, deprecated, fact, epsilon_0, hbar, Planck, fluxQ, nck, \
                         divide_diagonal_by_2, print_matrix, DataFrame_col_diff, get_instance_vars,\
                         sort_df_col, sort_Series_idx
from .toolbox_circuits import Calcs_basic
from .toolbox_plotting import cmap_discrete, legend_translucent
from .numeric_diag import bbq_hmt, make_dispersive
import matplotlib as mpl
from .toolbox_report import plot_convergence_f_vspass, plot_convergence_max_df, plot_convergence_solved_elem, plot_convergence_maxdf_vs_sol


class Project_Info(object):
    """
    Class containing options and information about the manipulation and analysis in HFSS.

    Junction info:
    -----------------------
        self.junctions : OrderedDict()

        A Josephson tunnel junction has to have its parameters specified here for the analysis.
        Each junction is given a name and is specified by a dictionary.
        It has the following properties:

            1. `Lj_variable` : Name of HFSS variable that specifies junction inductance Lj defined on the boundary condition in HFSS. DO NOT USE Global names that start with $.
            2. `rect`        : Name of HFSS rectangle on which lumped boundary condition is specified.
            3. `line`        : Name of HFSS polyline which spans the length of the recntalge. Used to define the voltage across the junction.  Used to define the current orientation for each junction. Used to define sign of ZPF.
            4. `length`      : Length in HFSS of the junction rectangle and line (specified in meters).

        Example definition:

        ..code-block python

            # Define a single junction
            pinfo = Project_Info('')
            pinfo.junctions['j1'] = {'Lj_variable' : 'Lj1',
                         'rect'        : 'JJrect1',
                         'line'        : 'JJline1',
                         'length'      : parse_units('50um')} # Length is in meters

            # Specify multiple junctions in HFSS model
            n_junctions = 5
            for i in range(1, 1+n_junctions):
                pinfo.junctions[f'j{i}'] = {'Lj_variable' : f'Lj{i}',
                                            'rect'        : f'JJrect{i}',
                                            'line'        : f'JJline{i}',
                                            'length'      : parse_units('50um')}

    HFSS app connection settings
    -----------------------
        project_path  : str
            Directory path to the hfss project file. Should be the directory, not the file.
            default = None: Assumes the project is open, and thus gets the project based on `project_name`
        project_name  : str,  None
            Name of the project within the project_path. "None" will get the current active one.
        design_name   : str,  None
            Name of the design within the project. "None" will get the current active one.
        setup_name    : str,  None
            Name of the setup within the design. "None" will get the current active one.


    Additional init setting:
    -----------------------
        do_connect : True by default. Connect to HFSS


    HFSS desgin settings
    -----------------------
    describe junction parameters
        junc_rects    = None
            Name of junction rectangles in HFSS
        junc_lines    = None
            Name of lines in HFSS used to define the current orientation for each junction
        junc_LJ_names = None
            Name of junction inductance variables in HFSS.
            Note, DO NOT USE Global names that start with $.
        junc_lens     = None
            Junciton rect. length, measured in meters.

    """

    class _Dissipative:
        #TODO: remove and turn to dict

        def __init__(self):
            self.dielectrics_bulk    = None
            self.dielectric_surfaces = None
            self.resistive_surfaces  = None
            self.seams               = None

    def __init__(self, project_path=None, project_name=None, design_name=None,
<<<<<<< HEAD
                 setup_name= None, do_connect = True):
=======
                setup_name=None, do_connect = True):
>>>>>>> origin/master

        self.project_path  = str(Path(project_path)) if not (project_path is None) else None # Path: format path correctly to system convention
        self.project_name  = project_name
        self.design_name   = design_name
<<<<<<< HEAD
        self.setup_name    =  setup_name
=======
        self.setup_name    = setup_name
>>>>>>> origin/master

        ## HFSS desgin: describe junction parameters
        # TODO: introduce modal labels
        self.junctions     = OrderedDict() # See above for help
        self.ports         = OrderedDict()

        ## Dissipative HFSS volumes and surfaces
        self.dissipative   = self._Dissipative()
        self.options       = config.options_hfss

        # Conected to HFSS variable
        self.app           = None
        self.desktop       = None
        self.project       = None
        self.design        = None
        self.setup         = None

        if do_connect:
            self.connect()


    _Forbidden = ['app', 'design', 'desktop', 'project',
                  'dissipative', 'setup', '_Forbidden', 'junctions']

    def save(self, hdf):
        '''
            hdf : pd.HDFStore
        '''
        hdf['project_info']           = pd.Series(get_instance_vars(self, self._Forbidden))
        hdf['project_info_dissip']    = pd.Series(get_instance_vars(self.dissipative))
        hdf['project_info_options']   = pd.Series(get_instance_vars(self.options))
        hdf['project_info_junctions'] = pd.DataFrame(self.junctions)
        hdf['project_info_ports']     = pd.DataFrame(self.ports)

    @deprecated
    def connect_to_project(self):
        return self.connect()

    def connect(self):
        '''
        Connect to HFSS design.
        '''
        #logger.info('Connecting to HFSS ...')

        self.app, self.desktop, self.project = hfss.load_ansys_project(
            self.project_name, self.project_path)
        self.project_name = self.project.name
        self.project_path = self.project.get_path()

        # Design
        if self.design_name is None:
            self.design = self.project.get_active_design()
            self.design_name = self.design.name
            logger.info(f'\tOpened active design\n\tDesign:    {self.design_name} [Solution type: {self.design.solution_type}]')
        else:
            try:
                self.design = self.project.get_design(self.design_name)
            except Exception as e:
                tb = sys.exc_info()[2]
                logger.error(f"Original error: {e}\n")
                raise(Exception(' Did you provide the correct design name? Failed to pull up design.').with_traceback(tb))

        #if not ('Eigenmode' == self.design.solution_type):
        #    logger.warning('\tWarning: The design tpye is not Eigenmode. Are you sure you dont want eigenmode?')

        # Setup
        try:
            n_setups = len(self.design.get_setup_names())
            if n_setups == 0:
                logger.warning('\tNo design setup detected.')
                if self.design.solution_type == 'Eigenmode':
                    logger.warning('\tCreating eigenmode default setup one.')
                    setup = self.design.create_em_setup()
                    self.setup_name = setup.name
                elif self.design.solution_type == 'DrivenModal':
                    setup = self.design.create_dm_setup() # adding a driven modal design
                    self.setup_name = setup.name
            self.get_setup(self.setup_name) # get the actual setup if there is one

        except Exception as e:
            tb = sys.exc_info()[2]
            logger.error(f"Original error: {e}\n")
            raise(Exception(' Did you provide the correct setup name? Failed to pull up setup.').with_traceback(tb))

        # Finalize
        self.project_name = self.project.name
        self.design_name  = self.design.name

        logger.info('\tConnected successfully.\t :)\t :)\t :)\t\n')

        return self

    def get_setup(self, name):
        """
        Connects to a specific setup for the design. 
        Sets  self.setup and self.setup_name.
        If Name is 
        """

        if name is None:
            return None
        else:
            self.setup = self.design.get_setup(name=self.setup_name)
            if self.setup is None:
                 logger.error(f"Could not retrieve setup: {self.setup_name}\n \
                     Did you give the right name? Does it exist?")
            self.setup_name = self.setup.name
            logger.info(f'\tOpened setup: name={self.setup_name}  type={type(self.setup)}')
            return self.setup
        

    def check_connected(self):
        """Checks if fully connected including setup
        """
        return\
        (self.setup   is not None) and\
        (self.design  is not None) and\
        (self.project is not None) and\
        (self.desktop is not None) and\
        (self.app     is not None)

    def disconnect(self):
        '''
        Disconnect from existing HFSS design.
        '''
        assert self.check_connected(
        ) is True, "it does not appear that you have connected to HFSS yet. use connect()"
        self.project.release()
        self.desktop.release()
        self.app.release()
        hfss.release()

    ### UTILITY FUNCTIONS

    def get_dm(self):
        '''
        Get the design and modeler

        .. code-block:: python
            oDesign, oModeler = projec.get_dm()

        '''
        oDesign = self.design
        oModeler = oDesign.modeler
        return oDesign, oModeler

    def get_all_variables_names(self):
        """Returns array of all project and local design names."""
        return self.project.get_variable_names() + self.design.get_variable_names()

    def get_all_object_names(self):
        """Returns array of strings"""
        oObjects = []
        for s in ["Non Model", "Solids", "Unclassified", "Sheets", "Lines"]:
            oObjects += self.design.modeler.get_objects_in_group(s)
        return oObjects

    def validate_junction_info(self):
        """ Validate that the user has put in the junction info correctly.
            Do no also forget to check the length of the rectangles/line of
            the junction if you change it.
        """
        all_variables_names = self.get_all_variables_names()
        all_object_names    = self.get_all_object_names()
        for jjnm, jj in self.junctions.items():
            assert jj['Lj_variable'] in all_variables_names, "pyEPR project_info user error found: Seems like for junction `%s` you specified a design or project variable for `Lj_variable` that does not exist in HFSS by the name: `%s` " % (
                jjnm, jj['Lj_variable'])
            for name in ['rect', 'line']:
                assert jj[name] in all_object_names, "pyEPR project_info user error found: Seems like for junction `%s` you specified a %s that does not exist in HFSS by the name: `%s` " % (
                    jjnm, name, jj[name])
        #TODO: Check the length of the rectnagle


#==============================================================================
#%% Main compuation class & interface with HFSS
#==============================================================================
class pyEPR_HFSS(object):
    """
    This class defines a pyEPR_HFSS object which calculates and saves
    Hamiltonian parameters from an HFSS simulation.
    Further, it allows one to calcualte dissipation, etc
    """

    def __init__(self, *args, **kwargs):
        '''
        Parameters:
        -------------------
            project_info    : Project_Info
                Suplpy the project info or the parameters to create pinfo

        Example use:
        -------------------
        '''
        if (len(args) == 1) and (args[0].__class__.__name__ == 'Project_Info'): #isinstance(args[0], Project_Info): # fails on module repload with changes
            project_info = args[0]
        else:
            assert len(args) == 0, 'Since you did not pass a Project_info object as a arguemnt, we now assuem you are trying to create a project info object here by apassing its arguments. See Project_Info. It does not take any arguments, only kwargs.'
            project_info = Project_Info(*args, **kwargs)


        # Input
        self.pinfo = project_info
        if self.pinfo.check_connected() is False:
            self.pinfo.connect()

        self.verbose          = True #TODO: change verbose to logger. remove verbose flags
        self.append_analysis  = False #TODO

        # hfss connect module
        self.fields           = None
        self.solutions        = None
        if self.setup:
            self.fields       = self.setup.get_fields()
            self.solutions    = self.setup.get_solutions()

        # Variations - the following get updated in update_variation_information
        self.nmodes           = int(1)
        self.listvariations   = ("",)
        self.nominalvariation = '0'
        self.nvariations      = 0
        self.update_variation_information()

        self.hfss_variables   = OrderedDict()  # container for eBBQ list of varibles


        if self.verbose:
            print('Design \"%s\" info:'%self.design.name)
            print('\t%-15s %d\n\t%-15s %d' %('# eigenmodes', self.nmodes, '# variations', self.nvariations))

        # Setup data saving
        self.setup_data()
        self.latest_h5_path = None # #self.get_latest_h5()
        '''         #TODO: to be implemented to use old files
        if self.latest_h5_path is not None and self.append_analysis:
            latest_bbq_analysis = pyEPR_Analysis(self.latest_h5_path)
            if self.verbose:
                print( 'Varied variables and values : ', latest_bbq_analysis.get_swept_variables(), \
                       'Variations : ', latest_bbq_analysis.variations)
            '''

    @property
    def setup(self):
        return self.pinfo.setup

    @property
    def design(self):
        return self.pinfo.design

    @property
    def project(self):
        return self.pinfo.project

    @property
    def desktop(self):
        return self.pinfo.desktop

    @property
    def app(self):
        return self.pinfo.app

    @property
    def junctions(self):
        return self.pinfo.junctions

    @property
    def ports(self):
        return self.pinfo.ports

    @property
    def options(self):
        return self.pinfo.options

    def get_latest_h5(self):
        '''
            No longer used. Could be added back in.
        '''
        dirpath = self.data_dir

        entries1 = (os.path.join(dirpath, fn) for fn in os.listdir(dirpath))     # get all entries in the directory w/ stats
        entries2 = ((os.stat(path), path) for path in entries1)
        entries3 = ((stat[ST_CTIME], path)                                       # leave only regular files, insert creation date
                   for stat, path in entries2 if S_ISREG(stat[ST_MODE]) and path[-4:]=='hdf5')
        #NOTE: on Windows `ST_CTIME` is a creation date but on Unix it could be something else
        #NOTE: use `ST_MTIME` to sort by a modification date

        paths_sorted = []
        for cdate, path in sorted(entries3):
            paths_sorted.append(path)
            #print time.ctime(cdate), os.path.basename(path)

        if len(paths_sorted) > 0:
            self.latest_h5_path = paths_sorted[-1]
            if self.verbose:
                print('This simulations has been analyzed, latest data in ' + self.latest_h5_path)

        else:
            self.latest_h5_path = None
            if self.verbose:
                print('This simulation has never been analyzed')

    def setup_data(self):
        '''
        Set up folder paths for saving data to.
        '''
        data_dir = Path(config.root_dir) / \
            Path(self.project.name)/Path(self.design.name)

        #if self.verbose:
        #    print("\nResults will be saved to:\n" +'-  '*20+'\n\t'+ str(data_dir)+'\n'+'-  '*20+'\n')
        if len(self.design.name) > 50:
            print_color('WARNING!   DESING FILENAME MAY BE TOO LONG! ')

        if not data_dir.is_dir():
            data_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir = str(data_dir)
        self.data_filename = str(
            data_dir / (time.strftime('%Y-%m-%d %H-%M-%S', time.localtime()) + '.hdf5'))

    """
    @deprecated
    def calc_p_j(self, modes=None, variation=None):
        '''
        Calculates the p_j for all the modes.
        Requires a calculator expression called P_J.
        '''
        lv = self.get_lv(variation)
        if modes is None:
            modes = range(self.nmodes)

        pjs = OrderedDict()
        for ii, m in enumerate(modes):
            print('Calculating p_j for mode ' + str(m) + ' (' + str(ii) + '/' + str(np.size(modes)-1) + ')')
            self.solutions.set_mode(m+1, 0)
            self.fields = self.setup.get_fields()
            P_J = self.fields.P_J
            pjs['pj_'+str(m)] = P_J.evaluate(lv=lv)
        self.pjs = pjs
        if self.verbose:
            print(pjs)
        return pjs
    """

    def calc_p_junction_single(self, mode):
        '''
        This function is used in the case of a single junction only.
        For multiple junctions, see `calc_p_junction`.

        Assumes no lumped capacitive elements.
        '''
        pj = OrderedDict()
        pj_val = (self.U_E-self.U_H)/self.U_E
        pj['pj_'+str(mode)] = np.abs(pj_val)
        print('    p_j_' + str(mode) + ' = ' + str(pj_val))
        return pj

    #TODO: replace this method with the one below, here because osme funcs use it still
    def get_freqs_bare(self, variation):
        #str(self.get_lv(variation))
        freqs_bare_vals = []
        freqs_bare_dict = OrderedDict()
        freqs, kappa_over_2pis = self.solutions.eigenmodes(
            self.get_lv_EM(variation))
        for m in range(self.nmodes):
            freqs_bare_dict['freq_bare_'+str(m)] = 1e9*freqs[m]
            freqs_bare_vals.append(1e9*freqs[m])
            if kappa_over_2pis is not None:
                freqs_bare_dict['Q_'+str(m)] = freqs[m]/kappa_over_2pis[m]
            else:
                freqs_bare_dict['Q_'+str(m)] = 0
        self.freqs_bare = freqs_bare_dict
        self.freqs_bare_vals = freqs_bare_vals
        return freqs_bare_dict, freqs_bare_vals

    def get_freqs_bare_pd(self, variation):
        '''
            Retun pd.Sereis of modal freq and qs for given variation
        '''
        freqs, kappa_over_2pis = self.solutions.eigenmodes(
            self.get_lv_EM(variation))
        if kappa_over_2pis is None:
            kappa_over_2pis = np.zeros(len(freqs))
        freqs = pd.Series(freqs, index=range(len(freqs)))  # GHz
        Qs    = freqs / pd.Series(kappa_over_2pis, index=range(len(freqs)))
        return freqs, Qs

    def get_lv(self, variation=None):
        '''
        List of variation variables.
        Returns list of var names and var values.
        Such as ['Lj1:=','13nH', 'QubitGap:=','100um']

        Parameters
        -----------
            variation :  string number such as '0' or '1' or ...
        '''

        if variation is None:
            lv = self.nominalvariation
            lv = self.parse_listvariations(lv)
        else:
            lv = self.listvariations[ureg(variation)]
            lv = self.parse_listvariations(lv)
        return lv

    def get_lv_EM(self, variation):
        if variation is None:
            lv = self.nominalvariation
            #lv = self.parse_listvariations_EM(lv)
        else:
            lv = self.listvariations[ureg(variation)]
            #lv = self.parse_listvariations_EM(lv)
        return str(lv)

    def parse_listvariations_EM(self, lv):
        lv = str(lv)
        lv = lv.replace("=", ":=,")
        lv = lv.replace(' ', ',')
        lv = lv.replace("'", "")
        lv = lv.split(",")
        return lv

    def parse_listvariations(self, lv):
        lv = str(lv)
        lv = lv.replace("=", ":=,")
        lv = lv.replace(' ', ',')
        lv = lv.replace("'", "")
        lv = lv.split(",")
        return lv

    def get_variables(self, variation=None):
        lv = self.get_lv(variation)
        variables = OrderedDict()
        for ii in range(int(len(lv)/2)):
            variables['_'+lv[2*ii][:-2]] = lv[2*ii+1]
        self.variables = variables
        return variables


    def calc_energy_electric(self,
                             variation=None,
                             volume='AllObjects',
                             smooth=False):
        r'''
        Calculates two times the peak electric energy, or 4 times the RMS, :math:`4*\mathcal{E}_{\mathrm{elec}}`
        (since we do not divide by 2 and use the peak phasors).

        .. math::
            \mathcal{E}_{\mathrm{elec}}=\frac{1}{4}\mathrm{Re}\int_{V}\mathrm{d}v\vec{E}_{\text{max}}^{*}\overleftrightarrow{\epsilon}\vec{E}_{\text{max}}


        volume : string | 'AllObjects'
        smooth : bool | False
            Smooth the electric field or not when performing calculation

        Example use to calcualte the energy participation of a substrate

        .. code-block python
            ℰ_total  = epr_hfss.calc_energy_electric(volume='AllObjects')
            ℰ_substr = epr_hfss.calc_energy_electric(volume='Box1')
            print(f'Energy in substrate = {100*ℰ_substr/ℰ_total:.1f}%')

        '''

        calcobject = CalcObject([], self.setup)

        vecE = calcobject.getQty("E")
        if smooth:
            vecE = vecE.smooth()
        A = vecE.times_eps()
        B = vecE.conj()
        A = A.dot(B)
        A = A.real()
        A = A.integrate_vol(name=volume)

        lv = self.get_lv(variation)
        return A.evaluate(lv=lv)

    def calc_energy_magnetic(self,
                 variation=None,
                 volume='AllObjects',
                 smooth=True):
        '''
        See calc_energy_electric
        '''

        calcobject = CalcObject([], self.setup)

        vecH = calcobject.getQty("H")
        if smooth:
            vecH = vecH.smooth()
        A = vecH.times_mu()
        B = vecH.conj()
        A = A.dot(B)
        A = A.real()
        A = A.integrate_vol(name=volume)

        lv = self.get_lv(variation)
        return A.evaluate(lv=lv)

    def calc_p_electric_volume(self,
                    name_dielectric3D,
                    relative_to='AllObjects',
                    E_total=None
                   ):
        r'''
        Calculate the dielectric energy-participatio ratio
        of a 3D object (one that has volume) relative to the dielectric energy of
        a list of object objects.

        This is as a function relative to another object or all objects.

        When all objects are specified, this does not include any energy
        that might be stored in any lumped elements or lumped capacitors.

        Returns:
        ---------
            ℰ_object/ℰ_total, (ℰ_object, _total)
        '''

        if E_total is None:
            logger.debug('Calculating ℰ_total')
            ℰ_total  = self.calc_energy_electric(volume=relative_to)
        else:
            ℰ_total = E_total

        logger.debug('Calculating ℰ_object')
        ℰ_object = self.calc_energy_electric(volume=name_dielectric3D)

        return ℰ_object/ℰ_total, (ℰ_object, ℰ_total)


    def calc_current(self, fields, line):
        '''
        Function to calculate Current based on line. Not in use
        line : integration line between plates - name
        '''
        self.design.Clear_Field_Clac_Stack()
        comp = fields.Vector_H
        exp  = comp.integrate_line_tangent(line)
        I    = exp.evaluate(phase = 90)
        self.design.Clear_Field_Clac_Stack()
        return I

    def calc_avg_current_J_surf_mag(self, variation, junc_rect, junc_line):
        ''' Peak current I_max for mdoe J in junction J
            The avg. is over the surface of the junction. I.e., spatial. '''
        lv = self.get_lv(variation)

        jl, uj = self.get_junc_len_dir(variation, junc_line)

        uj = ConstantVecCalcObject(uj, self.setup)
        calc = CalcObject([], self.setup)
        #calc = calc.getQty("Jsurf").mag().integrate_surf(name = junc_rect)
        calc = (((calc.getQty("Jsurf")).dot(uj)).imag()
                ).integrate_surf(name=junc_rect)
        I = calc.evaluate(lv=lv) / jl  # phase = 90
        #self.design.Clear_Field_Clac_Stack()
        return I

    def calc_current_line_voltage(self, variation, junc_line_name, junc_L_Henries):
        '''
        Peak current I_max for prespecified mode calculating line voltage across junction.

        Parameters:
        ------------------------------------------------
            variation: variation number
            junc_line_name: name of the HFSS line spanning the junction
            junc_L_Henries: junction inductance in henries

            TODO: Smooth?
        '''
        lv = self.get_lv(variation)
        v_calc_real = CalcObject([], self.setup).getQty(
            "E").real().integrate_line_tangent(name=junc_line_name)
        v_calc_imag = CalcObject([], self.setup).getQty(
            "E").imag().integrate_line_tangent(name=junc_line_name)
        V = np.sqrt(v_calc_real.evaluate(lv=lv)**2 +
                    v_calc_imag.evaluate(lv=lv)**2)
        freq = CalcObject(
            [('EnterOutputVar', ('Freq', "Complex"))], self.setup).real().evaluate()
        return V/(2*np.pi*freq*junc_L_Henries)  # I=V/(wL)s

    def calc_line_current(self, variation, junc_line_name):
        lv = self.get_lv(variation)
        calc = CalcObject([], self.setup)
        calc = calc.getQty("H").imag().integrate_line_tangent(
            name=junc_line_name)
        #self.design.Clear_Field_Clac_Stack()
        return calc.evaluate(lv=lv)

    def get_junc_len_dir(self, variation, junc_line):
        '''
        Return the length and direction of a junction defined by a line

        Inputs: variation: simulation variation
                junc_line: polyline object

        Outputs: jl (float) junction length
                 uj (list of 3 floats) x,y,z coordinates of the unit vector
                 tangent to the junction line
        '''
        #
        lv = self.get_lv(variation)
        u = []
        for coor in ['X', 'Y', 'Z']:
            calc = CalcObject([], self.setup)
            calc = calc.line_tangent_coor(junc_line, coor)
            u.append(calc.evaluate(lv=lv))

        jl = float(np.sqrt(u[0]**2+u[1]**2+u[2]**2))
        uj = [float(u[0]/jl), float(u[1]/jl), float(u[2]/jl)]
        return jl, uj

    def get_Qseam(self, seam, mode, variation):
        r'''
        Caculate the contribution to Q of a seam, by integrating the current in
        the seam with finite conductance: set in the config file
        ref: http://arxiv.org/pdf/1509.01119.pdf
        '''

        lv = self.get_lv(variation)
        Qseam = OrderedDict()
        print('Calculating Qseam_' + seam + ' for mode ' + str(mode) +
              ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')
        # overestimating the loss by taking norm2 of j, rather than jperp**2
        j_2_norm = self.fields.Vector_Jsurf.norm_2()
        int_j_2 = j_2_norm.integrate_line(seam)
        int_j_2_val = int_j_2.evaluate(lv=lv, phase=90)
        yseam = int_j_2_val/self.U_H/self.omega

        Qseam['Qseam_'+seam+'_' +
              str(mode)] = config.Dissipation_params.gseam/yseam

        print('Qseam_' + seam + '_' + str(mode) + str(' = ') +
              str(config.Dissipation_params.gseam/config.Dissipation_params.yseam))

        return Series(Qseam)

    def get_Qseam_sweep(self, seam, mode, variation, variable, values, unit, pltresult=True):
        # values = ['5mm','6mm','7mm']
        # ref: http://arxiv.org/pdf/1509.01119.pdf

        self.solutions.set_mode(mode+1, 0)
        self.fields = self.setup.get_fields()
        freqs_bare_dict, freqs_bare_vals = self.get_freqs_bare(variation)
        self.omega = 2*np.pi*freqs_bare_vals[mode]
        print(variation)
        print(type(variation))
        print(ureg(variation))
        self.U_H = self.calc_energy_magnetic(variation)

        lv = self.get_lv(variation)
        Qseamsweep = []
        print('Calculating Qseam_' + seam + ' for mode ' + str(mode) +
              ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')
        for value in values:
            self.design.set_variable(variable, str(value)+unit)

            # overestimating the loss by taking norm2 of j, rather than jperp**2
            j_2_norm = self.fields.Vector_Jsurf.norm_2()
            int_j_2 = j_2_norm.integrate_line(seam)
            int_j_2_val = int_j_2.evaluate(lv=lv, phase=90)
            yseam = int_j_2_val/self.U_H/self.omega
            Qseamsweep.append(config.Dissipation_params.gseam/yseam)
#        Qseamsweep['Qseam_sweep_'+seam+'_'+str(mode)] = gseam/yseam
            #Cprint 'Qseam_' + seam + '_' + str(mode) + str(' = ') + str(gseam/yseam)

        if pltresult:
            _, ax = plt.subplots()
            ax.plot(values, Qseamsweep)
            ax.set_yscale('log')
            ax.set_xlabel(variable+' ('+unit+')')
            ax.set_ylabel('Q'+'_'+seam)

        return Qseamsweep

    def get_Qdielectric(self, dielectric, mode, variation):
        Qdielectric = OrderedDict()
        print('Calculating Qdielectric_' + dielectric + ' for mode ' +
              str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')

        U_dielectric = self.calc_energy_electric(variation, volume=dielectric)
        p_dielectric = U_dielectric/self.U_E
        #TODO: Update make p saved sep. and get Q for diff materials, indep. specify in pinfo
        Qdielectric['Qdielectric_'+dielectric+'_' +
                    str(mode)] = 1/(p_dielectric*config.Dissipation_params.tan_delta_sapp)
        print('p_dielectric'+'_'+dielectric+'_' +
              str(mode)+' = ' + str(p_dielectric))
        return Series(Qdielectric)

    def get_Qsurface_all(self, mode, variation):
        '''
        caculate the contribution to Q of a dieletric layer of dirt on all surfaces
        set the dirt thickness and loss tangent in the config file
        ref: http://arxiv.org/pdf/1509.01854.pdf
        '''
        lv = self.get_lv(variation)
        Qsurf = OrderedDict()
        print('Calculating Qsurface for mode ' + str(mode) +
              ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')
#        A = self.fields.Mag_E**2
#        A = A.integrate_vol(name='AllObjects')
#        U_surf = A.evaluate(lv=lv)
        calcobject = CalcObject([], self.setup)
        vecE = calcobject.getQty("E")
        A = vecE
        B = vecE.conj()
        A = A.dot(B)
        A = A.real()
        A = A.integrate_surf(name='AllObjects')
        U_surf = A.evaluate(lv=lv)
        U_surf *= config.Dissipation_params.th*epsilon_0*config.Dissipation_params.eps_r
        p_surf = U_surf/self.U_E
        Qsurf['Qsurf_'+str(mode)] = 1 / \
            (p_surf*config.Dissipation_params.tan_delta_surf)
        print('p_surf'+'_'+str(mode)+' = ' + str(p_surf))
        return Series(Qsurf)

    def calc_Q_external(self, variation, freq_GHz, U_E):
        '''
        Calculate the coupling Q of mode m with each port p
        Expected that you have specified the mode before calling this
        '''

        Qp = pd.Series({})

        freq = freq_GHz * 1e9  # freq in Hz
        for port_nm, port in self.pinfo.ports.items():
            I_peak = self.calc_avg_current_J_surf_mag(variation, port['rect'],
                                                      port['line'])
            U_dissip = 0.5 * port['R'] * I_peak**2 * 1 / freq
            p = U_dissip / (U_E/2)  # U_E is 2x the peak electrical energy
            kappa = p * freq
            Q = 2 * np.pi * freq / kappa
            Qp['Q_' + port_nm] = Q

        return Qp

    def calc_p_junction(self, variation, U_H, U_E, Ljs):
        '''
            Expected that you have specified the mode before calling this, `self.set_mode(num)`

            Expected to precalc U_H and U_E for mode, will retunr pandas series object
                junc_rect = ['junc_rect1', 'junc_rect2'] name of junc rectangles to integrate H over
                junc_len  = [0.0001]   specify in SI units; i.e., meters
                LJs       = [8e-09, 8e-09] SI units
                calc_sign = ['junc_line1', 'junc_line2']


            This function assumes there are no lumped capacitors in model.

            Potential errors:  If you dont have a line or rect by the right name you will prob get an erorr o the type:
                com_error: (-2147352567, 'Exception occurred.', (0, None, None, None, 0, -2147024365), None)
        '''

        Pj = pd.Series({})
        Sj = pd.Series({})

        for junc_nm, junc in self.pinfo.junctions.items():

            logger.debug(f'Calculating participation for {(junc_nm, junc)}')

            # Get peak current though junction I_peak
            if self.pinfo.options.method_calc_P_mj is 'J_surf_mag':
                I_peak = self.calc_avg_current_J_surf_mag(
                    variation, junc['rect'], junc['line'])

            elif self.pinfo.options.method_calc_P_mj is 'line_voltage':
                I_peak = self.calc_current_line_voltage(
                    variation, junc['line'], Ljs[junc_nm])

            else:
                raise NotImplementedError(
                    'Other calculation methods (self.pinfo.options.method_calc_P_mj) are possible but not implemented here. ')

            Pj['p_' + junc_nm] = Ljs[junc_nm] * \
                I_peak**2 / U_E
            #  divie by U_E: participation normed to be between 0 and 1 by the total capacitive energy
            #  which should be the total inductive energy

            # Sign bit
            Sj['s_' + junc_nm] = + \
                1 if (self.calc_line_current(
                    variation, junc['line'])) > 0 else -1

            if self.verbose:
                print('\t{:<15} {:>8.6g} {:>5s}'.format(
                    junc_nm,
                    Pj['p_' + junc_nm],
                    '+' if Sj['s_' + junc_nm] > 0 else '-'))

        return Pj, Sj

    def do_EPR_analysis(self,
                        variations=None,
                        modes=None):
        """
        Main analysis routine

        Load results with pyEPR_Analysis
        ..code-block python
            pyEPR_Analysis(self.data_filename, variations=variations) ```

        Optional Parameters:
        ------------------------
            variations : list | None
                Example list of variations is ['0', '1']
                A variation is a combination of project/design variables in an optimetric sweep

            modes : list | None
                Modes to analyze
                for example  modes = [0, 2, 3]

        HFSS Notes:
        ------------------------
            Assumptions:
                    Low dissipation (high-Q).
                    Right now, we assume that there are no lumped capcitors to simply calculations. Not required.
                    We assume that there are only lumped inductors, so that U_tot = U_E+U_H+U_L    and U_C =0, so that U_tot = 2*U_E;
        """

        self._run_time = time.strftime('%Y%m%d_%H%M%S', time.localtime())

        self.update_variation_information()

        if modes is None:
            modes = range(self.nmodes)

        if variations is None:
            variations = self.variations

        # Setup save and save pinfo
        #TODO:  The pd.HDFStore  is used to store the pandas sereis and dataframe, but is otherwise cumbersome.
        #       We should move to a better saving paradigm
        if self.latest_h5_path is not None and self.append_analysis:
            shutil.copyfile(self.latest_h5_path, self.data_filename)
        hdf = pd.HDFStore(self.data_filename)
        self.pinfo.save(hdf)  # This will save only 1 globalinstance

        ###  Main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for ii, variation in enumerate(variations):
            #TODO: Move this into a funciton calle self.analyze_variation

            print(f'\nVariation {variation}  [{ii+1}/{len(variations)}]')

            # Previously analyzed?
            if (variation+'/hfss_variables') in hdf.keys() and self.append_analysis:
                print_NoNewLine('  previously analyzed ...\n')
                continue

            self.lv = self.get_lv(variation)
            time.sleep(0.4)

            if self.has_fields() == False:
                logger.error(f" Error: HFSS does not have field solution for mode={ii}.\
                                Skipping this mode in the analysis")
                continue

            freqs_bare_GHz, Qs_bare = self.get_freqs_bare_pd(variation)

            self.hfss_variables[variation] = pd.Series(
                self.get_variables(variation=variation))

            Ljs = pd.Series({})
            for junc_name, val in self.junctions.items():  # junction nickname
                Ljs[junc_name] = ureg.Quantity(
                    self.hfss_variables[variation]['_'+val['Lj_variable']]).to_base_units().magnitude

            # TODO: add this as pass and then set an attribute that specifies which pass is the last pass.
            # so we can save vs pass
            hdf['v'+variation+'/freqs_bare_GHz'] = freqs_bare_GHz
            hdf['v'+variation+'/Qs_bare'] = Qs_bare
            hdf['v'+variation+'/hfss_variables'] = self.hfss_variables[variation]
            hdf['v'+variation+'/Ljs'] = Ljs

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # This is crummy now. Maybe use xarray.

            Om = OrderedDict()  # Matrix of angular frequency (of analyzed modes)
            Pm = OrderedDict()  # Participation P matrix
            Sm = OrderedDict()  # Sign          S matrix
            Qm_coupling = OrderedDict()  # Quality factor matrix
            SOL = OrderedDict()  # other results

            for mode in modes:  # integer of mode number [0,1,2,3,..]

                # Mode setup & load fields
                
                self.set_mode(mode)

                #  Get hfss solved frequencie
                _Om = pd.Series({})
                temp_freq = freqs_bare_GHz[mode]
                _Om['freq_GHz'] = temp_freq  # freq
                Om[mode] = _Om
                print(f'  Mode {mode} = {"%.2f" % temp_freq} GHz   [{mode+1}/{self.nmodes}]')
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # EPR Hamiltonian calculations

                # Calculation global energies  and report
                print('    Calculating ℰ_magnetic', end=',')
                try:
                    self.U_H = self.calc_energy_magnetic(variation)
                  #  self.U_E = self.calc_energy_electric(variation)
                except Exception as e:
                    tb = sys.exc_info()[2]
                    print("\n\nError:\n", e)
                    raise(Exception(' Did you save the field solutions?\n  Failed during calculation of the total magnetic energy. This is the first calculation step, and is indicative that there are no field solutions saved. ').with_traceback(tb))
         #       print('\n \n U_E = {} \n \n'.format(self.U_E))
            #    time.sleep(2.0)
                print('ℰ_electric')
                self.U_E = self.calc_energy_electric(variation)
                #try:
                #except Exception as e:
                 #   tb = sys.exc_info()[2]
                  #  print("\n\nError:\n", e)
                   # raise(Exception(' Did you save the field solutions?\n  Failed during calculation of the total magnetic energy. This is the first calculation step, and is indicative that there are no field solutions saved. ').with_traceback(tb))
          #      print('\n \n U_H = {} \n \n'.format(self.U_H))
                
                sol = Series({'U_H': self.U_H, 'U_E': self.U_E})

                print(f"""     {'(ℰ_E-ℰ_H)/ℰ_E':>15s} {'ℰ_E':>9s} {'ℰ_H':>9s}
    {100*(self.U_E - self.U_H)/self.U_E:>15.1f}%  {self.U_E:>9.4g} {self.U_H:>9.4g}\n""")

                # Calcualte EPR for each of the junctions
                print(f'    Calculating junction EPR [method={self.pinfo.options.method_calc_P_mj}]')
                print(f"\t{'junction':<15s} EPR p_{mode}j   sign s_{mode}j")
                Pm[mode], Sm[mode] = self.calc_p_junction(
                    variation, self.U_H, self.U_E, Ljs)

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # EPR Dissipative calculations

                # TODO: this should really be passed as argument  to the functions rather than a property of the calss I would say
                self.omega = 2*np.pi*freqs_bare_GHz[mode]

                Qm_coupling[mode] = self.calc_Q_external(variation,
                                                        freqs_bare_GHz[mode],
                                                        self.U_E)

                if self.pinfo.dissipative.seams is not None:           # get seam Q
                    for seam in self.pinfo.dissipative.seams:
                        sol = sol.append(self.get_Qseam(seam, mode, variation))

                if self.pinfo.dissipative.dielectrics_bulk is not None:     # get Q dielectric
                    for dielectric in self.pinfo.dissipative.dielectrics_bulk:
                        sol = sol.append(self.get_Qdielectric(
                            dielectric, mode, variation))

                if self.pinfo.dissipative.resistive_surfaces is not None:
                    if self.pinfo.dissipative.resistive_surfaces is 'all':             # get Q surface
                        sol = sol.append(
                            self.get_Qsurface_all(mode, variation))
                    else:
                        raise NotImplementedError(
                            "Join the team, by helping contribute this piece of code.")

                if self.pinfo.dissipative.resistive_surfaces is not None:
                    raise NotImplementedError(
                        "Join the team, by helping contribute this piece of code.")

                SOL[mode] = sol

            # Save
            self._save_variation(hdf, variation, Om, Pm, Sm, Qm_coupling, SOL)

        hdf.close()
        print('\nANALYSIS DONE. Data saved to:\n\n' + self.data_filename+'\n\n')

        return self.data_filename, variations

    def _save_variation(self, hdf, variation, Om, Pm, Sm, Qm_coupling, SOL):
        '''
        Save variation
        '''
        hdf['v'+variation+'/O_matrix'] = pd.DataFrame(Om)
        # raw, not normalized
        hdf['v'+variation+'/P_matrix'] = pd.DataFrame(Pm).transpose()
        hdf['v'+variation+'/S_matrix'] = pd.DataFrame(Sm).transpose()
        hdf['v'+variation +
            '/Q_coupling_matrix'] = pd.DataFrame(Qm_coupling).transpose()
        hdf['v'+variation+'/pyEPR_sols'] = pd.DataFrame(SOL).transpose()

        if self.options.save_mesh_stats:
            mesh = self.get_mesh_statistics(variation)
            if mesh is not None:
                hdf['v'+variation+'/mesh_stats'] = mesh  # dataframe

            conv = self.get_convergence(variation)
            if conv is not None:
                hdf['v'+variation+'/convergence'] = conv  # dataframe

            convergence_f = self.hfss_report_f_convergence(variation) # dataframe
            if (not (conv is None)) and (not (conv is [])):
                hdf['v'+variation+'/convergence_f_pass'] = convergence_f  # dataframe


    def get_mesh_statistics(self, variation='0'):
        '''
        Input:
            variation='0' ,'1', ...

        Returns dataframe:
        ```
                Name	    Num Tets	Min edge    length	    Max edge length	RMS edge length	Min tet vol	Max tet vol	Mean tet vol	Std Devn (vol)
            0	Region	    909451	    0.000243	0.860488	0.037048	    6.006260e-13	0.037352	0.000029	6.268190e-04
            1	substrate	1490356	    0.000270	0.893770	0.023639	    1.160090e-12	0.031253	0.000007	2.309920e-04
        ```
        '''
        variation = self.listvariations[ureg(variation)]
        return self.setup.get_mesh_stats(variation)

    def get_convergence(self, variation='0'):
        '''
        Input:
            variation='0' ,'1', ...

        Returns dataframe:
        ```
        	Solved Elements	Max Delta Freq. % Pass Number
            1   	    128955	        NaN
            2       	167607	        11.745000
            3       	192746	        3.208600
            4       	199244	        1.524000
        ```
        '''
        variation = self.listvariations[ureg(variation)]
        df, _ = self.setup.get_convergence(variation)
        return df

    def get_convergence_vs_pass(self, variation='0'):
        '''
        Returns a convergence vs pass number of the eignemode freqs.
        Makes a plot in HFSS that return a pandas dataframe:
            ```
        	re(Mode(1)) [g]	re(Mode(2)) [g]	re(Mode(3)) [g]
            Pass []
            1	4.643101	4.944204	5.586289
            2	5.114490	5.505828	6.242423
            3	5.278594	5.604426	6.296777
            ```
        '''
        return self.hfss_report_f_convergence(variation)

    def set_mode(self, mode_num, phase=0):
        '''
        Set source excitations should be used for fields post processing.
        Counting modes from 0 onward
        '''
        assert self.setup, "ERROR: There is no 'setup' connected."

        if mode_num < 0:
            logger.error('Too small a mode number')

        self.solutions.set_mode(mode_num + 1, phase)

        if self.has_fields() == False:
            logger.warning(f" Error: HFSS does not have field solution for mode={mode_num}.\
                                    Skipping this mode in the analysis")

        self.fields = self.setup.get_fields()

    def get_variation_nominal(self):
        return self.design.get_nominal_variation()
    def get_num_variation_nominal(self):
        try:
            return str(self.listvariations.index(self.nominalvariation))
        except:
            return '0'
    def get_variations_all(self):
        self.update_variation_information()
        return self.listvariations

    def update_variation_information(self):
        ''''
        Updates all information about the variations.
        nmodes, listvariations, nominalvariation, nvariations

        variations = ['0','1','2'] or [] for empty
        '''
        if self.setup:
            self.listvariations   = self.design._solutions.ListVariations(str(self.setup.solution_name))
            self.nominalvariation = self.design.get_nominal_variation()
            self.nvariations      = np.size(self.listvariations)
            self.variations       = [str(i) for i in range(self.nvariations)]
            self.num_variation_nominal = self.get_num_variation_nominal()
            if self.design.solution_type == 'Eigenmode':
                self.nmodes       = int(self.setup.n_modes)
            else:
                self.nmodes       = 0
    def has_fields(self, variation=None):
        '''
        Determine if fields exist for a particular solution.

        variation : str | None
        If None, gets the nominal variation
        '''
        if self.solutions:
            return self.solutions.has_fields(variation)
        else:
            False

    def hfss_report_f_convergence(self, variation= '0'):
        '''
        Create  a report in HFSS to plot the converge of freq and style it

        Returns a convergence vs pass number of the eignemode freqs.
        Returns a pandas dataframe:
            ```
        	re(Mode(1)) [g]	re(Mode(2)) [g]	re(Mode(3)) [g]
            Pass []
            1	4.643101	4.944204	5.586289
            2	5.114490	5.505828	6.242423
            3	5.278594	5.604426	6.296777
            ```
        '''
        #TODO: Move to class for reporter ?
        if not self.setup:
            logger.error('NO SETUP PRESENT - hfss_report_f_convergence.')
            return None

        if not (self.design.solution_type == 'Eigenmode'):
            return None

        oDesign = self.design
        variation = self.get_lv(variation)
        report = oDesign._reporter

        # Create report
        ycomp = [f"re(Mode({i}))" for i in range(1,1+self.nmodes)]
        params = ["Pass:=", ["All"]]+variation
        report_name = "Freq. vs. pass"
        if report_name in report.GetAllReportNames():
            report.DeleteReports([report_name])
        self.solutions.create_report(report_name, "Pass", ycomp, params, pass_name='AdaptivePass')

        # Properties of lines
        curves = [f"{report_name}:re(Mode({i})):Curve1" for i in range(1,1+self.nmodes)]
        set_property(report, 'Attributes', curves, 'Line Width', 3)
        set_property(report, 'Scaling', f"{report_name}:AxisY1", 'Auto Units', False)
        set_property(report, 'Scaling', f"{report_name}:AxisY1", 'Units', 'g')
        set_property(report, 'Legend', f"{report_name}:Legend", 'Show Solution Name', False)

        if 1: # Save
            try:
                path = Path(self.data_dir)/'hfss_eig_f_convergence.csv'
                report.ExportToFile(report_name,path)
                return pd.read_csv(path, index_col= 0)
            except Exception as e:
                logger.error(f"Error could not save and export hfss plot to {path}. Is the plot made in HFSS with the correct name. Check the HFSS error window. \t Error =  {e}")

        return None

    def hfss_report_full_convergence(self, fig=None,_display=True):

        if fig is None:
            fig = plt.figure(figsize=(11,3.))
            fig.clf()
        gs = mpl.gridspec.GridSpec(1, 3, width_ratios=[1.2, 1.5, 1])#, height_ratios=[4, 1], wspace=0.5
        axs = [fig.add_subplot(gs[i]) for i in range(3)]

        for variation in self.variations:
            print(variation)
            convergence_t = self.get_convergence()
            convergence_f = self.hfss_report_f_convergence()
            ax0t = axs[1].twinx()
            plot_convergence_f_vspass(axs[0], convergence_f)
            plot_convergence_max_df(axs[1], convergence_t.iloc[:,1])
            plot_convergence_solved_elem(ax0t, convergence_t.iloc[:,0])
            plot_convergence_maxdf_vs_sol(axs[2], convergence_t.iloc[:,1],
                                                  convergence_t.iloc[:,0])


        fig.tight_layout(w_pad=0.1)#pad=0.0, w_pad=0.1, h_pad=1.0)
        if _display:
            from IPython.display import display
            display(fig)

        return fig



#%%==============================================================================
### ANALYSIS FUNCTIONS
#==============================================================================

def pyEPR_ND(freqs, Ljs, ϕzpf,
             cos_trunc=8,
             fock_trunc=9,
             use_1st_order=False,
             return_H=False):
    '''
    Numerical diagonalizaiton for pyEPR.

    :param fs: (GHz, not radians) Linearized model, H_lin, normal mode frequencies in Hz, length M
    :param ljs: (Henries) junction linerized inductances in Henries, length J
    :param fzpfs: (reduced) Reduced Zero-point fluctutation of the junction fluxes for each mode across each junction, shape MxJ

    :return: Hamiltonian mode freq and dispersive shifts. Shifts are in MHz. Shifts have flipped sign so that down shift is positive.
    '''

    freqs, Ljs, ϕzpf = map(np.array, (freqs, Ljs, ϕzpf))
    assert(all(freqs < 1E6)), "Please input the frequencies in GHz"
    assert(all(Ljs < 1E-3)), "Please input the inductances in Henries"

    Hs = bbq_hmt(freqs * 1E9, Ljs.astype(np.float), fluxQ*ϕzpf,
                 cos_trunc, fock_trunc, individual=use_1st_order)
    f_ND, χ_ND, _, _ = make_dispersive(
        Hs, fock_trunc, ϕzpf, freqs, use_1st_order=use_1st_order)
    χ_ND = -1*χ_ND * 1E-6  # convert to MHz, and flip sign so that down shift is positive

    return (f_ND, χ_ND, Hs) if return_H else (f_ND, χ_ND)

#==============================================================================
# ANALYSIS BBQ
#==============================================================================

class Results_Hamiltonian(OrderedDict):
    
    '''
         Class to store and process results from the analysis of H_nl.
    '''
        
    def __init__(self, dict_file=None,Data_dir=None):
        """ input: 
           dict file - 1. ethier None to create an empty results hamilitoninan as 
                       as was done in the original code
                          
                       2. or a string with the name of the file where the file of the
                       previously saved Results_Hamiltonian instatnce we wish
                       to load 
                      
                       3. or an existing instance of a dict class which will be 
                       upgraded to the Results_Hamiltonian class
            Data_dir -  the directory in which the file is to be saved or loaded 
                        from, defults to the config.root_dir
                        """
        super().__init__()
        
        if Data_dir is None:
            Data_dir = Path(config.root_dir)
        
        if dict_file is None:    
            self.file_name = str(Data_dir)+'\\Results_Hamiltonian.npz'
            
        elif isinstance(dict_file,str): 
            try:
                self.file_name = str(Data_dir)+'\\' +dict_file
                self.load_from_npz()
            except:
                self.file_name = dict_file
                self.load_from_npz()
        
        elif isinstance(dict_file,dict):
            self.inject_dic(dict_file)
            self.file_name = str(Data_dir)+'\\Results_Hamiltonian.npz'
        else:
            raise ValueError('type dict_file is of type {}'.format(type(dict_file)))
                #load file
    #TODO: make this savable and loadable
    def save_to_npz(self,filename= None):
        
        if filename is None: filename = self.file_name
        
        np.savez(filename, Res_Hamil=dict(self))
        return filename
    
    def load_from_npz(self,filename= None):
        
        if filename is None: filename = self.file_name
        
        self.inject_dic(extract_dic(file_name=filename)[0])
        return filename
    
    def inject_dic(self,add_dic):
        Init_number_of_keys = len(self.keys())
        for key,val in add_dic.items():
###TODO remove all copies of same data            
          #  if key in self.keys():
                #raise ValueError('trying to overwrite an exsiting varation')
            self[str(int(key)+Init_number_of_keys)] = val
        return 1
    
    def get_vs_variation(self, quantity,lv = None,vs='variation'):
        res = OrderedDict()
        if lv is None:
            for k, r in self.items():
                res[k] = r[quantity]
        else:
   #         try:
            for key in lv:
                if vs=='variation':
                    res[key] = self[key][quantity]
                else:
                    res[str(ureg.Quantity(self[key]['hfss_variables']['_'+vs]).magnitude)] = self[key][quantity]
    #        except:
            #    print(' No such variation as ' +key)
            
        return res
    #def get_vs_variable()
    def get_frequencies_HFSS(self,lv =None,vs='variation'):
        z = sort_df_col(pd.DataFrame(self.get_vs_variation('f_0',lv =lv,vs=vs))).sort_index(axis=1)
        z.index.name   = 'eigenmode'
        z.columns.name = vs
        return z

    def get_frequencies_O1(self,lv =None,vs='variation'):
        z = sort_df_col(pd.DataFrame(self.get_vs_variation('f_1',lv =lv,vs=vs))).sort_index(axis=1)
        z.index.name   = 'eigenmode'
        z.columns.name = vs
        return z

    def get_frequencies_ND(self,lv =None,vs='variation'):
        z = sort_df_col(pd.DataFrame(self.get_vs_variation('f_ND',lv =lv,vs=vs))).sort_index(axis=1)
        z.index.name   = 'eigenmode'
        z.columns.name = vs
        return z

    def get_chi_O1(self,lv =None,vs='variation'):
        z = self.get_vs_variation('chi_O1',lv =lv,vs=vs)
        return z

    def get_chi_ND(self,lv =None,vs='variation'):
        z = self.get_vs_variation('chi_ND',lv =lv,vs=vs)
        return z


class pyEPR_Analysis(object):
    #TODO make away to propgrate to the data analysis capablites of or to merge more than a single HDF file 
    '''
        Defines an analysis object which loads and plots data from a h5 file
        This data is obtained using pyEPR_HFSS
    '''
    def __init__(self, data_filename, variations=None, do_print_info = True, Res_hamil_filename= None):

        self.data_filename = data_filename
        self.results       = Results_Hamiltonian(dict_file=Res_hamil_filename)

        with HDFStore(data_filename, mode = 'r') as hdf:  # = h5py.File(data_filename, 'r')

            self.junctions            = hdf['/project_info_junctions']
            self.project_info         = Series(hdf['project_info'])
            self.project_info_dissip  = Series(hdf['/project_info_dissip'])
            self.project_info_options = Series(hdf['/project_info_options'])

            if variations is None:
                import re
                variations = []
                for key in hdf.keys():
                    if 'hfss_variables' in key:
                        variations += [re.findall(r'\bv\d+\b', key)[0][1:]]     # result should be ['v0']  etc
                                                                                # I need to strip the v

            self.variations     = variations
            self.hfss_variables = OrderedDict()
            self.freqs_hfss     = OrderedDict()
            self.Qs             = OrderedDict()
            self.Ljs            = OrderedDict()
            self.OM             = OrderedDict()
            self.PM             = OrderedDict() # participation matrices
            self.SM             = OrderedDict() # sign matrices
            self.QM_coupling    = OrderedDict() # sign matrices
            self.sols           = OrderedDict()
            self.mesh_stats     = OrderedDict()
            self.convergence    = OrderedDict()

            for variation in self.variations:
                #try:
                    self.hfss_variables[variation] = hdf['v'+variation+'/hfss_variables']
                    self.OM[variation]             = hdf['v'+variation+'/O_matrix']
                    self.Ljs[variation]            = hdf['v'+variation+'/Ljs']
                    self.PM[variation]             = hdf['v'+variation+'/P_matrix']
                    self.SM[variation]             = hdf['v'+variation+'/S_matrix']
                    self.QM_coupling[variation]    = hdf['v'+variation+'/Q_coupling_matrix']
                    self.freqs_hfss[variation]     = hdf['v'+variation+'/freqs_bare_GHz']
                    self.Qs[variation]             = hdf['v'+variation+'/Qs_bare']
                    self.sols[variation]           = hdf['v'+variation+'/pyEPR_sols']   # contains U_E and U_H
                    self.mesh_stats[variation]     = hdf['v'+variation+'/mesh_stats']   # could be made to panel
                    self.convergence[variation]    = hdf['v'+variation+'/convergence']  # could be made to panel
                #except Exception  as e:
                #    print('\t!! ERROR in variation ' + str(variation)+ ':  ' + e)

        self.freqs_hfss           = sort_df_col(DataFrame(self.freqs_hfss))
        self.Qs                   = sort_df_col(DataFrame(self.Qs))
        self.Ljs                  = sort_df_col(DataFrame(self.Ljs))
        self.hfss_variables       = sort_df_col(DataFrame(self.hfss_variables))
        self.nmodes               = self.sols[variations[0]].shape[0]
        self._renorm_pj           = True
        dum                       = DataFrame_col_diff(self.hfss_variables)
        self.hfss_vars_diff_idx   = dum if not (dum.any() == False) else []
        try:
            self.Num_hfss_vars_diff_idx =len(self.hfss_vars_diff_idx[self.hfss_vars_diff_idx==True])
        except:
            e = sys.exc_info()[0]
            warnings.warn( "<p>Error: %s</p>" % e )
            self.Num_hfss_vars_diff_idx= 0
        if do_print_info:
            self.print_info()

    def print_info(self):
            print("\t Differences in variations:" )
            if len(self.hfss_vars_diff_idx) > 0:
                print(self.hfss_variables[self.hfss_vars_diff_idx])
            print('\n')

    def get_variable_vs(self, swpvar,lv=None):
        """ lv is list of variations (example ['0', '1']), if None it takes all variations
            swpvar is the variable by which to orginize
            
            return:
            ret -ordered dicitonary of key which is the variation number and the magnitude of swaver as the item
        """
        ret = OrderedDict()
        if lv is None:
            for key, varz in self.hfss_variables.items():
                ret[key] = ureg.Quantity(varz['_'+swpvar]).magnitude
        else:
            try:
                for key in lv:
                    ret[key] = ureg.Quantity(self.hfss_variables[key]['_'+swpvar]).magnitude
            except:
                print(' No such variation as ' +key)
        return ret
    def get_variable_value(self,swpvar,lv = None):
        
        var =self.get_variable_vs(swpvar,lv=lv)    
        return [var[key] for key in var.keys()]
    
    def get_variations_of_variable_value(self, swpvar,value,lv =None):
        """A function to return all the variations in which one of the variables has a specific value
            lv is list of variations (example ['0', '1']), if None it takes all variations
            swpvar is a string and the name of the variable we wish to filter
            value is the value of swapvr in which we are intrested
            
            returns lv - a list of the variations for which swavr==value
            """
            
        if lv is None: lv = self.variations
        
        ret = self.get_variable_vs(swpvar,lv=lv)
        
        lv = np.array(list(ret.keys()))[np.array(list(ret.values()))==value] 
        
        #lv = lv_temp if not len(lv_temp) else lv
        if  not (len(lv)):
            raise ValueError('No variations have the variable-' +swpvar +'= {}'.format(value))
            
        return lv    
            
    def get_variation_of_multiple_variables_value(self, Var_dic,lv =None):
        """
                SEE get_variations_of_variable_value
            A function to return all the variations in which one of the variables has a specific value
            lv is list of variations (example ['0', '1']), if None it takes all variations
            Var_dic is a dic with the name of the variable as key and the value to filter as item
            """
            
        if lv is None: lv = self.variations
        
        var_str = None
        for key, var in Var_dic.items():
            lv = self.get_variations_of_variable_value(key,var,lv)
            if var_str is None:
                var_str = key +'= {}'.format(var)
            else:
                var_str = var_str +' & '+ key +'= {}'.format(var)
        return lv,var_str
        
        

    def get_convergences_Max_Tets(self):
        ''' Index([u'Pass Number', u'Solved Elements', u'Max Delta Freq. %' ])  '''
        ret = OrderedDict()
        for key, df in self.convergence.items():
            ret[key] = df['Solved Elements'].iloc[-1]
        return ret

    def get_convergences_Tets_vs_pass(self):
        ''' Index([u'Pass Number', u'Solved Elements', u'Max Delta Freq. %' ])  '''
        ret = OrderedDict()
        for key, df in self.convergence.items():
            s = df['Solved Elements']
            #s.index = df['Pass Number']
            ret[key] = s
        return ret

    def get_convergences_MaxDeltaFreq_vs_pass(self):
        ''' Index([u'Pass Number', u'Solved Elements', u'Max Delta Freq. %' ])  '''
        ret = OrderedDict()
        for key, df in self.convergence.items():
            s = df['Max Delta Freq. %']
            #s.index = df['Pass Number']
            ret[key] = s
        return ret


    def get_mesh_tot(self):
        ret = OrderedDict()
        for key, m in self.mesh_stats.items():
            ret[key] = m['Num Tets  '].sum()
        return ret

    '''
    def get_solution_column(self, col_name, swp_var, sort = True):
        # sort by variation -- must be numeri
        Qs, swp = [], []
        for key, sol in self.sols.items():
            Qs  += [ sol[col_name] ]
            varz  = self.hfss_variables[key]
            swp += [ ureg.Quantity(varz['_'+swp_var]).magnitude ]
        Qs  = DataFrame(Qs, index = swp)
        return Qs if not sort else Qs.sort_index()
    '''

    def get_Qs_vs_swp(self, swp_var, sort = True):
        raise NotImplementedError()
        return self.get_solution_column('modeQ', swp_var, sort)

    def get_Fs_vs_swp(self, swp_var, sort = True):
        raise NotImplementedError()
        ''' this returns the linear frequencies that HFSS gives'''
        return self.get_solution_column('freq', swp_var, sort)

    def get_Ejs(self, variation):
        ''' EJs in GHz '''
        Ljs = self.Ljs[variation]
        Ejs = fluxQ**2/Ljs/Planck*10**-9
        return Ejs

    def analyze_all_variations(self,
                               variations =None,#None returns all_variations otherwis this is a list with number as strings ['0', '1']
                               **kwargs):
        '''
            See analyze_variation
        '''
        result = OrderedDict()
        if variations is None:  variations = self.variations
        for variation in variations:
            result[variation] = self.analyze_variation(variation,**kwargs)
        return result

    def get_Pmj(self, variation, _renorm_pj=None, print_=False):
        '''
            Get normalized Pmj Matrix

            Return DataFrame object for PJ
        '''
        if _renorm_pj is None:
            _renorm_pj = self._renorm_pj

        Pm = self.PM[variation].copy()   # EPR matrix from Jsurf avg, DataFrame

        if self._renorm_pj:  # Renormalize
            s          = self.sols[variation]
            Pm_glb_sum = (s['U_E'] - s['U_H'])/s['U_E']     # sum of participations as calculated by global UH and UE
            Pm_norm    = Pm_glb_sum/Pm.sum(axis = 1)
            # Should we still do this when Pm_glb_sum is very small
            if print_:
                print("Pm_norm = %s " % str(Pm_norm))
            Pm = Pm.mul(Pm_norm, axis=0)

        else:
            Pm_norm     = 1
            if print_:
                print('NO renorm!')

        if np.any(Pm < 0.0):
            print_color("  ! Warning:  Some p_mj was found <= 0. This is probably a numerical error, or a super low-Q mode.  We will take the abs value.  Otherwise, rerun with more precision, inspect, and do due dilligence.)")
            print(Pm,'\n')
            Pm = np.abs(Pm)

        return {'PJ':Pm, 'Pm_norm':Pm_norm}

    def get_matrices(self, variation, _renorm_pj=None, print_=False):
        '''
            All as matrices
            :PJ: Participatuion matrix, p_mj
            :SJ: Sign matrix, s_mj
            :Om: Omega_mm matrix (in GHz) (\hbar = 1) Not radians.
            :EJ: E_jj matrix of Josephson energies (in same units as hbar omega matrix)
            :PHI_zpf: ZPFs in units of \phi_0 reduced flux quantum

            Return all as *np.array*
                PM, SIGN, Om, EJ, Phi_ZPF
        '''
        #TODO: superseed by Convert.ZPF_from_EPR

        PJ = self.get_Pmj(variation, _renorm_pj=_renorm_pj, print_=print_)
        PJ = np.array(PJ['PJ'])
        SJ = np.array(self.SM[variation])                # DataFrame
        Om = np.diagflat(self.OM[variation].values)      # GHz. Frequencies of HFSS linear modes. Input in dataframe but of one line. Output nd array
        EJ = np.diagflat(self.get_Ejs(variation).values) # GHz

        print(PJ, SJ, Om, EJ)

        PHI_zpf = Calcs_basic.epr_to_zpf(PJ, SJ, Om, EJ)

        return PJ, SJ, Om, EJ, PHI_zpf                   # All as np.array


    def analyze_variation(self,
                          variation,
                          cos_trunc     = None,
                          fock_trunc    = None,
                          print_result  = True,
                          junctions     = None,
                          modes         = None):
##TODO avoide analyzing a previously analyzed variation
        '''
        Can also print results neatly.
        Args:
            junctions: list or slice of junctions to include in the analysis. None defaults to analysing all junctions
            modes: list or slice of modes to include in the analysis. None defaults to analysing all modes


        Returns
        ----------------------------
            f_0 [MHz]    : Eigenmode frequencies computed by HFSS; i.e., linear freq returned in GHz
            f_1 [MHz]    : Dressed mode frequencies (by the non-linearity; e.g., Lamb shift, etc. ). If numerical diagonalizaiton is run, then we return the numerically diagonalizaed frequencies, otherwise, use 1st order pertuirbation theory on the 4th order expansion of the cosine.
            f_1 [MHz]    : Numerical diagonalizaiton

            chi_O1 [MHz] : Analytic expression for the chis based on a cos trunc to 4th order, and using 1st order perturbation theory. Diag is anharmonicity, off diag is full cross-Kerr.
            chi_ND [MHz] : Numerically diagonalized chi matrix. Diag is anharmonicity, off diag is full cross-Kerr.
        '''
        junctions = (junctions,) if type(junctions) == int else junctions # ensuring proper matrix dimensionality when slicing
        modes = (modes,) if type(modes) == int else modes # ensuring proper matrix dimensionality when slicing


        if (fock_trunc == None) or (cos_trunc == None):
            fock_trunc = cos_trunc = None

        if print_result:
            print('\n', '. '*40)
            print('Variation %s\n' % variation)
        else:
            print('%s, ' % variation, end='')

        # Get matrices
        PJ, SJ, Om, EJ, PHI_zpf = self.get_matrices(variation)
        freqs_hfss = self.freqs_hfss[variation].values
        Ljs = self.Ljs[variation].values

        # reduce matrices to only include certain modes/junctions
        if junctions is not None:
            Ljs = Ljs[junctions,]
            PJ = PJ[:,junctions]
            SJ = SJ[:,junctions]
            EJ = EJ[:,junctions][junctions,:]
            PHI_zpf = PHI_zpf[:,junctions]
        if modes is not None:
            freqs_hfss = freqs_hfss[modes,]
            PJ = PJ[modes,:]
            SJ = SJ[modes,:]
            Om = Om[modes,:][:,modes]
            PHI_zpf = PHI_zpf[modes,:]

        # Analytic 4-th order
        CHI_O1 = 0.25* Om @ PJ @ inv(EJ) @ PJ.T @ Om * 1000. # MHz
        f1s    = np.diag(Om) - 0.5*np.ndarray.flatten( np.array(CHI_O1.sum(1))) / 1000.                  # 1st order PT expect freq to be dressed down by alpha
        CHI_O1 = divide_diagonal_by_2(CHI_O1)   # Make the diagonals alpha

        # Numerical diag
        if cos_trunc is not None:
            f1_ND, CHI_ND = pyEPR_ND(freqs_hfss,
                                     Ljs,
                                     PHI_zpf,
                                     cos_trunc     = cos_trunc,
                                     fock_trunc    = fock_trunc)
        else:
            f1_ND, CHI_ND = None, None

        result            = OrderedDict()
        result['f_0']     = self.freqs_hfss[variation]*1E3  # MHz - obtained directly from HFSS
        result['f_1']     = pd.Series(f1s)*1E3     # MHz
        result['f_ND']    = pd.Series(f1_ND)*1E-6  # MHz
        result['chi_O1']  = pd.DataFrame(CHI_O1)
        result['chi_ND']  = pd.DataFrame(CHI_ND)   # why dataframe?
        result['ZPF']     = PHI_zpf
        result['Pm_normed'] = PJ
        result['_Pm_norm']  = self.get_Pmj(variation, _renorm_pj=self._renorm_pj,
                                           print_=print_result)['Pm_norm'] # calling again
        result['hfss_variables'] = self.hfss_variables[variation] # just propagate
        result['Ljs']            = self.Ljs[variation]
        result['Q_coupling']     = self.QM_coupling[variation]
        result['Qs']             = self.Qs[variation]
        result['fock_trunc']     = fock_trunc
        result['cos_trunc']      = cos_trunc

        self.results[variation]  = result

        if print_result:
            self.print_variation(variation)
            self.print_result(result)

        return result

    def print_variation(self, variation):
        if len(self.hfss_vars_diff_idx) > 0:
            print( '\n*** Different parameters'  )
            print(self.hfss_variables[self.hfss_vars_diff_idx][variation], '\n')

        print( '*** P (participation matrix, not normlz.)'  )
        print(self.PM[variation])

        print( '\n*** S (sign-bit matrix)'  )
        print(self.SM[variation])

    def print_result(self, result):
        pritm = lambda x, frmt="{:9.2g}": print_matrix(x, frmt = frmt) #TODO: actually make into dataframe with mode labela and junction labels

        print( '*** P (participation matrix, normalized.)'  )
        pritm(result['Pm_normed'])

        print( '\n*** Chi matrix O1 PT (MHz)\n    Diag is anharmonicity, off diag is full cross-Kerr.'  )
        pritm(result['chi_O1'], "{:9.3g}")

        print( '\n*** Chi matrix ND (MHz) '  )
        pritm(result['chi_ND'], "{:9.3g}")

        print( '\n*** Frequencies O1 PT (MHz)'  )
        print(result['f_1'])

        print( '\n*** Frequencies ND (MHz)'  )
        print(result['f_ND'])

        print( '\n*** Q_coupling'  )
        print(result['Q_coupling'])
        
    def plotting_dic_x(self,Var_dic,var_name):
        dic = {}
        
        if (len(Var_dic.keys())+1)== self.Num_hfss_vars_diff_idx:
            lv,lv_str = self.get_variation_of_multiple_variables_value(Var_dic)
            dic['label']= lv_str
            dic['x_label']= var_name
            dic['x']= self.get_variable_value(var_name,lv=lv)
        else:
            raise ValueError('more than one hfss variablae changes each time')
        
        return lv, dic
    def plotting_dic_data(self,Var_dic,var_name,data_name):
        lv,dic =self.plotting_dic_x()
        dic['y_label']= data_name
        
    def plot_results(self, result,Y_label, variable,X_label ,variations =None):
        
        import matplotlib.pyplot as plt
        
        
        
    
        
    def plot_Hresults(self, variable=None, fig=None,variations =None):
        '''
            versus variations
        '''
        import matplotlib.pyplot as plt

        epr = self # lazyhack
        if not(variable is None):
            var_val = self.get_variable_value(variable,lv = variations)
        else:
            variable ='variation'
        fig, axs = plt.subplots(2,2, figsize=(10,6)) if fig is None else (fig, fig.axes)

        ax = axs[0,0]
        ax.set_title('Modal frequencies (MHz)')
        f0 = self.results.get_frequencies_HFSS(lv =variations,vs = variable)
        f1 = self.results.get_frequencies_O1(lv =variations,vs = variable)
        f_ND = self.results.get_frequencies_ND(lv =variations,vs = variable)
        mode_idx = list(f0.index)
        nmodes   = len(mode_idx)
        cmap     = cmap_discrete(nmodes)

        if f_ND.empty:
            plt_me_line = f1
            markerf1    = 'o'
        else:
            plt_me_line = f_ND
            markerf1    = '.'
            f_ND.transpose().plot(ax = ax, lw=0, marker='o',ms=4, legend=False, zorder =30, color = cmap)

        f0.transpose().plot(ax = ax, lw=0, marker='x',ms=2, legend=False, zorder = 10, color = cmap)
        f1.transpose().plot(ax = ax, lw=0, marker=markerf1,ms=4, legend=False, zorder = 20, color = cmap)
        plt_me_line.transpose().plot(ax = ax, lw=1, alpha = 0.2, color = 'grey', legend=False)

        ax = axs[1,0]
        ax.set_title('Quality factors')
        Qs = self.Qs if variations is None else self.Qs[variations] 
        Qs.transpose().plot(ax = ax, lw=0, marker=markerf1, ms=4, legend=True, zorder = 20, color = cmap)
        Qs.transpose().plot(ax = ax, lw=1, alpha = 0.2, color = 'grey', legend=False)
        ax.set_yscale('log')


        axs[0][1].set_title('Anharmonicities (MHz)')
        axs[1][1].set_title('Cross-Kerr frequencies (MHz)')
        def plot_chi_alpha(chi, primary):
            for i, m in enumerate(mode_idx):
                ax = axs[0,1]
                z = sort_Series_idx(pd.Series({k: chim.loc[m,m] for k, chim in chi.items()})).sort_index()
                z.plot(ax = ax, lw=0, ms=4, label = m, color = cmap[i], marker='o' if primary else 'x')
                if primary:
                    z.plot(ax = ax, lw=1, alpha = 0.2, color = 'grey', label = '_nolegend_')
                for i, n in enumerate(mode_idx):
                    if int(n) > int(m):
                        # plot chi
                        ax = axs[1,1]
                        z = sort_Series_idx(pd.Series({k: chim.loc[m,n] for k, chim in chi.items()}))
                        z.plot(ax = ax, lw=0, ms=4, label = str(m)+','+str(n), color = cmap[i], marker='o' if primary else 'x')
                        if primary:
                            z.plot(ax = ax, lw=1, alpha = 0.2, color = 'grey', label = '_nolegend_')

        def do_legends():
            legend_translucent(axs[0][1],  leg_kw = dict(fontsize = 7, title = 'Mode'))
            legend_translucent(axs[1][1],  leg_kw = dict(fontsize = 7))

        chiND  = self.results.get_chi_ND(lv =variations,vs = variable)
        chiO1  = self.results.get_chi_O1(lv =variations,vs = variable)
        use_ND = not (np.any([r['fock_trunc'] == None for k, r in self.results.items()]))
        if use_ND:
            plot_chi_alpha(chiND, True)
            do_legends()
            plot_chi_alpha(chiO1, False)
        else:
            plot_chi_alpha(chiO1, True)
            do_legends()

        fig.tight_layout()

        return fig, axs


def extract_dic(name=None, file_name=None):
    """#name is the name of the dictionry as saved in the npz file if it is None, 
    the function will return a list of all dictionaries in the npz file
    file name is the name of the npz file"""
    with np.load(file_name,allow_pickle=True) as f:
        if name is None:
            return  [f[i][()] for i in f.keys()]
        return [f[name][()]]