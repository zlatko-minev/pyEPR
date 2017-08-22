'''
 Copyright Zlatko Minev and Zaki Leghtas
 2015, 2016, 2017, 2018
'''
from __future__ import print_function    # Python 2.7 and 3 compatibility
import os
import time
import shutil
#import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Standard imports
from stat        import S_ISREG, ST_CTIME, ST_MODE
from pandas      import HDFStore, Series, DataFrame
from pint        import UnitRegistry
from collections import OrderedDict

# pyEPR custom imports
from . import hfss
from . import config
from .hfss        import CalcObject
from .toolbox     import print_NoNewLine, print_color, deprecated, pi, fact, epsilon_0, hbar, fluxQ, nck, \
                         divide_diagonal_by_2, print_matrix, DataFrame_col_diff, isint, get_instance_vars
from .numeric_diag import bbq_hmt, make_dispersive

### Definitions
ureg  = UnitRegistry(system='mks')


class Project_Info(object):
    """
    Container for HFSS project info.
    HFSS app connection settings
    -----------------------
        project_path  : str
            Directory path to the hfss project file. Should be the directory, not the file.
        project_name  : str,  None
            Name of the project within the project_path. "None" will get the current active one.
        design_name   : str,  None
            Name of the design within the project. "None" will get the current active one.
        setup_name    : str,  None
            Name of the setup within the design. "None" will get the current active one.

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
        def __init__(self):
            self.dielectrics_bulk    = None
            self.dielectric_surfaces = None
            self.resistive_surfaces  = None
            self.seams               = None

    class _pyEPR_Options:
        '''
        Pj_from_current:
            Multi-junction calculation of energy participation ratio matrix based on <I_J>. Current is integrated average of J_surf by default: (zkm 3/29/16)
            Will calculate the Pj matrix for the selected modes for the given junctions junc_rect array & length of juuncs

        pJ_method  : 'J_surf_mag'
            takes the avg. Jsurf over the rect. Make sure you have seeded lots of tets here. i recommend starting with 4 across smallest dimension.
        '''
        def __init__(self):
            self.Pj_from_current  = True
            self.pJ_method        = 'J_surf_mag'
            self.save_mesh_stats  = True

    def __init__(self, project_path):
        '''
        HFSS app connection settings
        -----------------------
        project_path  : str
            Directory path to the hfss project file. Should be the directory, not the file.
        '''
        self.project_path  = project_path
        self.project_name  = None
        self.design_name   = None
        self.setup_name    = None

        ## HFSS desgin: describe junction parameters
        self.junc_rects    = None   # Name of junction rectangles in HFSS
        self.junc_lines    = None   # Name of lines in HFSS used to define the current orientation for each junction
        self.junc_LJ_names = None   # Name of junction inductance variables in HFSS. DO NOT USE Global names that start with $.
        self.junc_lens     = None   # Junciton rect. length, measured in meters.

        ## Dissipative HFSS volumes and surfaces
        self.dissipative   = self._Dissipative()
        self.options       = self._pyEPR_Options()

        # Conected to HFSS variable
        self.app           = None
        self.desktop       = None
        self.project       = None
        self.design        = None
        self.setup         = None

    _Forbidden = ['app', 'design', 'desktop', 'project',
                  'dissipative', 'setup', '_Forbidden']
    def save(self, hdf):
        '''
            hdf : pd.HDFStore
        '''
        hdf['project_info']         = pd.Series(get_instance_vars(self, self._Forbidden))
        hdf['project_info_dissip']  = pd.Series(get_instance_vars(self.dissipative))
        hdf['project_info_options'] = pd.Series(get_instance_vars(self.options))


    def connect_to_project(self):
        '''
        Connect to HFSS design.
        '''
        print('\n\n\n')
        #print('* '*40)
        print('Connecting to HFSS ...')
        assert self.project_path is not None

        self.app, self.desktop, self.project = hfss.load_HFSS_project(self.project_name, self.project_path)
        self.design  = self.project.get_design(self.design_name) if self.design_name != None else self.project.get_active_design()
        self.setup   = self.design.get_setup(name=self.setup_name)

        self.project_name = self.project.name
        self.design_name  = self.design.name
        self.setup_name   = self.setup.name

        print('\tConnected successfully.\n\t :)\t :)\t :)\t\n')
        return self.app, self.desktop, self.project, self.design, self.setup

    def check_connected(self):
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
        assert self.check_connected() is True, "it does not appear that you have connected to HFSS yet. use connect_to_project()"
        self.project.release()
        self.desktop.release()
        self.app.release()
        hfss.release()




#==============================================================================
#%% Main compuation class & interface with HFSS
#==============================================================================
class pyEPR_HFSS(object):
    """
    This class defines a pyEPR_HFSS object which calculates and saves
    Hamiltonian parameters from an HFSS simulation.
    Further, it allows one to calcualte dissipation, etc
    """

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

    def __init__(self, project_info, verbose=True, append_analysis=False):
        '''
        Parameters:
        -------------------
            project_info    : Project_Info
                Suplpy the project info
            append_analysis : True / False
                Wheather to sppend or overwrite the save filed data.
                This feature is still in the works.

        Example use:
        -------------------
        '''
        # Input
        self.pinfo   = project_info
        if self.pinfo.check_connected() is False:
            self.pinfo.connect_to_project()
        self.verbose          = verbose
        self.append_analysis  = append_analysis

        # hfss connect module
        self.fields           = self.setup.get_fields()
        self.solutions        = self.setup.get_solutions()

        # Variations
        self.nmodes           = int(self.setup.n_modes)
        self.listvariations   = self.design._solutions.ListVariations(str(self.setup.solution_name))
        self.nominalvariation = self.design.get_nominal_variation()
        self.nvariations      = np.size(self.listvariations)

        # Solutions
        self.hfss_variables   = OrderedDict()                             # container for eBBQ list of varibles
        self.sols             = OrderedDict()                             # container for eBBQ solutions; could make a Panel
        self.mesh_stats       = OrderedDict()                             # mesh statistics for each variation
        self.conv_stats       = OrderedDict()                             # convergence statistics for each variation

        if self.verbose:
            print('Design \"%s\" info:'%self.design.name)
            print('\t%-15s %d\n\t%-15s %d' %('# eigenmodes', self.nmodes, '# variations', self.nvariations))

        self.setup_data()

        self.latest_h5_path = None # #self.get_latest_h5()
        #TODO: to be implemented
        '''
        if self.latest_h5_path is not None and self.append_analysis:
            latest_bbq_analysis = pyEPR_Analysis(self.latest_h5_path)
            if self.verbose:
                print( 'Varied variables and values : ', latest_bbq_analysis.get_swept_variables(), \
                       'Variations : ', latest_bbq_analysis.variations)
            '''

    def get_latest_h5(self):
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
            Setups up the folder path
        '''
        data_dir = config.root_dir + '/' + self.project.name + '/' + self.design.name

        #if self.verbose:
        #    print("\nResults will be saved to:\n" +'-  '*20+'\n\t'+ str(data_dir)+'\n'+'-  '*20+'\n')
        if len(self.design.name) > 50:
            print_color('WARNING!   DESING FILENAME MAY BE TOO LONG! ')

        if not os.path.isdir(data_dir):
            os.makedirs(data_dir)
        self.data_dir = data_dir
        self.data_filename = self.data_dir + '/' + self.design.name + '_' + time.strftime('%Y%m%d_%H%M%S', time.localtime()) + '.hdf5'

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

    def get_p_j(self, mode):
        '''
            This is only for a single junction!
        '''
        pj = OrderedDict()
        pj_val = (self.U_E-self.U_H)/(2*self.U_E)
        pj['pj_'+str(mode)] = np.abs(pj_val)
        print('    p_j_' + str(mode) + ' = ' + str(pj_val))
        return pj

    def get_freqs_bare(self, variation):
        #str(self.get_lv(variation))
        freqs_bare_vals = []
        freqs_bare_dict = OrderedDict()
        freqs, kappa_over_2pis = self.solutions.eigenmodes(self.get_lv_EM(variation))
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


    def get_lv(self, variation):
        ''' variation is a string #; e.g., '0'
            returns array of var names and var values '''
        if variation is None:
            lv = self.nominalvariation
            lv = self.parse_listvariations(lv)
        else:
            lv = self.listvariations[ ureg(variation) ]
            lv = self.parse_listvariations(lv)
        return lv

    def get_lv_EM(self, variation):
        if variation is None:
            lv = self.nominalvariation
            #lv = self.parse_listvariations_EM(lv)
        else:
            lv = self.listvariations[ ureg(variation) ]
            #lv = self.parse_listvariations_EM(lv)
        return str(lv)

    def parse_listvariations_EM(self,lv):
        lv = str(lv)
        lv = lv.replace("=",":=,")
        lv = lv.replace(' ',',')
        lv = lv.replace("'","")
        lv = lv.split(",")
        return lv

    def parse_listvariations(self,lv):
        lv = str(lv)
        lv = lv.replace("=",":=,")
        lv = lv.replace(' ',',')
        lv = lv.replace("'","")
        lv = lv.split(",")
        return lv

    def get_variables(self,variation=None):
        lv = self.get_lv(variation)
        variables=OrderedDict()
        for ii in range(int(len(lv)/2)):
            variables['_'+lv[2*ii][:-2]]=lv[2*ii+1]
        self.variables = variables
        return variables

    def get_Qseam(self, seam, mode, variation):
        '''
        caculate the contribution to Q of a seam, by integrating the current in
        the seam with finite conductance: set in the config file
        ref: http://arxiv.org/pdf/1509.01119.pdf
        '''
        lv = self.get_lv(variation)
        Qseam = OrderedDict()
        print('Calculating Qseam_'+ seam +' for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')
        j_2_norm = self.fields.Vector_Jsurf.norm_2() # overestimating the loss by taking norm2 of j, rather than jperp**2
        int_j_2 = j_2_norm.integrate_line(seam)
        int_j_2_val = int_j_2.evaluate(lv=lv, phase=90)
        yseam = int_j_2_val/self.U_H/self.omega
        Qseam['Qseam_'+seam+'_'+str(mode)] = config.Dissipation_params.gseam/yseam
        print('Qseam_' + seam + '_' + str(mode) + str(' = ') + str(config.Dissipation_params.gseam/config.Dissipation_params.yseam))
        return Series(Qseam)

    def get_Qseam_sweep(self, seam, mode, variation, variable, values, unit, pltresult=True):
        # values = ['5mm','6mm','7mm']
        # ref: http://arxiv.org/pdf/1509.01119.pdf

        self.solutions.set_mode(mode+1, 0)
        self.fields = self.setup.get_fields()
        freqs_bare_dict, freqs_bare_vals = self.get_freqs_bare(variation)
        self.omega = 2*np.pi*freqs_bare_vals[mode]
        print( variation)
        print( type(variation))
        print( ureg(variation))
        self.U_H = self.calc_U_H(variation)
        lv = self.get_lv(variation)
        Qseamsweep = []
        print('Calculating Qseam_'+ seam +' for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')
        for value in values:
            self.design.set_variable(variable, str(value)+unit)

            j_2_norm = self.fields.Vector_Jsurf.norm_2() # overestimating the loss by taking norm2 of j, rather than jperp**2
            int_j_2 = j_2_norm.integrate_line(seam)
            int_j_2_val = int_j_2.evaluate(lv=lv, phase=90)
            yseam = int_j_2_val/self.U_H/self.omega
            Qseamsweep.append(config.Dissipation_params.gseam/yseam)
#        Qseamsweep['Qseam_sweep_'+seam+'_'+str(mode)] = gseam/yseam
            #Cprint 'Qseam_' + seam + '_' + str(mode) + str(' = ') + str(gseam/yseam)
        if pltresult:
            fig, ax = plt.subplots()
            ax.plot(values,Qseamsweep)
            ax.set_yscale('log')
            ax.set_xlabel(variable+' ('+unit+')')
            ax.set_ylabel('Q'+'_'+seam)
        return Qseamsweep

    def get_Qdielectric(self, dielectric, mode, variation):
        Qdielectric = OrderedDict()
        print('Calculating Qdielectric_'+ dielectric +' for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')

        U_dielectric = self.calc_U_E(variation, volume=dielectric)
        p_dielectric = U_dielectric/self.U_E
        Qdielectric['Qdielectric_'+dielectric+'_'+str(mode)] = 1/(p_dielectric*config.Dissipation_params.tan_delta_sapp)
        print('p_dielectric'+'_'+dielectric+'_'+str(mode)+' = ' + str(p_dielectric))
        return Series(Qdielectric)

    def get_Qsurface_all(self, mode, variation):
        '''
        caculate the contribution to Q of a dieletric layer of dirt on all surfaces
        set the dirt thickness and loss tangent in the config file
        ref: http://arxiv.org/pdf/1509.01854.pdf
        '''
        lv = self.get_lv(variation)
        Qsurf = OrderedDict()
        print('Calculating Qsurface for mode ' + str(mode) + ' (' + str(mode) + '/' + str(self.nmodes-1) + ')')
#        A = self.fields.Mag_E**2
#        A = A.integrate_vol(name='AllObjects')
#        U_surf = A.evaluate(lv=lv)
        calcobject=CalcObject([],self.setup)
        vecE=calcobject.getQty("E")
        A=vecE
        B=vecE.conj()
        A=A.dot(B)
        A=A.real()
        A=A.integrate_surf(name='AllObjects')
        U_surf = A.evaluate(lv=lv)
        U_surf *= config.Dissipation_params.th*epsilon_0*config.Dissipation_params.eps_r
        p_surf = U_surf/self.U_E
        Qsurf['Qsurf_'+str(mode)] = 1/(p_surf*config.Dissipation_params.tan_delta_surf)
        print('p_surf'+'_'+str(mode)+' = ' + str(p_surf))
        return Series(Qsurf)

    @deprecated
    def get_Hparams(self, freqs, pjs, lj):
        '''
            Outdated.
            This is for 1 junction,  and for N=4th first-order  PT on EPR
        '''
        Hparams = OrderedDict()
        fzpfs = []

        # calculate Kerr and fzpf
        for m in self.modes:
            omega = 2*pi*freqs[m]
            ej = fluxQ**2/lj
            pj = pjs['pj_'+str(m)]
            fzpf = np.sqrt(pj*hbar*omega/ej)
            fzpfs.append(fzpf)
            Hparams['fzpf_'+str(m)] = fzpf
            alpha = 2*ej/fact(4)*nck(4,2)*(fzpf**4)/hbar
            Hparams['alpha_'+str(m)] = alpha
            Hparams['freq_'+str(m)]=(omega-alpha)/2/pi

        # calculate chi
        for m in self.modes:
            for n in self.modes: # will fail
                if n<m:
                    chi_mn = ej/hbar*(fzpfs[m]*fzpfs[n])**2
                    Hparams['chi_'+str(m)+'_'+str(n)] = chi_mn

        return Hparams

    def calc_U_E(self, variation, volume=None):
        ''' This is 2 * the peak electric energy.(since we do not divide by 2, and use the peak phasors) '''
        lv = self.get_lv(variation)
        if volume is None:
            volume = 'AllObjects'
        else:
            pass
        calcobject=CalcObject([],self.setup)
        vecE=calcobject.getQty("E")
        A=vecE.times_eps()
        B=vecE.conj()
        A=A.dot(B)
        A=A.real()
        A=A.integrate_vol(name=volume)
        return A.evaluate(lv=lv)

    def calc_U_H(self, variation, volume=None):
        lv = self.get_lv(variation)
        if volume is None:
            volume = 'AllObjects'
        else:
            pass
        calcobject=CalcObject([],self.setup)
        vecH=calcobject.getQty("H")
        A=vecH.times_mu()
        B=vecH.conj()
        A=A.dot(B)
        A=A.real()
        A=A.integrate_vol(name=volume)
        return A.evaluate(lv=lv)

    def calc_current(self, fields, line ):
        '''Function to calculate Current based on line. Not in use
            line = integration line between plates - name
        '''
        self.design.Clear_Field_Clac_Stack()
        comp = fields.Vector_H
        exp  = comp.integrate_line_tangent(line)
        I    = exp.evaluate(phase = 90)
        self.design.Clear_Field_Clac_Stack()
        return I

    def calc_avg_current_J_surf_mag(self, variation, junc_rect, junc_len):
        ''' Peak current I_max for mdoe J in junction J
            The avg. is over the surface of the junction. I.e., spatial. '''
        lv   = self.get_lv(variation)
        calc = CalcObject([],self.setup)
        calc = calc.getQty("Jsurf").mag().integrate_surf(name = junc_rect)
        I    = calc.evaluate(lv=lv) / junc_len #phase = 90
        #self.design.Clear_Field_Clac_Stack()
        return  I

    def calc_line_current(self, variation, junc_line_name):
        lv   = self.get_lv(variation)
        calc = CalcObject([],self.setup)
        calc = calc.getQty("H").imag().integrate_line_tangent(name = junc_line_name)
        #self.design.Clear_Field_Clac_Stack()
        return calc.evaluate(lv=lv)

    def calc_Pjs_from_I_for_mode(self,variation, U_H,U_E, LJs, junc_rects,junc_lens, method = 'J_surf_mag' ,
                                 freq = None, calc_sign = None):
        ''' Expected that you have specified the mode before calling this
            Expected to precalc U_H and U_E for mode, will retunr pandas series object
                junc_rect = ['junc_rect1', 'junc_rect2'] name of junc rectangles to integrate H over
                junc_len  = [0.0001]   specify in SI units; i.e., meters
                LJs       = [8e-09, 8e-09] SI units
                calc_sign = ['junc_line1', 'junc_line2']    used to define sign of ZPF
            Potential errors:  If you dont have a line or rect by the right name you will prob get an erorr o the type:
                com_error: (-2147352567, 'Exception occurred.', (0, None, None, None, 0, -2147024365), None)
        '''
        dat = OrderedDict()
        for i, junc_rect in enumerate(junc_rects):
            print_NoNewLine('     ' + junc_rect)
            if method is 'J_surf_mag':
                I_peak = self.calc_avg_current_J_surf_mag(variation, junc_rect, junc_lens[i])
            else:
                print('Not yet implemented.')
            if LJs is None:
                print_color(' -----> ERROR: Why is LJs passed as None!?')
            #dat['I_'  +junc_rect] = I_peak # stores the phase information as well
            dat['pJ_' +junc_rect] = LJs[i] * I_peak**2 / (2*U_E)
            if calc_sign is not None:
                Idum = self.calc_line_current(variation, calc_sign[i])
                dat['sign_'+junc_rect] = +1 if Idum > 0 else -1
                print(  '  %+.5f' %(dat['pJ_' +junc_rect] * dat['sign_'+junc_rect] ))
            else: print( '  %0.5f' %(dat['pJ_' +junc_rect]))
        return pd.Series(dat)

    def do_EPR_analysis(self,
                        variations       = None,
                        modes            = None):
        """
        Optional Parameters:
        ------------------------
            variations : list | none
                Example list of variations is ['0', '1']
                A variation is a combination of project/design variables in an optimetric sweep

        HFSS Notes:
        ------------------------
            Assumptions:
                    Low dissipation (high-Q).
                    Right now, we assume that there are no lumped capcitors to simply calculations. Not required.
                    We assume that there are only lumped inductors, so that U_tot = U_E+U_H+U_L    and U_C =0, so that U_tot = 2*U_E;
        """

        self._run_time = time.strftime('%Y%m%d_%H%M%S', time.localtime())

        if variations      is None:
            variations = (['-1'] if self.listvariations == (u'',)  else [str(i) for i in range(self.nvariations)] )

        if modes           is None:
            modes = range(self.nmodes)


        self._variations    = variations # for debug
        junc_rect           = self.pinfo.junc_rects
        junc_LJ_names       = self.pinfo.junc_LJ_names
        assert type(junc_LJ_names) == list, "Please pass junc_LJ_names as a list; e.g., pinfo.junc_LJ_names = ['junc1'] "
        assert type(junc_rect)     == list, "Please pass junc_rect as a list; e.g., pinfo.junc_rects = ['junc1']"


        # Setup save and save pinfo
        #TODO:  The pd.HDFStore  is used to store the pandas sereis and dataframe, but is otherwise cumbersome.
        #       We should move to a better saving paradigm
        if self.latest_h5_path is not None and self.append_analysis:
            shutil.copyfile(self.latest_h5_path, self.data_filename)
        hdf = pd.HDFStore(self.data_filename)
        self.pinfo.save(hdf)  # This will save only 1 globalinstance

        ###  Main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #TODO: Organize the results data in a better way
        #      & give better naming convntions
        for ii, variation in enumerate(variations):
            # Get variation, see if analyzed previously
            print_color('variation : ' + variation + ' / ' + str(self.nvariations-1), bg = 44)
            if (variation+'/hfss_variables') in hdf.keys() and self.append_analysis:
                print_NoNewLine('  previously analyzed ...\n')
                continue

            self.lv = self.get_lv(variation)
            time.sleep(0.4)

            freqs_bare_dict, freqs_bare_vals     = self.get_freqs_bare(variation)   # get bare freqs from HFSS
            self.hfss_variables[variation]       = pd.Series(self.get_variables(variation=variation))
            hdf['v'+variation+'/hfss_variables'] = self.hfss_variables[variation]

            self.LJs                  = [ureg.Quantity(self.hfss_variables[variation]['_'+LJvar_nm]).to_base_units().magnitude  for LJvar_nm in junc_LJ_names]
            hdf['v'+variation+'/Ljs'] = pd.Series(dict(zip(junc_LJ_names, self.LJs)))

            SOL = []
            for mode in modes:
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Mode setup & load fields

                print(' Mode  \x1b[0;30;46m ' +  str(mode) + ' \x1b[0m / ' + str(self.nmodes-1)+'  calculating:')

                sol = Series({'freq'    : freqs_bare_vals[mode]*10**-9,
                              'modeQ'   : freqs_bare_dict['Q_'+str(mode)]
                             })

                self.solutions.set_mode(mode+1, 0)
                self.fields = self.setup.get_fields()

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # EPR calculations

                print_NoNewLine('   U_H ...')
                self.U_H = self.calc_U_H(variation)
                print_NoNewLine('   U_E')
                self.U_E = self.calc_U_E(variation)
                print(  "   =>   U_L = %.3f%%" %( (self.U_E - self.U_H )/(2*self.U_E)) )
                sol['U_H'] = self.U_H
                sol['U_E'] = self.U_E

                if self.pinfo.options.Pj_from_current:
                    #print('   Calculating p_mj from fields')
                    sol_PJ = self.calc_Pjs_from_I_for_mode(variation,
                                                           self.U_H,
                                                           self.U_E,
                                                           self.LJs,
                                                           junc_rect,
                                                           self.pinfo.junc_lens,
                                                           method    = self.pinfo.options.pJ_method,
                                                           freq      = freqs_bare_vals[mode]*10**-9,
                                                           calc_sign = self.pinfo.junc_lines
                                                           )
                    sol = sol.append(sol_PJ)

                if len(junc_rect) == 1:             # Single-junction method using global U_H and U_E
                    sol['pj1'] = self.get_p_j(mode)
                    self.pjs.update(sol['pj1'])     # convinience function for single junction case,TODO: maybe this should be removed

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Dissipative EPR calculations

                self.omega  = 2*np.pi*freqs_bare_vals[mode]    #TODO: this should really be passed as argument  to the functions rather than a property of the calss I would say
                if self.pinfo.dissipative.seams is not None:           # get seam Q
                    for seam in self.pinfo.dissipative.seams:
                        sol = sol.append(self.get_Qseam(seam,mode,variation))

                if self.pinfo.dissipative.dielectrics_bulk is not None:     # get Q dielectric
                    for dielectric in self.pinfo.dissipative.dielectrics_bulk:
                        sol = sol.append(self.get_Qdielectric(dielectric, mode, variation))

                if self.pinfo.dissipative.resistive_surfaces is not None:
                    if self.pinfo.dissipative.resistive_surfaces is 'all':             # get Q surface
                        sol = sol.append( self.get_Qsurface(mode, variation) )
                    else:
                        raise NotImplementedError("Join the team, by helping contribute this piece of code.")

                if self.pinfo.dissipative.resistive_surfaces is not None:
                    raise NotImplementedError("Join the team, by helping contribute this piece of code.")

                SOL += [sol]

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Save
            if self.pinfo.options.save_mesh_stats:
                self._save_mesh_conv_stats(hdf, variation)

            self.sols[variation] = pd.DataFrame(SOL, index = modes)
            hdf['v'+variation+'/pyEPR_solution']  = self.sols[variation]


        hdf.close()
        print('ANALYSIS DONE.' + '. '*40 + '\nData saved to:\n\n' + self.data_filename+'\n\n')

        self.bbq_analysis = pyEPR_Analysis(self.data_filename, variations=variations)
        return self.bbq_analysis

    def _save_mesh_conv_stats(self, hdf, variation):
        msh = self.setup.get_mesh_stats(self.listvariations[ureg(variation)])
        if msh is not None:
            hdf['v'+variation+'/mesh_stats']  = msh   # returns dataframe

        conv = self.setup.get_convergence(self.listvariations[ureg(variation)])  # returns dataframe
        if conv is not None:
            hdf['v'+variation+'/convergence'] = conv

        self.mesh_stats[variation] = msh
        self.conv_stats[variation] = conv



#==============================================================================
### ANALYSIS FUNCTIONS
#==============================================================================


def pyEPR_ND(freqs, PJ, Om, EJ, LJs, SIGN, cos_trunc = 6, fock_trunc  = 7, use_1st_order = False):
    '''
        #TODO: MAKE THE input arguments nicer?
        numerical diagonalizaiton for energy BBQ
        fzpfs: reduced zpf  ( in units of \phi_0 )
    '''

    assert(all(freqs<1E6)), "Please input the frequencies in GHz"
    assert(all(LJs  <1E-3)),"Please input the inductances in Henries"
    assert((PJ   >0).any()),"ND -- PJs are not all > 0; \n %s" % (PJ)

    fzpfs = np.zeros(PJ.T.shape)
    for junc in range(fzpfs.shape[0]):
        for mode in range(fzpfs.shape[1]):
            fzpfs[junc, mode] = np.sqrt(PJ[mode,junc] * Om[mode,mode] /  EJ[junc,junc] ) #*0.001
    fzpfs = fzpfs * SIGN.T

    Hs = bbq_hmt(freqs*10**9, LJs.astype(np.float), fluxQ*fzpfs, cos_trunc, fock_trunc, individual = use_1st_order)
    f1s, CHI_ND, fzpfs, f0s  = make_dispersive(Hs, fock_trunc, fzpfs, freqs,use_1st_order = use_1st_order)  # f0s = freqs
    CHI_ND = -1*CHI_ND *1E-6;

    return f1s, CHI_ND, fzpfs, f0s


def pyEPR_Pmj_to_H_params(s,
                         meta_data,
                         cos_trunc     = None,
                         fock_trunc    = None,
                         _renorm_pj    = True,
                         use_1st_order = False):
    '''
    Parameters:
    ---------------
        fock_trunc    : None.   If not None, used for numerical diagonalizaiton of the Hamiltonian
        fock_trunc    : None.   If not None, used for numerical diagonalizaiton of the Hamiltonian
        _renorm_pj    : True.   If True, Use the difference in W_E and W_H to converge faster on the multi-junction pmjs
        use_1st_order : False.  If True use 1st O PT  to identify correct eigenvectors in ND

    Returns:
    ---------------
    Returns the CHIs as MHz with anharmonicity alpha as the diagonal  (with - sign)

        f0s [GHz]: Eigenmode frequencies computed by HFSS; i.e., linear freq returned in GHz
        f1s [GHz]: Dressed mode frequencies (by the non-linearity; e.g., Lamb shift, etc. ). If numerical diagonalizaiton is run, then we return the numerically diagonalizaed frequencies, otherwise, use 1st order pertuirbation theory on the 4th order expansion of the cosine.

        CHI_O1 [MHz] : Analytic expression for the chis based on a cos trunc to 4th order, and using 1st order perturbation theory.
        CHI_ND [MHz] : Numerically diagonalized chi matrix.

        PJ       : Participation matrix
        Om [GHz] : Diagnoal matrix of of linear mode (HFSS) frequencies
        EJ [GHz] : Diagonal matrix of junction energies, in GHz.

        ask Zlatko for more info.
    '''
    import  scipy
    Planck     = scipy.constants.Planck

    f0s        = np.array( s['freq'] )
    Qs         = s['modeQ']
    LJ_nms     = meta_data['junc_LJ_names']                           # ordered
    LJs        = np.array([meta_data['LJs'][nm] for nm in LJ_nms])    # LJ in Henries, must make sure these are given in the right order
    EJs        = (fluxQ**2/LJs/Planck*10**-9).astype(np.float)        # EJs in GHz
    PJ_Jsu     = s.loc[:,s.keys().str.contains('pJ')]                 # EPR from Jsurf avg
    PJ_Jsu_sum = PJ_Jsu.apply(sum, axis = 1)                          # sum of participations as calculated by avg surf current
    PJ_glb_sum = (s['U_E'] - s['U_H'])/(2*s['U_E'])                   # sum of participations as calculated by global UH and UE
    diff       = (PJ_Jsu_sum-PJ_glb_sum)/PJ_glb_sum*100               # debug

    if _renorm_pj:  # Renormalize
        PJs = PJ_Jsu.divide(PJ_Jsu_sum, axis=0).mul(PJ_glb_sum,axis=0)
    else:
        PJs = PJ_Jsu
        print('NO renorm')

    if (PJs < 0).any().any() == True:
        print("\n\n**************\n\n")
        print_color("Warning,  caution!  Some p_mj was found <= 0. This is probably a numerical error, or a super low-Q mode.  We will take the abs value.  Otherwise, rerun with more precision, inspect, and do due dilligence.)")
        print(PJs)
        print("\n\n**************\n\n")
        PJs = np.abs(PJs)

    SIGN  = s.loc[:,s.keys().str.contains('sign_')]
    PJ    = np.mat(PJs.values)
    Om    = np.mat(np.diagflat(f0s))
    EJ    = np.mat(np.diagflat(EJs))
    CHI_O1= Om * PJ * EJ.I * PJ.T * Om * 1000.      # MHz
    CHI_O1= divide_diagonal_by_2(CHI_O1)            # Make the diagonals alpha
    f1s   = f0s - np.diag(CHI_O1/1000.)             # 1st order PT expect freq to be dressed down by alpha

    if cos_trunc is not None:
        f1s, CHI_ND, fzpfs, f0s = pyEPR_ND(f0s, PJ, Om, EJ, LJs, SIGN, cos_trunc = cos_trunc, fock_trunc = fock_trunc, use_1st_order = use_1st_order)
    else:
        CHI_ND, fzpfs = None, None

    return CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs
    # the return could be made clener, or dictionary


#==============================================================================
# ANALYSIS BBQ
#==============================================================================

def sort_df_col(df):
    '''         sort by numerical int order    '''
    col_names = df.columns
    if np.all(col_names.map(isint)):
        return df[col_names.astype(int).sort_values().astype(str)]
    else:
        return df

class pyEPR_Analysis(object):
    ''' defines an analysis object which loads and plots data from a h5 file
    This data is obtained using pyEPR_HFSS

    '''
    def __init__(self, data_filename, variations=None, do_print_info = True):

        self.data_filename = data_filename
        with HDFStore(data_filename, mode = 'r') as hdf:  # = h5py.File(data_filename, 'r')

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
            self.sols           = OrderedDict()
            self.Ljs            = OrderedDict()
            self.mesh_stats     = OrderedDict()
            self.convergence    = OrderedDict()

            for variation in self.variations:
                try:
                    self.hfss_variables[variation] = hdf['v'+variation+'/hfss_variables']
                    self.Ljs[variation]            = hdf['v'+variation+'/Ljs']
                    self.sols[variation]           = hdf['v'+variation+'/pyEPR_solution']
                    self.mesh_stats[variation]     = hdf['v'+variation+'/mesh_stats']
                    self.convergence[variation]    = hdf['v'+variation+'/convergence']
                except Exception  as e:
                    print('\t!! ERROR in variation ' + str(variation)+ ':  ' + e)

        self.hfss_variables       = sort_df_col(DataFrame(self.hfss_variables))
        self.nmodes               = self.sols[variations[0]].shape[0]
        self._renorm_pj           = True

        if do_print_info:
            self.print_info()

    def print_info(self):
            print('. '*40)
            print("\t Differences in variations:" )
            print(self.hfss_variables[DataFrame_col_diff(self.hfss_variables)])
            print('\n')

    def get_variable_vs(self, swpvar):
        ret = OrderedDict()
        for key, varz in self.hfss_variables.items():
            ret[key] = ureg.Quantity(varz['_'+swpvar]).magnitude
        return ret

    def get_convergences_max_tets(self):
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

    def get_solution_column(self, col_name, swp_var, sort = True):
        ''' sort by variation -- must be numeric '''
        Qs, swp = [], []
        for key, sol in self.sols.items():
            Qs  += [ sol[col_name] ]
            varz  = self.hfss_variables[key]
            swp += [ ureg.Quantity(varz['_'+swp_var]).magnitude ]
        Qs  = DataFrame(Qs, index = swp)
        return Qs if not sort else Qs.sort_index()

    def get_Qs(self, swp_var, sort = True):
        return self.get_solution_column('modeQ', swp_var, sort)

    def get_Fs(self, swp_var, sort = True):
        ''' this returns the linear frequencies that HFSS gives'''
        return self.get_solution_column('freq', swp_var, sort)


    def get_junc_rect_names(self):
        return self.meta_data.loc['junc_rect',:]

    def analyze_variation(self,
                          variation,
                          cos_trunc     = None,
                          fock_trunc    = None,
                          print_results = True,
                          frmt          = "{:7.2f}"):
        '''
        Container function to call eBBQ_Pmj_to_H_params
        Can also print results neatly.

        Returns
        ----------------------------
        f0s [GHz]: Eigenmode frequencies computed by HFSS; i.e., linear freq returned in GHz
        f1s [GHz]: Dressed mode frequencies (by the non-linearity; e.g., Lamb shift, etc. ). If numerical diagonalizaiton is run, then we return the numerically diagonalizaed frequencies, otherwise, use 1st order pertuirbation theory on the 4th order expansion of the cosine.

        CHI_O1 [MHz] : Analytic expression for the chis based on a cos trunc to 4th order, and using 1st order perturbation theory.
        CHI_ND [MHz] : Numerically diagonalized chi matrix.

        PJ       : Participation matrix
        Om [GHz] : Diagnoal matrix of of linear mode (HFSS) frequencies
        EJ [GHz] : Diagonal matrix of junction energies, in GHz.
        '''

        print_color('. '*40, bg = 42, style = 2)
        print('\t\t Analysing single variation: %s'% variation)

        s         = self.sols[variation];
        meta_data = self.meta_data[variation]
        varz      = self.hfss_variables[variation]

        CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs = \
            pyEPR_Pmj_to_H_params(s,
                                 meta_data,
                                 cos_trunc = cos_trunc,
                                 fock_trunc = fock_trunc,
                                 _renorm_pj = self._renorm_pj)

        if print_results: ##TODO: generalize to more modes

            print( '\nPJ=\t(renorm.)'  )
            print_matrix(PJ*SIGN, frmt = "{:7.4f}")
            print( "\n","* "*5, "CHI matrix (MHz)", "* "*5)

            if cos_trunc is not None:
                print( '\nCHI_ND=\t PJ O(%d) [alpha diag]'%(cos_trunc))
                print_matrix(CHI_ND, append_row ="MHz", frmt = frmt)
            else:
                print( '\nCHI_O1=\t [alpha diag]')
                print_matrix(CHI_O1, append_row ="MHz", frmt = frmt)

            if len(f0s) == 3:
                print( '\nf0={:6.2f} {:7.2f} {:7.2f} GHz'.format(*f0s))
                print( '\nf1={:6.2f} {:7.2f} {:7.2f} GHz'.format(*(f1s*1E-9))   )
                print( 'Q={:8.1e} {:7.1e} {:6.0f}'.format(*(Qs)))
            else:
                print( "\n","* "*5, "Eigen (Linear) vs Dressed Frequencies MHz", "* "*5)
                print( pd.DataFrame(np.array([f0s*1E3,f1s*1E3]).transpose(), columns = ['Linear', 'Dressed']))
                #print( "\n", "* "*5, "Dressed freqs Frequencies MHz", "* "*5  # these are the ND if ND was used, else it is the O1PT)

                print( "\n","* "*5, "Eigen (linear) Qs ", "* "*5)
                print( pd.Series(Qs))  # Q =0 means no dissipation used in sim.

        #TODO: use dictonary in the future
        return CHI_O1, CHI_ND, PJ, Om, EJ, diff, LJs, SIGN, f0s, f1s, fzpfs, Qs, varz
    #TODO: add print feature

    '''
    @deprecated
    def get_swept_variables(self):
        #TODO: needs to be updated to new standard; currently borken
        swept_variables_names = []
        swept_variables_values = []
        for name in self.h5data[self.variations[0]].keys():
            if '_'==name[0]: # design variables all start with _
                variables = []
                for variation in self.variations:
                    variables.append(self.h5data[variation][name].value)
                if len(set(variables))>1:
                    swept_variables_names.append(name)
                    swept_variables_values.append(list(set(variables)))
            else:
                pass
        return swept_variables_names, swept_variables_values

    @deprecated
    def get_variable_variations(self, variablename):
        variables = []
        for variation in self.variations:
            variables.append(self.h5data[variation][variablename].value)
        return np.asarray(variables)

    @deprecated
    def get_float_units(self, variable_name, variation='0'):
        variable_value = self.h5data[variation][variable_name].value
        n = 1
        try:
            float(variable_value)
            return float(variable_value), ''
        except ValueError:
            while True:
                try:
                    float(variable_value[:-n])
                    return float(variable_value[:-n]), variable_value[len(variable_value)-n:]
                except:
                    n+=1
    @deprecated
    def print_Hparams(self, variation=None, modes=None):
        #TODO: needs to be updated to new standard; currently borken
        if modes==None:
            modes = range(self.nmodes)
        else:
            pass
        if variation == None:
            variation = self.variations[-1]
        else:
            pass
        swept_variables_names, swept_variables_values = self.get_swept_variables()

        for vname in swept_variables_names:
            print( vname + ' = ' + self.h5data[variation][vname].value)
        for ii, m in enumerate(modes):
            freq_m = 'freq_'+str(m)
            Kerr_m = 'alpha_'+str(m)
            Q_m = 'Q_'+str(m)
            if freq_m not in self.h5data[variation].keys():
                freq_m = 'freq_bare_'+str(m)
            else:
                pass
            if Kerr_m in self.h5data[variation].keys():
                print( Kerr_m + ' = ' +str(self.h5data[variation][Kerr_m].value/2/pi/1e6) + ' MHz')
            else:
                pass

            print( freq_m +' = ' + str(self.h5data[variation][freq_m].value/1e9) + ' GHz'   )
            if Q_m in self.h5data[variation].keys():
                print( Q_m  + ' = ' + str(self.h5data[variation][Q_m].value))
            else:
                pass

            for n in modes[0:ii]:
                chi_m_n = 'chi_'+str(m)+'_'+str(n)
                if chi_m_n in self.h5data[variation].keys():
                    print( chi_m_n + ' = ' + str(self.h5data[variation][chi_m_n].value/2/pi/1e6) + ' MHz')

    @deprecated
    def plot_Hparams(self, variable_name=None, modes=None):
        #TODO: needs to be updated to new standard; currently borken
        fig, ax = plt.subplots(2,2, figsize=(24,10))

        if variable_name == None:
            xaxis = self.variations
        else:
            xaxis = []
            for variation in self.variations:
                xaxis.append(self.get_float_units(variable_name, variation)[0])

        if modes==None:
            modes = range(self.nmodes)
        else:
            pass

        for ii, m in enumerate(modes):
            freq_m = 'freq_'+str(m)
            Kerr_m = 'alpha_'+str(m)
            Q_m = 'Q_'+str(m)
            Qsurf_m = 'Qsurf_'+str(m)

            if freq_m not in self.h5data[self.variations[0]].keys():
                freq_m = 'freq_bare_'+str(m)
            else:
                pass
            if Kerr_m in self.h5data[self.variations[0]].keys():
                ax[0][1].plot(xaxis, self.get_variable_variations(Kerr_m)/2/pi/1e6, 'o', label = str(m))
            else:
                pass

            ax[0][0].plot(xaxis, self.get_variable_variations(freq_m)/1e9, 'o', label=str(m))

            if Q_m in self.h5data[self.variations[0]].keys():
                ax[1][1].plot(xaxis, self.get_variable_variations(Q_m), 'o', label = Q_m)
            else:
                pass

            if Qsurf_m in self.h5data[self.variations[0]].keys():
                ax[1][1].plot(xaxis, self.get_variable_variations(Qsurf_m), 'o', label = Qsurf_m)
            else:
                pass

            if 'seams' in self.h5data[self.variations[0]].keys():
                for seam in self.h5data[self.variations[0]]['seams'].value:
                    Qseam_m = 'Qseam_'+seam+'_'+str(m)
                    if Qseam_m in self.h5data[self.variations[0]].keys():
                        ax[1][1].plot(xaxis, self.get_variable_variations(Qseam_m), 'o', label = Qseam_m)
                    else:
                        pass

            if 'dielectrics' in self.h5data[self.variations[0]].keys():
                for dielectric in self.h5data[self.variations[0]]['dielectrics'].value:
                    Qdielectric_m = 'Qdielectric_'+dielectric+'_'+str(m)
                    if Qdielectric_m in self.h5data[self.variations[0]].keys():
                        ax[1][1].plot(xaxis, self.get_variable_variations(Qdielectric_m), 'o', label = Qdielectric_m)
                    else:
                        pass

            for n in modes[0:ii]:
                chi_m_n = 'chi_'+str(m)+'_'+str(n)
                if chi_m_n in self.h5data[self.variations[0]].keys():
                    ax[1][0].plot(xaxis, self.get_variable_variations(chi_m_n)/2/pi/1e6, 'o', label=str(m)+','+str(n))

        ax[0][0].legend()
        ax[0][0].set_ylabel('freq (GHz)')

        ax[0][1].legend()
        ax[0][1].set_ylabel('Kerr/2pi (MHz)')
        ax[0][1].set_yscale('log')

        ax[1][0].legend()
        ax[1][0].set_ylabel('Chi/2pi (MHz)')
        ax[1][0].set_yscale('log')

        ax[1][1].legend()
        ax[1][1].set_ylabel('Q')
        ax[1][1].set_yscale('log')

        if variable_name == None:
            swept_variables_names, swept_variables_values = self.get_swept_variables()
            xticks = []
            for variation in xaxis:
                xtick = ''
                for name in swept_variables_names:
                    xtick += name[1:] + ' = ' + self.h5data[variation][name].value + '\n'
                xticks.append(str(xtick))
            ax[1][0].set_xticks([int(v) for v in xaxis])
            ax[1][0].set_xticklabels(xticks, rotation='vertical')
            ax[1][1].set_xticks([int(v) for v in xaxis])
            ax[1][1].set_xticklabels(xticks, rotation='vertical')

            ax[0][0].set_xticklabels([])
            ax[0][1].set_xticklabels([])
        else:
            xlabel = variable_name + ' (' + self.get_float_units(variable_name, self.variations[0])[1] + ')'
            ax[1][0].set_xlabel(xlabel)
            ax[1][1].set_xlabel(xlabel)

        fig.subplots_adjust(bottom=0.3)
        fig.suptitle(self.data_filename)
        fig.savefig(self.data_filename[:-5]+'.jpg')

        return fig, ax



#        for variable in swept_variables_names:
#            fig1 = plt.subplots()
#            ax1 = fig1.add_subplot(221)
#            ax.scatter()
#        return
    '''