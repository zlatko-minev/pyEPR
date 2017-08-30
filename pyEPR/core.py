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
from .toolbox     import print_NoNewLine, print_color, deprecated, pi, fact, epsilon_0, hbar, Planck, fluxQ, nck, \
                         divide_diagonal_by_2, print_matrix, DataFrame_col_diff, get_instance_vars,\
                         sort_df_col, sort_Series_idx
from .toolbox_plotting import cmap_discrete, legend_translucent
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
            self.p_mj_method      = 'J_surf_mag'
            self.save_mesh_stats  = True

    def __init__(self, project_path):
        '''
        HFSS app connection settings
        -----------------------
        project_path  : str
            Directory path to the hfss project file. Should be the directory, not the file.

        junctions     : OrderedDict
            The key of this dict give the junction nickname in pyEPR.
            Each junction is given the following 4 parameters:
            * rect - Name of junction rectangles in HFSS
            * line - Name of lines in HFSS used to define the current orientation for each junction. Used to define sign of ZPF.
            * Lj_variable -Name of junction inductance variables in HFSS. DO NOT USE Global names that start with $.
            * length - of Junciton rect. length, measured in meters.
        '''
        self.project_path  = project_path
        self.project_name  = None
        self.design_name   = None
        self.setup_name    = None

        ## HFSS desgin: describe junction parameters
        self.junctions     = OrderedDict()
        # TODO: introduce modal labels

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
                  'dissipative', 'setup', '_Forbidden', 'junctions']
    def save(self, hdf):
        '''
            hdf : pd.HDFStore
        '''
        hdf['project_info']           = pd.Series(get_instance_vars(self, self._Forbidden))
        hdf['project_info_dissip']    = pd.Series(get_instance_vars(self.dissipative))
        hdf['project_info_options']   = pd.Series(get_instance_vars(self.options))
        hdf['project_info_junctions'] = pd.DataFrame(self.junctions)


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

    @property
    def junctions(self):
        return self.pinfo.junctions

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
        pj_val = (self.U_E-self.U_H)/self.U_E
        pj['pj_'+str(mode)] = np.abs(pj_val)
        print('    p_j_' + str(mode) + ' = ' + str(pj_val))
        return pj

    #TODO: replace this method with the one below, here because osme funcs use it still
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

    def get_freqs_bare_pd(self, variation):
        '''
            Retun pd.Sereis of modal freq and qs for given variation
        '''
        freqs, kappa_over_2pis = self.solutions.eigenmodes(self.get_lv_EM(variation))
        if kappa_over_2pis is None:
            kappa_over_2pis = np.zeros(len(freqs))
        freqs = pd.Series(freqs, index = range(len(freqs))) # GHz
        Qs    = freqs / pd.Series(kappa_over_2pis, index = range(len(freqs)))
        return freqs, Qs


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

    def calculate_p_mj(self, variation, U_H, U_E, Ljs):
        ''' Expected that you have specified the mode before calling this

            Expected to precalc U_H and U_E for mode, will retunr pandas series object
                junc_rect = ['junc_rect1', 'junc_rect2'] name of junc rectangles to integrate H over
                junc_len  = [0.0001]   specify in SI units; i.e., meters
                LJs       = [8e-09, 8e-09] SI units
                calc_sign = ['junc_line1', 'junc_line2']

            Potential errors:  If you dont have a line or rect by the right name you will prob get an erorr o the type:
                com_error: (-2147352567, 'Exception occurred.', (0, None, None, None, 0, -2147024365), None)
        '''

        Pj = pd.Series({})
        Sj = pd.Series({})

        for junc_nm, junc in self.pinfo.junctions.items():

            if self.pinfo.options.p_mj_method is 'J_surf_mag':
                #print(' Integrating rectangle: ' + junc['rect'])
                I_peak = self.calc_avg_current_J_surf_mag(variation, junc['rect'], junc['length'])
            else:
                raise NotImplementedError('Other calculation methods are possible but not implemented here. ')

            Pj['p_' + junc_nm] = Ljs[junc_nm] * I_peak**2 / U_E  # participation normed to 1
            Sj['s_' + junc_nm] = +1 if (self.calc_line_current(variation, junc['line'])) > 0 else -1 # sign bit

            if self.verbose:
                _str =  '(+)' if Sj['s_' + junc_nm] > 0 else '(-)'
                print('\t{:<15} {:>9.2E} '.format(junc_nm, Pj['p_' + junc_nm]) + _str)

        return Pj, Sj

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

        # Setup save and save pinfo
        #TODO:  The pd.HDFStore  is used to store the pandas sereis and dataframe, but is otherwise cumbersome.
        #       We should move to a better saving paradigm
        if self.latest_h5_path is not None and self.append_analysis:
            shutil.copyfile(self.latest_h5_path, self.data_filename)
        hdf = pd.HDFStore(self.data_filename)
        self.pinfo.save(hdf)  # This will save only 1 globalinstance

        ###  Main loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        for ii, variation in enumerate(variations):
            # Get variation, see if analyzed previously
            print('\nVariation:  ' + variation + ' / ' + str(self.nvariations-1))
            if (variation+'/hfss_variables') in hdf.keys() and self.append_analysis:
                print_NoNewLine('  previously analyzed ...\n')
                continue

            self.lv = self.get_lv(variation)
            time.sleep(0.4)
            freqs_bare_GHz, Qs_bare        = self.get_freqs_bare_pd(variation)
            self.hfss_variables[variation] = pd.Series(self.get_variables(variation=variation))
            Ljs = pd.Series({})
            for junc_name, val in self.junctions.items(): # junction nickname
                Ljs[junc_name] = ureg.Quantity(self.hfss_variables[variation]['_'+val['Lj_variable']]).to_base_units().magnitude
            hdf['v'+variation+'/freqs_bare_GHz'] = freqs_bare_GHz
            hdf['v'+variation+'/Qs_bare']        = Qs_bare
            hdf['v'+variation+'/hfss_variables'] = self.hfss_variables[variation]
            hdf['v'+variation+'/Ljs']            = Ljs

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Pm  = OrderedDict() # P matrix
            Sm  = OrderedDict() # S matrix
            SOL = OrderedDict() # otehr results
            for mode in modes:
                # Mode setup & load fields
                print('  Mode ' +  str(mode) + ' / ' + str(self.nmodes-1))
                self.solutions.set_mode(mode+1, 0)
                self.fields = self.setup.get_fields()

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # EPR calculations
                print_NoNewLine('\tU_H')
                self.U_H = self.calc_U_H(variation)
                print_NoNewLine(', U_E')
                self.U_E = self.calc_U_E(variation)
                print(  "; U_L=  {:>9.2E}" .format( (self.U_E - self.U_H )/self.U_E) )
                sol = Series({'U_H':self.U_H, 'U_E':self.U_E})
                # calcualte for each of the junctions
                Pm[mode], Sm[mode] = self.calculate_p_mj(variation, self.U_H, self.U_E, Ljs)

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Dissipative EPR calculations

                self.omega  = 2*np.pi*freqs_bare_GHz[mode]    #TODO: this should really be passed as argument  to the functions rather than a property of the calss I would say
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

                SOL[mode] = sol

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # Save
            hdf['v'+variation+'/P_matrix']   = pd.DataFrame(Pm).transpose()
            hdf['v'+variation+'/S_matrix']   = pd.DataFrame(Sm).transpose()
            hdf['v'+variation+'/pyEPR_sols'] = pd.DataFrame(SOL).transpose()

            if self.pinfo.options.save_mesh_stats:
                self._save_mesh_conv_stats(hdf, variation)


        hdf.close()
        print('\nANALYSIS DONE. Data saved to:\n\n' + self.data_filename+'\n\n')

        #pyEPR_Analysis(self.data_filename, variations=variations)
        return self.data_filename, variations

    def _save_mesh_conv_stats(self, hdf, variation):
        msh = self.setup.get_mesh_stats(self.listvariations[ureg(variation)])
        if msh is not None:
            hdf['v'+variation+'/mesh_stats']  = msh   # returns dataframe

        conv = self.setup.get_convergence(self.listvariations[ureg(variation)])  # returns dataframe
        if conv is not None:
            hdf['v'+variation+'/convergence'] = conv



#==============================================================================
### ANALYSIS FUNCTIONS
#==============================================================================


def pyEPR_ND(freqs, PJ, Om, EJ, LJs, SIGN,
             cos_trunc     = 6,
             fock_trunc    = 7,
             use_1st_order = False):
    '''
        #TODO: MAKE THE input arguments nicer?
        numerical diagonalizaiton for EPR
        fzpfs: reduced zpf  ( in units of \phi_0 )
    '''

    assert(all(freqs<1E6)), "Please input the frequencies in GHz"
    assert(all(LJs  <1E-3)),"Please input the inductances in Henries"
    assert((PJ   >0).any()),"ND -- PJs are not all > 0; \n %s" % (PJ)

    fzpfs = np.zeros(PJ.T.shape)
    for junc in range(fzpfs.shape[0]):
        for mode in range(fzpfs.shape[1]):
            fzpfs[junc, mode] = np.sqrt(PJ[mode,junc] * Om[mode,mode] /  (2*EJ[junc,junc]) ) #*0.001
    fzpfs = fzpfs * SIGN.T

    Hs = bbq_hmt(freqs*10**9, LJs.astype(np.float), fluxQ*fzpfs, cos_trunc, fock_trunc, individual = use_1st_order)
    f1s, CHI_ND, fzpfs, f0s  = make_dispersive(Hs, fock_trunc, fzpfs, freqs, use_1st_order = use_1st_order)  # f0s = freqs
    CHI_ND = -1*CHI_ND *1E-6

    return f1s, CHI_ND, fzpfs, f0s

#==============================================================================
# ANALYSIS BBQ
#==============================================================================

class Results_Hamiltonian(OrderedDict):
    '''
         Class to store and process results from the analysis of H_nl.
    '''
    #TODO: make this savable and loadable

    def get_vs_variation(self, quantity):
        res = OrderedDict()
        for k, r in self.items():
            res[k] = r[quantity]
        return res

    def get_frequencies_HFSS(self):
        z = sort_df_col(pd.DataFrame(self.get_vs_variation('f_0')))
        z.index.name   = 'eigenmode'
        z.columns.name = 'variation'
        return z

    def get_frequencies_O1(self):
        z = sort_df_col(pd.DataFrame(self.get_vs_variation('f_1')))
        z.index.name   = 'eigenmode'
        z.columns.name = 'variation'
        return z

    def get_frequencies_ND(self):
        z = sort_df_col(pd.DataFrame(self.get_vs_variation('f_ND')))
        z.index.name   = 'eigenmode'
        z.columns.name = 'variation'
        return z

    def get_chi_O1(self):
        z = self.get_vs_variation('chi_O1')
        return z

    def get_chi_ND(self):
        z = self.get_vs_variation('chi_ND')
        return z


class pyEPR_Analysis(object):
    '''
        Defines an analysis object which loads and plots data from a h5 file
        This data is obtained using pyEPR_HFSS
    '''
    def __init__(self, data_filename, variations=None, do_print_info = True):

        self.data_filename = data_filename
        self.results       = Results_Hamiltonian()

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
            self.PM             = OrderedDict() # participation matrices
            self.SM             = OrderedDict() # sign matrices
            self.sols           = OrderedDict()
            self.mesh_stats     = OrderedDict()
            self.convergence    = OrderedDict()

            for variation in self.variations:
                #try:
                    self.hfss_variables[variation] = hdf['v'+variation+'/hfss_variables']
                    self.Ljs[variation]            = hdf['v'+variation+'/Ljs']
                    self.PM[variation]             = hdf['v'+variation+'/P_matrix']
                    self.SM[variation]             = hdf['v'+variation+'/S_matrix']
                    self.freqs_hfss[variation]     = hdf['v'+variation+'/freqs_bare_GHz']
                    self.Qs[variation]             = hdf['v'+variation+'/Qs_bare']
                    self.sols[variation]           = hdf['v'+variation+'/pyEPR_sols']
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

        if do_print_info:
            self.print_info()

    def print_info(self):
            print("\t Differences in variations:" )
            if len(self.hfss_vars_diff_idx) > 0:
                print(self.hfss_variables[self.hfss_vars_diff_idx])
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
                               cos_trunc     = None,
                               fock_trunc    = None,
                               print_result  = True):
        '''
            See analyze_variation
        '''
        result = OrderedDict()
        for variation in self.variations:
            result[variation] = self.analyze_variation(variation, cos_trunc, fock_trunc, print_result)
        return result

    def analyze_variation(self,
                          variation,
                          cos_trunc     = None,
                          fock_trunc    = None,
                          print_result  = True):

        '''
        Can also print results neatly.

        Returns
        ----------------------------
            f_0 [MHz]    : Eigenmode frequencies computed by HFSS; i.e., linear freq returned in GHz
            f_1 [MHz]    : Dressed mode frequencies (by the non-linearity; e.g., Lamb shift, etc. ). If numerical diagonalizaiton is run, then we return the numerically diagonalizaed frequencies, otherwise, use 1st order pertuirbation theory on the 4th order expansion of the cosine.
            f_1 [MHz]    : Numerical diagonalizaiton

            chi_O1 [MHz] : Analytic expression for the chis based on a cos trunc to 4th order, and using 1st order perturbation theory. Diag is anharmonicity, off diag is full cross-Kerr.
            chi_ND [MHz] : Numerically diagonalized chi matrix. Diag is anharmonicity, off diag is full cross-Kerr.
        '''
        if (fock_trunc == None) or (cos_trunc == None):
            fock_trunc = cos_trunc = None

        if print_result:
            print('\n', '. '*40)
        else:
            print('%s, ' % variation, end='')

        Pm = self.PM[variation]   # EPR matrix from Jsurf avg
        if self._renorm_pj:  # Renormalize
            s          = self.sols[variation]
            Pm_glb_sum = (s['U_E'] - s['U_H'])/s['U_E']                       # sum of participations as calculated by global UH and UE
            Pm_norm    = Pm_glb_sum/Pm.sum(axis = 1)
            # should we still dothis when Pm_glb_sum is very small
            #print('\n*** NORM: ')
            print(Pm_norm)   # for debug
            Pm = Pm.mul(Pm_norm, axis=0)
        else:
            Pm_norm     = 1
            print('NO renorm!')

        if np.any(Pm < 0.0):
            print_color("  ! Warning:  Some p_mj was found <= 0. This is probably a numerical error, or a super low-Q mode.  We will take the abs value.  Otherwise, rerun with more precision, inspect, and do due dilligence.)")
            print(Pm,'\n')
            Pm = np.abs(Pm)



        # Analytic 4-th order
        f0s   = self.freqs_hfss[variation]
        PJ    = np.mat(Pm)
        Om    = np.mat(np.diagflat(f0s.values))
        EJ    = np.mat(np.diagflat(self.get_Ejs(variation).values )) # GHz
        CHI_O1= 0.25* Om * PJ * EJ.I * PJ.T * Om * 1000.             # MHz
        f1s   = f0s - 0.5*np.ndarray.flatten( np.array(CHI_O1.sum(1))) / 1000.                  # 1st order PT expect freq to be dressed down by alpha
        CHI_O1= divide_diagonal_by_2(CHI_O1)                  # Make the diagonals alpha

        # numerical diag
        if cos_trunc is not None:
            f1_ND, CHI_ND, fzpfs, f0p = pyEPR_ND(f0s, PJ, Om, EJ,
                                               self.Ljs[variation],
                                               self.SM[variation],
                                               cos_trunc     = cos_trunc,
                                               fock_trunc    = fock_trunc,
                                               use_1st_order = False)
        else:
            CHI_ND, f1_ND, fzpfs = None, None, None

        result            = OrderedDict()
        result['f_0']     = f0s*1E3                # MHz
        result['f_1']     = pd.Series(f1s)*1E3     # MHz
        result['f_ND']    = pd.Series(f1_ND)*1E-6  # MHz
        result['chi_O1']  = pd.DataFrame(CHI_O1)
        result['chi_ND']  = pd.DataFrame(CHI_ND)
        result['ZPF']     = pd.DataFrame(fzpfs)
        result['Pm_normed'] = Pm
        result['_Pm_norm']  = Pm_norm
        result['hfss_variables'] = self.hfss_variables[variation] # just propagate
        result['Ljs']            = self.Ljs[variation]
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
        print( '*** P (participation matrix, normalized.)'  )
        print(result['Pm_normed'])

        print( '\n*** Chi matrix O1 PT (MHz)\n    Diag is anharmonicity, off diag is full cross-Kerr.'  )
        print(result['chi_O1'])

        print( '\n*** Chi matrix ND (MHz) '  )
        print(result['chi_ND'])

        print( '\n*** Frequencies O1 PT (MHz)'  )
        print(result['f_1'])

        print( '\n*** Frequencies ND (MHz)'  )
        print(result['f_ND'])

    def plot_Hresults(self):
        '''
            versus varaitions
        '''
        import matplotlib.pyplot as plt

        epr = self # lazyhack

        fig, axs = plt.subplots(2,2, num=1, figsize=(10,6))

        ax = axs[0,0]
        ax.set_title('Modal frequencies (MHz)')
        f0 = epr.results.get_frequencies_HFSS()
        f1 = epr.results.get_frequencies_O1()
        f_ND = epr.results.get_frequencies_ND()
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
        Qs = epr.Qs
        Qs.transpose().plot(ax = ax, lw=0, marker=markerf1, ms=4, legend=True, zorder = 20, color = cmap)
        Qs.transpose().plot(ax = ax, lw=1, alpha = 0.2, color = 'grey', legend=False)
        ax.set_yscale('log')


        axs[0][1].set_title('Anharmonicities (MHz)')
        axs[1][1].set_title('Cross-Kerr frequencies (MHz)')
        def plot_chi_alpha(chi, primary):
            for i, m in enumerate(mode_idx):
                ax = axs[0,1]
                z = sort_Series_idx(pd.Series({k: chim.loc[m,m] for k, chim in chi.items()}))
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

        chiND  = epr.results.get_chi_ND()
        chiO1  = epr.results.get_chi_O1()
        use_ND = not (np.any([r['fock_trunc'] == None for k, r in epr.results.items()]))
        if use_ND:
            plot_chi_alpha(chiND, True)
            do_legends()
            plot_chi_alpha(chiO1, False)
        else:
            plot_chi_alpha(chiO1, True)
            do_legends()

        fig.tight_layout()

        return fig, axs






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