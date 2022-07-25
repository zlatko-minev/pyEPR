"""
Main interface module to use pyEPR.

Contains code that works on the analysis after hfss, ansys, etc. These can now be closed.

Copyright Zlatko Minev, Zaki Leghtas, and the pyEPR team
2015, 2016, 2017, 2018, 2019, 2020
"""
# pylint: disable=invalid-name
# todo remove this pylint hack later

from __future__ import print_function  # Python 2.7 and 3 compatibility

from typing import List
import pickle
import sys
import time
from collections import OrderedDict
from pathlib import Path
from .calcs.convert import Convert

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import Markdown, display
from numpy.linalg import inv

# pyEPR custom imports
from . import Dict, config, logger
from .ansys import ureg
from .calcs.back_box_numeric import epr_numerical_diagonalization
from .calcs.basic import CalcsBasic
from .calcs.constants import Planck, fluxQ
from .core_distributed_analysis import DistributedAnalysis
from .toolbox.plotting import cmap_discrete, legend_translucent
from .toolbox.pythonic import (DataFrame_col_diff, divide_diagonal_by_2,
                               print_color, print_matrix, sort_df_col,
                               sort_Series_idx, df_find_index, series_of_1D_dict_to_multi_df)

from .reports import (plot_convergence_max_df, plot_convergence_solved_elem)

class HamiltonianResultsContainer(OrderedDict):
    """
    The user should only use the QuantumAnalysis class interface.

    This class is largely for internal use.

    It is a dictionary based class to contain the results stored.
    """

    file_name_extra = ' HamiltonianResultsContainer.npz'

    def __init__(self, dict_file=None, data_dir=None):
        """ input:
           dict file - 1. ethier None to create an empty results hamiltonian as
                       as was done in the original code

                       2. or a string with the name of the file where the file of the
                       previously saved HamiltonianResultsContainer instance we wish
                       to load

                       3. or an existing instance of a dict class which will be
                       upgraded to the HamiltonianResultsContainer class

            data_dir -  the directory in which the file is to be saved or loaded
                        from, defaults to the config.root_dir
        """

        super().__init__()

        self.sort_index = True  # for retrieval

        if data_dir is None:
            data_dir = Path(config.root_dir) / 'temp' / \
                time.strftime('%Y-%m-%d %H-%M-%S', time.localtime())

        data_dir = Path(data_dir).resolve()
        file_name = data_dir.stem
        directory = data_dir.parents[0]
        if not directory.is_dir():
            directory.mkdir(parents=True, exist_ok=True)

        if dict_file is None:
            self.file_name = str(
                directory/(str(file_name)+self.file_name_extra))
            #logger.info(f'Filename hamiltonian params to {self.file_name }')

        elif isinstance(dict_file, str):
            try:
                self.file_name = str(data_dir)+'\\' + dict_file
                self.load()
            except:
                self.file_name = dict_file
                self.load()

        elif isinstance(dict_file, dict):
            # Depreciated
            self._inject_dic(dict_file)
            self.file_name = str(data_dir)+self.file_name_extra

        else:
            raise ValueError(
                'type dict_file is of type {}'.format(type(dict_file)))
            # load file

    def save(self, filename: str = None):
        """
        Uses numpy npz file.
        """

        if filename is None:
            filename = self.file_name

        np.savez(filename, Res_Hamil=dict(self))
        return filename

    def load(self, filename=None):
        """
        Uses numpy npz file.
        """
        if filename is None:
            filename = self.file_name

        self._inject_dic(extract_dic(file_name=filename)[0])
        return filename

    def _inject_dic(self, add_dic):
        Init_number_of_keys = len(self.keys())
        for key, val in add_dic.items():
            # TODO remove all copies of same data
            #  if key in self.keys():
                #raise ValueError('trying to overwrite an existing variation')
            self[str(int(key)+Init_number_of_keys)] = val
        return 1

    @staticmethod
    def _do_sort_index(z: pd.DataFrame):
        """Overwrite to sort by custom function

        Arguments:
            z {pd.DataFrame} -- Input

        Returns:
            Sorted DataFrame
        """
        if isinstance(z, pd.DataFrame):
            return z.sort_index(axis=1)
        else:
            return z

    def vs_variations(self,
                      quantity: str,
                      variations: list = None,
                      vs='variation',
                      to_dataframe=False):
        """

        QUANTITIES:
            `f_0`  : HFSS Frequencies
            `f_1`  : Analytical first order PT on the p=4 term of the cosine
            `f_ND` : Numerically diagonalized
            `chi_O1`: chi matrix from 1st order PT

        Arguments:
            quantity {[type]} -- [description]

        Keyword Arguments:
            variations {list of strings} -- Variations (default: {None} -- means all)
            vs {str} -- Swept against (default: {'variation'})
            to_dataframe {bool} -- convert or not the result to dataframe.
                         Make sure to call only if it can be converted to a DataFrame or can
                         be concatenated into a multi-index DataFrame

        Returns:
            [type] -- [description]
        """
        variations = variations or self.keys()

        res = OrderedDict()
        for key in variations:
            if vs == 'variation':
                res[key] = self[key][quantity]
            else:
                # convert the key to numeric if possible
                key_new = ureg.Quantity(
                    self[key]['hfss_variables']['_'+vs]).magnitude
                res[key_new] = self[key][quantity]

        # Convert to dataframe
        z = res
        if to_dataframe:  # only call if z can be converted to a dataframe
            z = sort_df_col(pd.DataFrame(z))
            if self.sort_index:
                z = self._do_sort_index(z)
            # z.index.name = 'eigenmode'
            z.columns.name = vs

        return z

    # Quick lookup function

    def get_frequencies_HFSS(self, variations: list = None, vs='variation'):
        '''See help for `vs_variations`'''
        return self.vs_variations('f_0', variations=variations, vs=vs, to_dataframe=True)

    def get_frequencies_O1(self, variations: list = None, vs='variation'):
        '''See help for `vs_variations`'''
        return self.vs_variations('f_1', variations=variations, vs=vs, to_dataframe=True)

    def get_frequencies_ND(self, variations: list = None, vs='variation'):
        '''See help for `vs_variations`'''
        return self.vs_variations('f_ND', variations=variations, vs=vs, to_dataframe=True)

    def get_chi_O1(self, variations: list = None, vs='variation'):
        return self.vs_variations('chi_O1', variations=variations, vs=vs)

    def get_chi_ND(self, variations: list = None, vs='variation'):
        return self.vs_variations('chi_ND', variations=variations, vs=vs)


class QuantumAnalysis(object):
    '''
    Defines an analysis object which loads and plots data from a h5 file
    This data is obtained using DistributedAnalysis
    '''

    def __init__(self, data_filename,
                 variations: list = None,
                 do_print_info=True,
                 Res_hamil_filename=None):

        self.data_filename = data_filename
        self.results = HamiltonianResultsContainer(dict_file=Res_hamil_filename,
                                                   data_dir=data_filename)

        with open(str(data_filename), 'rb') as handle:
            # Contain everything: project_info and results
            self.data = Dict(pickle.load(handle))

        # Reverse from variations on outside to on inside
        results = DistributedAnalysis.results_variations_on_inside(
            self.data.results)

        # Convenience functions
        self.variations = variations or list(self.data.results.keys())
        self._hfss_variables = results['hfss_variables']
        self.freqs_hfss = results['freqs_hfss_GHz']
        self.Qs = results['Qs']
        self.Qm_coupling = results['Qm_coupling']
        self.Ljs = results['Ljs']  # DataFrame
        self.Cjs = results['Cjs']  # DataFrame
        self.OM = results['Om']  # dict of dataframes
        self.PM = results['Pm']  # participation matrices - raw, unnormed here
        # participation matrices for capacitive elements
        self.PM_cap = results['Pm_cap']
        self.SM = results['Sm']  # sign matrices
        self.I_peak = results['I_peak']
        self.V_peak = results['V_peak']
        self.modes = results['modes']

        self.sols = results['sols']
        self.ansys_energies = results.get('ansys_energies', {})

        self.mesh_stats = results['mesh']
        self.convergence = results['convergence']
        self.convergence_f_pass = results['convergence_f_pass']

        self.n_modes = len(self.modes[self.variations[0]])
        self._renorm_pj = config.epr.renorm_pj

        # Unique variation params -- make a get function
        dum = DataFrame_col_diff(self._hfss_variables)
        self.hfss_vars_diff_idx = dum if not (dum.any() == False) else []
        try:
            self.Num_hfss_vars_diff_idx = len(
                self.hfss_vars_diff_idx[self.hfss_vars_diff_idx == True])
        except:
            e = sys.exc_info()[0]
            logger.warning("<p>Error: %s</p>" % e)
            self.Num_hfss_vars_diff_idx = 0

        if do_print_info:
            self.print_info()

    @property
    def project_info(self):
        return self.data.project_info

    def print_info(self):
        print("\t Differences in variations:")
        if len(self.hfss_vars_diff_idx) > 0:
            display(self._hfss_variables[self.hfss_vars_diff_idx])
        print('\n')

    def get_vs_variable(self, swp_var, attr: str):
        """
        Convert the index of a dictionary that is stored here from
        variation number to variable value.

        Args:
            swp_var (str) :name of sweep variable in ansys
            attr: name of local attribute, eg.., 'ansys_energies'
        """
        #from collections import OrderedDict
        variable = self.get_variable_vs(swp_var)
        return OrderedDict([(variable[variation], val)
                            for variation, val in getattr(self, attr).items()])

    def get_variable_vs(self, swpvar, lv=None):
        """ lv is list of variations (example ['0', '1']), if None it takes all variations
            swpvar is the variable by which to organize

            return:
            ordered dictionary of key which is the variation number and the magnitude
            of swaver as the item
        """
        ret = OrderedDict()
        if lv is None:
            for key, varz in self._hfss_variables.items():
                ret[key] = ureg.Quantity(varz['_'+swpvar]).magnitude
        else:
            try:
                for key in lv:
                    ret[key] = ureg.Quantity(
                        self._hfss_variables[key]['_'+swpvar]).magnitude
            except:
                print(' No such variation as ' + key)
        return ret

    def get_variable_value(self, swpvar, lv=None):

        var = self.get_variable_vs(swpvar, lv=lv)
        return [var[key] for key in var.keys()]

    def get_variations_of_variable_value(self, swpvar, value, lv=None):
        """A function to return all the variations in which one of the variables
            has a specific value lv is list of variations (example ['0', '1']),
            if None it takes all variations
            swpvar is a string and the name of the variable we wish to filter
            value is the value of swapvr in which we are interested

            returns lv - a list of the variations for which swavr==value
            """

        if lv is None:
            lv = self.variations

        ret = self.get_variable_vs(swpvar, lv=lv)

        lv = np.array(list(ret.keys()))[np.array(list(ret.values())) == value]

        #lv = lv_temp if not len(lv_temp) else lv
        if not (len(lv)):
            raise ValueError('No variations have the variable-' + swpvar +
                             '= {}'.format(value))

        return list(lv)

    def get_variation_of_multiple_variables_value(self, Var_dic, lv=None):
        """
                SEE get_variations_of_variable_value
            A function to return all the variations in which one of the variables has a specific value
            lv is list of variations (example ['0', '1']), if None it takes all variations
            Var_dic is a dic with the name of the variable as key and the value to filter as item
            """

        if lv is None:
            lv = self.variations

        var_str = None
        for key, var in Var_dic.items():
            lv = self.get_variations_of_variable_value(key, var, lv)
            if var_str is None:
                var_str = key + '= {}'.format(var)
            else:
                var_str = var_str + ' & ' + key + '= {}'.format(var)
        return lv, var_str

    def get_convergences_max_tets(self):
        ''' Index([u'Pass Number', u'Solved Elements', u'Max Delta Freq. %' ])  '''
        ret = OrderedDict()
        for key, df in self.convergence.items():
            ret[key] = df['Solved Elements'].iloc[-1]
        return ret

    def get_convergences_tets_vs_pass(self, as_dataframe=True):
        ''' Index([u'Pass Number', u'Solved Elements', u'Max Delta Freq. %' ])  '''
        ret = OrderedDict()
        for key, df in self.convergence.items():
            s = df['Solved Elements']
            s = s.reset_index().dropna().set_index('Pass Number')
            #s.index = df['Pass Number']
            ret[key] = s

        if as_dataframe:
            ret = pd.concat(ret)
            ret = ret.unstack(0)['Solved Elements']

        return ret

    def get_convergences_max_delta_freq_vs_pass(self, as_dataframe=True):
        ''' Index([u'Pass Number', u'Solved Elements', u'Max Delta Freq. %' ])  '''
        KEY = 'Max Delta Freq. %'

        ret = OrderedDict()
        for key, df in self.convergence.items():
            s = df[KEY]
            s = s.reset_index().dropna().set_index('Pass Number')
            #s.index = df['Pass Number']
            ret[key] = s

        if as_dataframe:
            ret = pd.concat(ret)
            ret = ret.unstack(0)[KEY]

        return ret

    def get_mesh_tot(self):
        ret = OrderedDict()
        for key, m in self.mesh_stats.items():
            ret[key] = m['Num Tets  '].sum()
        return ret

    def get_Ejs(self, variation):
        ''' EJs in GHz
        See calcs.convert
        '''
        Ljs = self.Ljs[variation]
        Ejs = fluxQ**2/Ljs/Planck*10**-9
        return Ejs

    def get_Ecs(self, variation):
        ''' ECs in GHz
        Returns as pandas series
        '''
        Cs = self.Cjs[variation]
        return Convert.Ec_from_Cs(Cs,  units_in='F', units_out='GHz')

    def analyze_all_variations(self,
                               variations: List[str] = None,
                               analyze_previous=False,
                               **kwargs):
        '''
        See analyze_variation for full documentation

        Args:
            variations: None returns all_variations otherwise this is a list with number as strings ['0', '1']
            analyze_previous: set to true if you wish to overwrite previous analysis
            **kwargs: Keyword arguments passed to :func:`~pyEPR.QuantumAnalysis.analyze_variation`.
        '''

        result = OrderedDict()

        if variations is None:
            variations = self.variations

        for variation in variations:
            if (not analyze_previous) and (variation in self.results.keys()):
                result[variation] = self.results[variation]
            else:
                result[variation] = self.analyze_variation(variation, **kwargs)


        self.results.save()

        return result

    def _get_ansys_total_energies(self, variation):
        res = {}
        for getkey in ['U_tot_cap', 'U_tot_ind', 'U_H', 'U_E', 'U_norm']:
            res[getkey] = pd.Series({mode: self.ansys_energies[variation][mode][getkey]
                                     for mode in self.ansys_energies[variation]})
        df = pd.DataFrame(res)
        df.index.name = 'modes'
        return df

    def _get_participation_normalized(self, variation, _renorm_pj=None, print_=False):
        '''
            Get normalized Pmj Matrix

            Return DataFrame object for PJ
        '''
        if _renorm_pj is None:
            _renorm_pj = self._renorm_pj

        # Columns are junctions; rows are modes
        Pm = self.PM[variation].copy()   # EPR matrix DataFrame

        # EPR matrix for capacitor DataFrame
        Pm_cap = self.PM_cap[variation].copy()

        if _renorm_pj:  # just non False

            # Renormalize
            # Should we still do this when Pm_glb_sum is very small
            #s = self.sols[variation]
            # sum of participation energies as calculated by global UH and UE
            # U_mode = s['U_E'] # peak mode energy; or U bar as i denote it sometimes
            # We need to add the capacitor here, and maybe take the mean of that

            energies = self._get_ansys_total_energies(variation)

            U_mode = (energies['U_tot_cap'] + energies['U_tot_ind'])/2.
            U_diff = abs(energies['U_tot_cap'] - energies['U_tot_ind'])/U_mode
            if np.any(U_diff > 0.15):
                logger.error(f"WARNING: U_tot_cap-U_tot_ind / mean = {np.max(np.abs(U_diff))*100:.1f}% is > 15%. \
                    \nIs the simulation converged? Proceed with caution")

            # global sums of participations
            Pm_glb_sum = abs((U_mode-energies['U_H'])/U_mode)
            Pm_cap_glb_sum = abs((U_mode-energies['U_E'])/U_mode)

            # norms
            Pm_norm = Pm_glb_sum/Pm.sum(axis=1)
            Pm_cap_norm = Pm_cap_glb_sum/Pm_cap.sum(axis=1)
            # this is not the correct scaling yet! WARNING. Factors of 2 laying around too
            # these numbers are a bit all over the place for now. very small

            if _renorm_pj == True or _renorm_pj == 1:
                idx = Pm > -1E6  # everywhere scale
                idx_cap = Pm_cap > -1E6
            elif _renorm_pj == 2:
                idx = Pm > 0.15  # Mask for where to scale
                idx_cap = Pm_cap > 0.15
            else:
                raise NotImplementedError(
                    "Unknown _renorm_pj argument or config values!")

            if print_:
                # \nPm_cap_norm=\n{Pm_cap_norm}")
                print(f"Pm_norm=\n{Pm_norm}\n")
                print(f"Pm_norm idx =\n{idx}")


            Pm[idx] = Pm[idx].mul(Pm_norm, axis=0)
            Pm_cap[idx_cap] = Pm_cap[idx_cap].mul(Pm_cap_norm, axis=0)
            #Pm = Pm.mul(Pm_norm, axis=0)
            #Pm_cap = Pm_cap.mul(Pm_cap_norm, axis=0)

        else:
            Pm_norm = 1
            Pm_cap_norm = 1
            idx = None
            idx_cap = None
            if print_:
                print('NO renorm!')

        if np.any(Pm < 0.0):
            print_color("  ! Warning:  Some p_mj was found <= 0. This is probably a numerical error,'\
                'or a super low-Q mode.  We will take the abs value.  Otherwise, rerun with more precision,'\
                'inspect, and do due diligence.)")
            print(Pm, '\n')
            Pm = np.abs(Pm)

        return {'PJ': Pm, 'Pm_norm': Pm_norm, 'PJ_cap': Pm_cap,
                'Pm_cap_norm': Pm_cap_norm,
                'idx': idx,
                'idx_cap': idx_cap}

    def get_epr_base_matrices(self, variation, _renorm_pj=None, print_=False):
        r'''
        Return the key matrices used in the EPR method for analytic calculations.

        All as matrices
            :PJ: Participation matrix, p_mj
            :SJ: Sign matrix, s_mj
            :Om: Omega_mm matrix (in GHz) (\hbar = 1) Not radians.
            :EJ: E_jj matrix of Josephson energies (in same units as hbar omega matrix)
            :PHI_zpf: ZPFs in units of \phi_0 reduced flux quantum
            :PJ_cap: capacitive participation matrix

            Return all as *np.array*
                PM, SIGN, Om, EJ, Phi_ZPF
        '''
        # TODO: supersede by Convert.ZPF_from_EPR

        res = self._get_participation_normalized(
            variation, _renorm_pj=_renorm_pj, print_=print_)

        PJ = np.array(res['PJ'])
        PJ_cap = np.array(res['PJ_cap'])

        # Sign bits
        SJ = np.array(self.SM[variation])                # DataFrame
        #  Frequencies of HFSS linear modes.
        #  Input in dataframe but of one line. Output nd array
        Om = np.diagflat(self.OM[variation].values)      # GHz
        # Junction energies
        EJ = np.diagflat(self.get_Ejs(variation).values)  # GHz
        Ec = np.diagflat(self.get_Ecs(variation).values)  # GHz

        for x in ("PJ", "SJ", "Om", "EJ"):
            logger.debug(f"{x}=")
            logger.debug(locals()[x])

        PHI_zpf = CalcsBasic.epr_to_zpf(PJ, SJ, Om, EJ)
        n_zpf = CalcsBasic.epr_cap_to_nzpf(PJ, SJ, Om, Ec)

        return PJ, SJ, Om, EJ, PHI_zpf, PJ_cap, n_zpf                 # All as np.array

    def analyze_variation(self,
                          variation: str,
                          cos_trunc: int = None,
                          fock_trunc: int = None,
                          print_result: bool = True,
                          junctions: List = None,
                          modes: List = None):
        # TODO avoid analyzing a previously analyzed variation
        '''
        Core analysis function to call!

        Args:
            junctions: list or slice of junctions to include in the analysis.
                None defaults to analysing all junctions
            modes: list or slice of modes to include in the analysis.
                None defaults to analysing all modes

        Returns:
            dict: Dictionary containing at least the following:
                * f_0 [MHz]: Eigenmode frequencies computed by HFSS; i.e., linear freq returned in GHz
                * f_1 [MHz]: Dressed mode frequencies (by the non-linearity; e.g., Lamb shift, etc. ).
                  Result based on 1st order perturbation theory on the 4th order expansion of the cosine.
                * f_ND [MHz]: Numerical diagonalization result of dressed mode frequencies.
                  only available if `cos_trunc` and  `fock_trunc` are set (non None).
                * chi_O1 [MHz]: Analytic expression for the chis based on a cos trunc to 4th order, and using 1st
                  order perturbation theory. Diag is anharmonicity, off diag is full cross-Kerr.
                * chi_ND [MHz]: Numerically diagonalized chi matrix. Diag is anharmonicity, off diag is full
                  cross-Kerr.
        '''

        # ensuring proper matrix dimensionality when slicing
        junctions = (junctions,) if type(junctions) is int else junctions

        if modes is None:
            modes = list(range(self.n_modes))
        
        tmp_n_modes = self.n_modes
        tmp_modes = self.modes[variation]
        self.n_modes = len(modes)
        self.modes[variation] = modes

        if (fock_trunc is None) or (cos_trunc is None):
            fock_trunc = cos_trunc = None

        if print_result:
            print('\n', '. '*40)
            print('Variation %s\n' % variation)
        else:
            print('%s, ' % variation, end='')

        # Get matrices
        PJ, SJ, Om, EJ, PHI_zpf, PJ_cap, n_zpf = self.get_epr_base_matrices(
            variation)
        freqs_hfss = self.freqs_hfss[variation].values[(modes)]
        Ljs = self.Ljs[variation].values

        # reduce matrices to only include certain modes/junctions
        if junctions is not None:
            Ljs = Ljs[junctions, ]
            PJ = PJ[:, junctions]
            SJ = SJ[:, junctions]
            EJ = EJ[:, junctions][junctions, :]
            PHI_zpf = PHI_zpf[:, junctions]
            PJ_cap = PJ_cap[:, junctions]

        if modes is not None:
            freqs_hfss = freqs_hfss[range(len(self.modes[variation])), ]
            PJ = PJ[range(len(modes)), :]
            SJ = SJ[range(len(modes)), :]
            Om = Om[range(len(modes)), :][:, range(len(modes))]
            PHI_zpf = PHI_zpf[range(len(modes)), :]
            PJ_cap = PJ_cap[:, junctions]

        # Analytic 4-th order
        CHI_O1 = 0.25 * Om @ PJ @ inv(EJ) @ PJ.T @ Om * 1000.  # MHz
        f1s = np.diag(Om) - 0.5*np.ndarray.flatten(np.array(CHI_O1.sum(1))) / \
            1000.                  # 1st order PT expect freq to be dressed down by alpha
        CHI_O1 = divide_diagonal_by_2(CHI_O1)   # Make the diagonals alpha

        # Numerical diag
        if cos_trunc is not None:
            f1_ND, CHI_ND = epr_numerical_diagonalization(freqs_hfss,
                                                          Ljs,
                                                          PHI_zpf,
                                                          cos_trunc=cos_trunc,
                                                          fock_trunc=fock_trunc)
        else:
            f1_ND, CHI_ND = None, None

        result = OrderedDict()
        result['f_0'] = self.freqs_hfss[variation][modes] * 1E3  # MHz - obtained directly from HFSS
        result['f_1'] = pd.Series(f1s)*1E3     # MHz
        result['f_ND'] = pd.Series(f1_ND)*1E-6  # MHz
        result['chi_O1'] = pd.DataFrame(CHI_O1)
        result['chi_ND'] = pd.DataFrame(CHI_ND)   # why dataframe?
        result['ZPF'] = PHI_zpf
        result['Pm_normed'] = PJ
        try:
            result['Pm_raw'] = self.PM[variation][self.PM[variation].columns[0]][modes]#TODO change the columns to junctions
        except:
             result['Pm_raw'] = self.PM[variation]
        _temp = self._get_participation_normalized(
            variation, _renorm_pj=self._renorm_pj, print_=print_result)
        result['_Pm_norm'] = _temp['Pm_norm'][modes]
        result['_Pm_cap_norm'] = _temp['Pm_cap_norm'][modes]

        # just propagate
        result['hfss_variables'] = self._hfss_variables[variation]
        result['Ljs'] = self.Ljs[variation]
        result['Cjs'] = self.Cjs[variation]
        try:
            result['Q_coupling'] = self.Qm_coupling[variation][self.Qm_coupling[variation].columns[junctions]][modes]#TODO change the columns to junctions
        except:
             result['Q_coupling'] = self.Qm_coupling[variation]
        
        try:
            result['Qs'] = self.Qs[variation][self.PM[variation].columns[junctions]][modes] #TODO change the columns to junctions
        except:
             result['Qs'] = self.Qs[variation][modes]
        result['fock_trunc'] = fock_trunc
        result['cos_trunc'] = cos_trunc

        self.results[variation] = result
        self.results.save()
        
        if print_result:
            self.print_variation(variation)
            self.print_result(result)
    
        self.n_modes = tmp_n_modes # TODO is this smart should consider defining the modes of interest in the initialisation of the quantum object
        self.modes[variation]=tmp_modes 
        return result

    def full_report_variations(self, var_list: list=None):
        """see full_variation_report"""
        if var_list is None: var_list =self.variations
        for variation in var_list: 
            self.full_variation_report(variation)
    
    def full_variation_report(self, variation):
        """
        prints the results and parameters of a specific variation

        Parameters
        ----------
        variation : int or str
            the variation to be printed .

        Returns
        -------
        None.

        """
        self.print_variation(variation)
        
        self.print_result(variation)

    def print_variation(self, variation):
        """
        Utility reporting function
        """
        if variation is int: variation = str(variation)

        if len(self.hfss_vars_diff_idx) > 0:
            print('\n*** Different parameters')
            display(self._hfss_variables[self.hfss_vars_diff_idx][variation])
            print('\n')

        print('*** P (participation matrix, not normlz.)')
        print(self.PM[variation])

        print('\n*** S (sign-bit matrix)')
        print(self.SM[variation])

    def print_result(self, result):
        """
        Utility reporting function
        """
        if type(result) is str or type(result) is int: result = self.results[str(result)]

        # TODO: actually make into dataframe with mode labels and junction labels
        pritm = lambda x, frmt="{:9.2g}": print_matrix(x, frmt=frmt)

        print('*** P (participation matrix, normalized.)')
        pritm(result['Pm_normed'])

        print('\n*** Chi matrix O1 PT (MHz)\n    Diag is anharmonicity, off diag is full cross-Kerr.')
        pritm(result['chi_O1'], "{:9.3g}")

        print('\n*** Chi matrix ND (MHz) ')
        pritm(result['chi_ND'], "{:9.3g}")

        print('\n*** Frequencies O1 PT (MHz)')
        print(result['f_1'])

        print('\n*** Frequencies ND (MHz)')
        print(result['f_ND'])

        print('\n*** Q_coupling')
        print(result['Q_coupling'])

    def plotting_dic_x(self, Var_dic, var_name):
        dic = {}

        if (len(Var_dic.keys())+1) == self.Num_hfss_vars_diff_idx:
            lv, lv_str = self.get_variation_of_multiple_variables_value(
                Var_dic)
            dic['label'] = lv_str
            dic['x_label'] = var_name
            dic['x'] = self.get_variable_value(var_name, lv=lv)
        else:
            raise ValueError('more than one hfss variable changes each time')

        return lv, dic

    # Does not seem used. What is Var_dic and var_name going to?
    # def plotting_dic_data(self, Var_dic, var_name, data_name):
    #   lv, dic = self.plotting_dic_x()
    #    dic['y_label'] = data_name

    def plot_results(self, result, Y_label, variable, X_label, variations: list = None):
        # TODO?
        pass

    def plot_hamiltonian_results(self,
                                 swp_variable: str = 'variation',
                                 variations: list = None,
                                 fig=None,
                                 x_label: str = None):
        """Plot results versus variation

        Keyword Arguments:
            swp_variable {str} -- Variable against which we swept. If none, then just
                                    take the variation index (default: {None})
            variations {list} -- [description] (default: {None})
            fig {[type]} -- [description] (default: {None})

        Returns:
            fig, axs
        """
        x_label = x_label or swp_variable

        # Create figure and axes
        if not fig:
            fig, axs = plt.subplots(2, 2, figsize=(10, 6))
        else:
            axs = fig.axs

        ############################################################################
        # Axis: Frequencies
        f0 = self.results.get_frequencies_HFSS(
            variations=variations, vs=swp_variable).transpose().sort_index(key=lambda x : x.astype(int))
        f1 = self.results.get_frequencies_O1(
            variations=variations, vs=swp_variable).transpose().sort_index(key=lambda x : x.astype(int))
        f_ND = self.results.get_frequencies_ND(
            variations=variations, vs=swp_variable).transpose().sort_index(key=lambda x : x.astype(int))
        # changed by Asaf from f0 as not all modes are always analyzed
        mode_idx = list(f1.columns)
        n_modes = len(mode_idx)

        ax = axs[0, 0]
        ax.set_title('Modal frequencies (MHz)')

        # TODO: should move these kwargs to the config
        cmap = cmap_discrete(n_modes)
        kw = dict(ax=ax, color=cmap, legend=False, lw=0, ms=0)

        # Choose which freq should have the solid line drawn with it. ND if present, else f1
        if f_ND.empty:
            plt_me_line = f1
            markerf1 = 'o'
        else:
            plt_me_line = f_ND
            markerf1 = '.'

            # plot the ND as points if present
            f_ND.plot(**{**kw, **dict(marker='o', ms=4, zorder=30)})

        f0.plot(**{**kw, **dict(marker='x', ms=2, zorder=10)})
        f1.plot(**{**kw, **dict(marker=markerf1, ms=4, zorder=20)})
        plt_me_line.plot(**{**kw, **dict(lw=1, alpha=0.6, color='grey')})

        ############################################################################
        # Axis: Quality factors
        Qs = self.get_quality_factors(swp_variable=swp_variable)
        Qs = Qs if variations is None else Qs[variations]
        Qs = Qs.transpose().sort_index(key=lambda x : x.astype(int))

        ax = axs[1, 0]
        ax.set_title('Quality factors')
        Qs.plot(ax=ax, lw=0, marker=markerf1, ms=4,
                legend=True, zorder=20, color=cmap)
        Qs.plot(ax=ax, lw=1, alpha=0.2, color='grey', legend=False)
        
        df_Qs = np.isinf(Qs)
        # pylint: disable=E1101 
        # Instance of 'ndarray' has no 'values' member (no-member)
        Qs_val = df_Qs.values
        Qs_inf = Qs_val.sum()
        if not (len(Qs) == 0 or Qs_inf > 0): 
          ax.set_yscale('log')

        ############################################################################
        # Axis: Alpha and chi

        axs[0][1].set_title('Anharmonicities (MHz)')
        axs[1][1].set_title('Cross-Kerr frequencies (MHz)')

        def plot_chi_alpha(chi, primary):
            """
            Internal function to plot chi and then also to plot alpha
            """
            idx = pd.IndexSlice
            kw1 = dict(lw=0, ms=4,  marker='o' if primary else 'x')
            kw2 = dict(lw=1, alpha=0.2, color='grey', label='_nolegend_')
            # ------------------------
            # Plot anharmonicity
            ax = axs[0, 1]
            for i, mode in enumerate(mode_idx):  # mode index number, mode index
                alpha = chi.loc[idx[:, mode], mode].unstack(1)
                alpha.columns = [mode]
                alpha.plot(ax=ax, label=mode, color=cmap[i], **kw1)
                if primary:
                    alpha.plot(ax=ax, **kw2)

            # ------------------------
            # Plot chi
            ax = axs[1, 1]
            for mode in mode_idx:  # mode index number, mode index
                # restart the color counter i; n= mode2
                for i, mode2 in enumerate(mode_idx):
                    if int(mode2) > int(mode):
                        chi_element = chi.loc[idx[:, mode], mode2].unstack(1)
                        chi_element.plot(ax=ax, label=f"{mode},{mode2}", color=cmap[i], **kw1)
                        if primary:
                            chi_element.plot(ax=ax, **kw2)

        def do_legends():
            legend_translucent(axs[0][1],  leg_kw=dict(fontsize=7, title='Mode'))
            legend_translucent(axs[1][1],  leg_kw=dict(fontsize=7))

        chiO1 = self.get_chis(variations=variations,
                              swp_variable=swp_variable, numeric=False)
        chiND = self.get_chis(variations=variations,
                              swp_variable=swp_variable, numeric=True)

        use_ND = not np.any(
            [r['fock_trunc'] == None for k, r in self.results.items()])
        if use_ND:
            plot_chi_alpha(chiND, True)
            do_legends()
            plot_chi_alpha(chiO1, False)
        else:
            plot_chi_alpha(chiO1, True)
            do_legends()

        for ax1 in axs:
            for ax in ax1:
                ax.set_xlabel(x_label)

        # Wrap up
        fig.tight_layout()

        return fig, axs

    # Below are functions introduced in v0.8 and newer

    def report_results(self, swp_variable='variation', numeric=True):
        """
        Report in table form the results in a markdown friendly way in Jupyter notebook
        using the pandas interface.
        """

        with pd.option_context('display.precision', 2):
            display(Markdown(("#### Mode frequencies (MHz)")))
            display(Markdown(("###### Numerical diagonalization")))
            display(self.get_frequencies(
                swp_variable=swp_variable, numeric=numeric))

            display(Markdown(("#### Kerr Non-linear coefficient table (MHz)")))
            display(Markdown(("###### Numerical diagonalization")))
            display(self.get_chis(swp_variable=swp_variable, numeric=numeric))

    def get_chis(self, swp_variable='variation', numeric=True, variations: list = None,
                 m=None, n=None):
        """return as multiindex data table

        If you provide m and n as integers or mode labels, then the chi between these modes will
        be returned as a pandas Series.
        """
        label = 'chi_ND' if numeric else 'chi_O1'
        df = pd.concat(self.results.vs_variations(
            label, vs=swp_variable, variations=variations),
            names=[swp_variable])

        if m is None and n is None:
            return df
        else:
            s = df.loc[pd.IndexSlice[:, m], n].unstack(1)[m]
            return s

    def get_frequencies(self, swp_variable='variation', numeric=True, variations: list = None):
        """return as multiindex data table
            index: eigenmode label
            columns: variation label
        """
        label = 'f_ND' if numeric else 'f_1'
        return self.results.vs_variations(label, vs=swp_variable, to_dataframe=True, variations=variations)

    def get_quality_factors(self, swp_variable='variation', variations: list = None):
        """return as pd.Series
            index: eigenmode label
            columns: variation label
        """
        return self.results.vs_variations('Qs', vs=swp_variable, to_dataframe=True, variations=variations)

    def get_participations(self, swp_variable='variation',
                           variations: list = None,
                           inductive=True,
                           _normed=True):
        """

            inductive (bool): EPR for junction inductance when True, else for capacitors

        Returns:
        ----------------
        Returns a multiindex dataframe:
            index 0: sweep variable
            index 1: mode number
            column: junction number

        Example use:
        ---------------
        Plot the participation ratio of all junctions for a given mode vs a sweep of Lj.

        .. code-block language:python

            df=epra.get_participations(swp_variable='Lj')
            df.loc[pd.IndexSlice[:,0],0].unstack(1).plot(marker='o')
        """

        if inductive:
            if _normed:
                getme = 'Pm_normed'
            else:
                getme = 'Pm_raw'
        else:
            if _normed:
                getme = 'Pm_cap'
            else:
                raise NotImplementedError(
                    'not inductive and not _normed not implemented')

        participations = self.results.vs_variations(getme, vs=swp_variable)

        p2 = OrderedDict()
        for key, val in participations.items():
            df = pd.DataFrame(val)
            df.index.name = 'mode'
            df.columns.name = 'junc_idx'
            p2[key] = df

        participations = pd.concat(p2, names=[swp_variable])

        return participations

    def _get_PM_as_DataFrame(self):
        """
        Pm = epra._get_PM_as_DataFrame()
        Pm.unstack(1).groupby(axis=1,level=1).plot()
        """
        Pm = pd.concat(self.PM)
        Pm.index.set_names(['variation', 'mode'], inplace=True)
        Pm.columns.set_names(['junction'], inplace=True)
        return Pm

    def get_ansys_energies(self, swp_var='variation'):
        """
        Return a multi-index dataframe of ansys energies vs swep_variable

        Args:
            swp_var (str) :
        """
        if swp_var == 'variation':
            energies = self.ansys_energies
        else:
            energies = self.get_vs_variable(swp_var, 'ansys_energies')

        df = pd.concat({k: pd.DataFrame(v).transpose()
                        for k, v in energies.items()})
        df.index.set_names([swp_var, 'mode'], inplace=True)

        return df

    def quick_plot_participation(self, mode, junction, swp_variable='variation', ax=None, kw=None):
        """Quick plot participation for one mode

            kw : extra plot arguments
        """
        df = self.get_participations(swp_variable=swp_variable)
        kw = kw or {}
        ax = ax or plt.gca()
        df.loc[pd.IndexSlice[:, mode], junction].unstack(
            1).plot(marker='o', ax=ax, **kw)
        ax.set_ylabel(f'p_({mode},{junction})')

    def quick_plot_frequencies(self, mode, swp_variable='variation', ax=None, kw=None, numeric=False):
        """Quick plot freq for one mode

            kw : extra plot arguments
        """
        kw = kw or {}
        ax = ax or plt.gca()

        s = self.get_frequencies(
            numeric=numeric, swp_variable=swp_variable).transpose()[mode]
        s.plot(marker='o', ax=ax, **kw)

        ax.set_ylabel(f'$\\omega_{mode}$ (MHz)')

    def quick_plot_chi_alpha(self, mode1, mode2, swp_variable='variation', ax=None, kw=None, numeric=False):
        """Quick plot chi between mode 1 and mode 2.

        If you select mode1=mode2, then you will plot the alpha

            kw : extra plot arguments
        """
        kw = kw or {}
        ax = ax or plt.gca()

        s = self.get_chis(swp_variable=swp_variable,
                          numeric=numeric).loc[pd.IndexSlice[:, mode1], mode2].unstack(1)
        s.plot(marker='o', ax=ax, **kw)

        if mode1 == mode2:
            ax.set_ylabel(f'$\\alpha$({mode1}) (MHz) [anharmonicity]')
        else:
            ax.set_ylabel(f'$\\chi$({mode1,mode2}) (MHz) [total split]')

    def quick_plot_mode(self, mode, junction, mode1=None, swp_variable='variation', numeric=False, sharex=True):
        r"""Create a quick report to see mode parameters for only a single mode and a
        cross-kerr coupling to another mode.
        Plots the participation and cross participation
        Plots the frequencie
        plots the anharmonicity

        The values are either for the numeric or the non-numeric results, set by `numeric`
        """

        fig, axs = plt.subplots(2, 2, figsize=(12*0.9, 7*0.9))
        self.quick_plot_frequencies(
            mode, swp_variable=swp_variable, numeric=numeric, ax=axs[0, 1])
        self.quick_plot_participation(
            mode, junction, swp_variable=swp_variable, ax=axs[0, 0])
        self.quick_plot_chi_alpha(mode, mode, numeric=numeric, swp_variable=swp_variable, ax=axs[1, 0],
                                  kw=dict(sharex=sharex))
        if mode1:
            self.quick_plot_chi_alpha(
                mode, mode1, numeric=numeric, swp_variable=swp_variable, ax=axs[1, 1])
            twinax = axs[0, 0].twinx()
            self.quick_plot_participation(mode1, junction, swp_variable=swp_variable, ax=twinax,
                                          kw=dict(alpha=0.7, color='maroon', sharex=sharex))

        for ax in np.ndarray.flatten(axs):
            ax.grid(alpha=0.2)

        axs[0, 1].set_title('Frequency (MHz)')
        axs[0, 0].set_title('Self- and cross-EPR')
        axs[1, 0].set_title('Anharmonicity')
        axs[1, 1].set_title('Cross-Kerr')

        fig.suptitle(f'Mode {mode}', y=1.025)
        fig.tight_layout()

    def quick_plot_convergence(self, ax = None):
        """
        Plot a report of the Ansys convergence vs pass number ona twin axis
        for the number of tets and the max delta frequency of the eignemode.
        """
        ax = ax or plt.gca()
        ax_t = ax.twinx()

        convergence_tets = self.get_convergences_tets_vs_pass()
        convergence_freq = self.get_convergences_max_delta_freq_vs_pass()

        convergence_freq.name = 'Δf'

        plot_convergence_max_df(ax, convergence_freq)
        plot_convergence_solved_elem(ax_t, convergence_tets)


def extract_dic(name=None, file_name=None):
    """#name is the name of the dictionary as saved in the npz file if it is None,
    the function will return a list of all dictionaries in the npz file
    file name is the name of the npz file"""
    with np.load(file_name, allow_pickle=True) as f:
        if name is None:
            return [f[i][()] for i in f.keys()]
        return [f[name][()]]
