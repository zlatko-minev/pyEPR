# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 16:26:21 2020


@author: Jules
"""
import numpy as np
# import pandas as pd
import os
import time
from collections import namedtuple

#from . import logger
import re
from pyEPR.ansys import Optimetrics
from pyEPR import DistributedAnalysis, QuantumAnalysis

# Bayesian optimization 
from bayes_opt import BayesianOptimization, UtilityFunction 
from bayes_opt import SequentialDomainReductionTransformer

### TODO : Creer les classes transmon, memory, purcell, readout, etc 
# class Transmon(): 
#     def __init__(self, freq, Q, Ker, line,...):
#         self.freq = freq

### TODO : compléter les doc

### TODO : ask raph if we should use another timestamp function
### Mettre 0 devant le chiffre si jamais il n'y a qu'un chiffre  
def _timestamp_name(name):
    secondsSinceEpoch = time.time()
    timeObj = time.localtime(secondsSinceEpoch)
    timestamp = '%d%d%d_%d%d%d' % (timeObj.tm_year,timeObj.tm_mon,timeObj.tm_mday,
                                   timeObj.tm_hour, timeObj.tm_min, timeObj.tm_sec)
    return timestamp+'_'+name


class Autotune():
    """ 
    
    pinfo = ProjectInfo("D:/Jules/HFSS_Jules/jules_project",
                        project_name='jules_project', 
                        design_name='autotune', 
                        setup_name='Setup1')
    ### The names used here ('transmon',..) must be exactly the same as thoses in 
    ### target_freqs etc...
    pinfo.junctions['transmon'] = {'rect' : 'ind_JJind_rect',
                                                          'line' : 'ind_JJ_junction_line',
                                                          'Lj_variable' : 'Lj'}
    ### The names used here ('memory',..) must be exactly the same as thoses in 
    ### target_freqs etc...
    pinfo.resonators['memory'] = {'line' : 'mem_line'}
    pinfo.resonators['readout'] = {'line' : 'ro_line'}
    pinfo.resonators['purcell'] = {'line' : 'pur_line'}
    
    pinfo.ports['port50ohm'] = {'rect':'port_50_purcell_track',
                                                 'line': 'port_50_purcell_line',
                                                 'R': 50}
    
    
    ### Must be the same names as in the pinfo.resonators and pinfo.junctions
    target_Qs = {'readout' : 30000, 'purcell' : 300}
    target_freqs = {'transmon' : 5500,  'readout' : 6500, 'purcell' : 6500, 'memory' : 8000}
    target_kers = {'transmon' : 200}
    target_chis = {('transmon', 'memory') : 2, ('transmon', 'readout') : 0.2}
    
    var_name = ['Lj', 'trm_length', 'capa_mem_length', 'capa_ro_length',
                  'shift_capa_mem_trm_Y', 'coupling_dist', 'connector_coupling',
                  'tune_ro']
    
    x = np.array([[7.18123456, 0.4778, 190, 170, 2500, 13, 520, 1.1],
                  [8.18123456, 0.4668, 190, 175, 2460, 13, 540, 1.1]])
    
    a = Autotune(pinfo, target_freqs, target_Qs, 
                     target_kers, target_chis, var_name)
    
    
    The main functions you may use are : 
        sorted_FQK = get_sorted_modes(x)
        cost = cost_function(x)
        x0 = optimize(args (not finished yet))

    
    """
    
    
    # """Class that contains the necessary functions to autotune a chip.
    #     A cost function should be written using thoses classe's instances.
    #     The parameters are X0, target and hfss_info
        
    #     var_name = ['Lj', 'trm_length', 'capa_mem_length', 'capa_ro_length',
    #                 'shift_capa_mem_trm_Y', 'coupling_dist',
    #                 'connector_coupling', 'tune_ro']
        
    #     x = 
        
    #     x0 = {} is the dictionnary of the pysical parameters of the chip, make sure
    #     to give the same names and units as the hfss's ones
        
    #     x0 = {'Lj': '8.5nH', 'trm_length': '0.48mm', 'capa_mem_length' : '190um',
    #           'capa_ro_length' : '170um', 'shift_capa_mem_trm_Y' : '2500um',
    #           'coupling_dist' : '13um', 'connector_coupling' : '520um',
    #           'tune_ro' : '1.1mm'}

        
    #     target = {'Freqs' :{}, 'Qs' : {}, 'Kers' : {}, 'Chis' : {}}
    #     (must contain those 4 keys written like that)
    #     is the dictionnary containing the values of frequencies, QS, Kers
    #     and Chis you are aiming to tune. (Kers and the self-anharmonicity 
    #     frequencies and Chis are the cross terms). The name of the resonators 
    #     and transmons is very important. If you name a resonator 'memory' for 
    #     instance, you must keep using 'memory' for this one in hfss_info.
        
    #     target = {'Freqs' :{'transmon' : 5500, 
    #                         'readout' : 6500,
    #                         'purcell' : 6500,
    #                         'memory' : 8000},
    #               'Qs' : {'readout' : 30000,
    #                       'purcell' : 300},
    #               'Kers' : {'transmon' : 200},
    #               'Chis' : {('transmon', 'memory') : 2,
    #                         ('transmon', 'readout') : 0.2}
    #               }
        
    #     hfss_info = {'hfss_project_path' : "", 'project_name' : '',
    #                  'design_name' : '','setup_name' : 'Setup1',
    #                  'resonators' : {}, 'junctions' : {}, 'ports' : {}}
    #     (must contain those 8 keys written like that)
    #     is the dictionnary providing all the necessary information 
    #     for HFSS to load your project in the right way. The name given to the 
    #     resonators and junctions must be the same here and in target.
        
    #     hfss_info = {'hfss_project_path' : "D:/Jules/HFSS_Jules/jules_project",
    #                  'project_name' : 'jules_project',
    #                  'design_name' : 'autotune',
    #                  'setup_name' : 'Setup1',
    #                  'resonators' : {'memory' : {'line' : 'mem_line'},
    #                                  'readout' : {'line' : 'ro_line'},
    #                                  'purcell' : {'line' : 'pur_line'},
    #                                 },
    #                  'junctions' : {'transmon' : {'rect' : 'ind_JJind_rect',
    #                                                   'line' : 'ind_JJ_junction_line',
    #                                                   'Lj_variable' : 'Lj'}
    #                                 }, 
                                    
    #                  'ports' : {'readout' : {'rect':'port_50_purcell_track',
    #                                          'line': 'port_50_purcell_line',
    #                                          'R': 50}}
    #                  }
               
        
    #     """
    
    
    
    ### TODO il y a toujours un bug qui fait que si on a pas lancé de variation
    ### avant sur HFSS directement et bah epr_hfss n'a pas de fields
    ### TODO : init

    def __init__(self, project_info, target_freqs, target_Qs, target_kers,
                 target_chis, var_name,
                 weight_freqs = None, weight_Qs = None, weight_kers = None,
                 weight_chis = None, cost_fn = None, cost_fn_args = None ):
        '''
        
        Args : 
            ------------------------------------------------------------
            ##### project_info : class from pyEPR.project_info #####
            project_info = ProjectInfo("D:/Jules/HFSS_Jules/jules_project",
                                        project_name='jules_project', 
                                        design_name='autotune', 
                                        setup_name='Setup1')
            
            Then define the junctions, resonators and ports : 
            transmon, memory, purcell and readout in this case
            
            project_info.junctions['transmon'] = {'rect'        : 'ind_JJind_rect',
                                                  'line'        : 'ind_JJ_junction_line',
                                                  'Lj_variable' : 'Lj'}
            
            pinfo.resonators['memory']         = {'line' : 'mem_line'}
            pinfo.resonators['readout']        = {'line' : 'ro_line'}
            pinfo.resonators['purcell']        = {'line' : 'pur_line'}
            
            pinfo.ports['port50ohm']           = {'rect' :'port_50_purcell_track',
                                                          'line': 'port_50_purcell_line',
                                                          'R': 50}

            ProjectInfo should be well initialized before calling Autotune()
            
            
            ##### target_freqs, target_Qs, target_kers, target_chis : dict #####
            target_xxx are 4 dicts containg the values that you want to obtain.
            The keys must be the same strings as those used in pinfo.junctions
            and pinfo.resonators. The keys of target_chis must be tuples of str.
            If you don't want to tune any chis for example, please set target_chis
            to an empty dict.
            target_freqs = {'transmon' : 5500,
                              'readout'  : 6500,
                              'purcell'  : 6500,
                              'memory'   : 8000}

            target_Qs = {'readout' : 30000,
                         'purcell' : 300}
            
            target_kers = {'transmon' : 200}
           
            target_chis = {('transmon', 'memory') : 2,
                           ('transmon', 'readout') : 0.2}
            
            ##### var_name : list of str #####
            It is the list of stringsof HFSS variable names that you wish to tune.
            Make sure that the strings are written the same as in the HFSS design
            
            var_name = ['Lj', 'trm_length', 'capa_mem_length', 'capa_ro_length',
                        'shift_capa_mem_trm_Y', 'coupling_dist', 'connector_coupling',
                        'tune_ro']
            
            ##### weight_freqs, weight_Qs, weight_kers, weight_chis : dict #####
            Optional, see target_xxx for details.
            Weights are used when calculated the cost
            weight_xxx must have the same keys as target_xxx.
            weight_freqs = {'readout' : 1, 'purcell' : 1, 'memory' : 3}
            weight_Qs    = target_Qs = {'readout' : 2, 'purcell' : 2}
            weight_kers  = {}
            weight_chis  = {}
            
            ##### cost_fn and cost_fn_args : function and whatever you want #####
            Optional
            If you want your own cost function, it must take 3 arguments : 
            sorted_FQK (return of sort_modes()), var_list (list of strings such
            as '0', '1', ...) and args which can contain anything.
            If they are None, the method used is the mean_square. 
            


        '''
        ### TODO : isinstance !!!!!!
        assert(isinstance(target_freqs, dict))
        
        self.project_info = project_info
   
    
        Target = namedtuple('Target', ['freqs', 'Qs', 'kers', 'chis'])
        self.target = Target(freqs = target_freqs, 
                             Qs    = target_Qs, 
                             kers  = target_kers, 
                             chis  = target_chis)
        
    
        Weight = namedtuple('Weight', ['freqs', 'Qs', 'kers', 'chis'])
        self.weight = Weight(freqs = weight_freqs, 
                             Qs    = weight_Qs, 
                             kers  = weight_kers, 
                             chis  = weight_chis)


        self.var_name = var_name
        # Get all HFSS variable (usefull to have the units)
        self._all_var = self.project_info.design.get_variables()       
        # Cool way to add the matching units to the HFSS variables
        self._units = [re.match('.*[0-9]([a-zA-Z]*)$',
                         self._all_var[key]).group(1) for key in self.var_name]
        
        
        # Create a directory path to save the parametrics.txt
        design_name = self.project_info.design_name   
        project_path = self.project_info.project_path
        self._dir_path = project_path +'/parametrics_'+ design_name
        os.makedirs(self._dir_path, exist_ok = True)        
   
        
        self.cost_fn = cost_fn if cost_fn else mean_square
        self.cost_fn_args = (self.target,
                             self.weight, cost_fn_args) if cost_fn_args else (self.target,
                             self.weight)
                            
        # self._unsorted_Qs    = pd.DataFrame()
        # self._unsorted_freqs = pd.DataFrame()
        # self._unsorted_chis  = 
        
  
    


    def run_HFSS(self, x) : 
        '''
        This function runs the HFSS variations that correspond to the x
        array. It uses the ansys.Optimetrics package that requires a 
        licence (I think).
        
        Args : 
            x (ndarray) of size (n_variations, n_variables) containing
            the values of the HFSS variables for each variation
               
        Returns : 
            Nothing 
        '''
        if len(x.shape) == 1 :
            x = np.array([x])
        # Number of variations sent to HFSS
        self._nvar = x.shape[0]
        
        # Create a unique partametric_name
        parametric_name = _timestamp_name('parametric')
        print(parametric_name)
        
        # Create a file .txt with the correct shape to import to HFSS
        # This file contains the var_names, the values of x and their units.
        # It is stored in self._dir_path 
        f = open(self._dir_path + '/' + parametric_name +'.txt', 'wt')
        for i in range(len(self.var_name)) : 
            f.write(self.var_name[i] + '; ') 
        f.write('\n')
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                f.write('{}{}; '.format(str(x[i,j]), self._units[j]))
            f.write('\n')
        f.close()

        # Load the Optimetrics module from the ansys package
        opti=Optimetrics(self.project_info.design)
            
        # Import the parametric setup with the list of variations 
        # This function is not pushed on the master branch yet (by manu)
        opti.import_setup(parametric_name, self._dir_path+"\%s.txt"%parametric_name,
                          calc_enable= False)
        
        # Solve the variation, will do 'analyze setup'
        opti.solve_setup(parametric_name)
        print('#####')
        print('HFSS passes are done')
        print('#####')
    
    
    def getFQK(self): 
        ''' 
        This function is post processing the HFSS data and it gets the 
        frequencies, Qs and chi matrix for all variations.
        
        Returns : 
            A namedtuple tuple of 3 panda.DataFrame results with the 
            frequencies, the Qs and the chi matrix. Each DataFrame contains 
            the results of every variation.
            To Use it :
                FQK = a.getFQK
                a.FQK.freqs[var][mode]
                
                where mode is the number of the mode 
                where var is a string such as '0', '1', ..
                (var are the HFSS variation string, you will find them in 
                self.var_list)
                
                you can also get the Qs and the chi matrix
                a.FQK.Qs[var][mode]
                a.FQK.chis[var][mode]
            
            Remark : 
                The chis are returned in an OrderedDict whose values are 
                pandas DataFrames
                        
        '''
        ### Load the results of the HFSS simulation
        self.epr_hfss = DistributedAnalysis(self.project_info)
        nvar = self._nvar
        #Get only the last variations
        self.var_list = self.epr_hfss.variations[-nvar:]
          
        if self.target.chis or self.target.kers: # There is at least 1 junction
            print('There should be at least one junction')
            self.epr_hfss.do_EPR_analysis(self.var_list)
            epr = QuantumAnalysis(self.epr_hfss.data_filename,self.var_list)
            epr.analyze_all_variations(self.var_list)
            
            ###Get the Chis
            unsorted_chis = epr.results.get_chi_O1()
            ### Get the frequencies after the quantum analysis
            unsorted_freqs = epr.get_frequencies()
            ### Get the Qs, if there are no losses, Q = inf
            unsorted_Qs = epr.Qs
           
           
        else : # There is no junction
            print('There must be no junction')
            self.project_info.options.calc_U_E = False
            self.epr_hfss.do_EPR_analysis(self.var_list)
            ### Get the freqs and Qs without the quantum analysis (from HFSS)
            data = self.epr_hfss.get_ansys_frequencies_all().T[self.var_list].T
            
            unsorted_freqs = data['Freq. (GHz)']*1000
            unsorted_freqs.name = 'Freq. (MHz)'
            unsorted_Qs    = data ['Quality Factor']
            unsorted_chis = {}
 
        
        Unsorted = namedtuple('Unsorted', ['freqs', 'Qs', 'chis'])
        self.unsorted_FQK = Unsorted(freqs = unsorted_freqs,
                                     Qs    = unsorted_Qs,
                                     chis  = unsorted_chis)
                    
        print('#####')
        print('The modes frequencies, Qs, Kers, and Chis are calculated but')
        print('unsorted yet. You will find them in the namedtuple : self.unsorted_FQK')
        print('Units : MHz for everything')       
        return self.unsorted_FQK
    ### TODO : verifier que si des trucs sont vides ca sort quand meme bien !!! 
    
    def sort_modes(self, freqs = None, Qs = None, chis = None) : 
        ''' 
        This function is sorting the modes from HFSS. It uses the fact that 
        each resonators should have a line to claculate their participation 
        matrixes.
        Args : 
            freqs panda.DataFrame : should be the result returned by 
            QuantumAnalysis.get_frequencies()
            
            Qs panda.DataFrame : should be the result returned by 
            QuantumAnalysis.Qs
       
            chis collections.OrderedDict : should be the result returned by 
            QuantumAnalysis.get_chi_01()
        
        
        Returns : 
            A namedtuple tuple of 3 panda.DataFrame results with the 
            the sorted results of every variation.
            To Use it :
                sorted_FQK = a.sort_modes
                
                a.sorted_FQK.freqs[var][mode]
                
                where mode is the number of the mode 
                where var is a string such as '0', '1', ..
                (var are the HFSS variation string, you will find them in 
                self.var_list)
                               
                You can also get the sorted Qs, kers and chis
                a.sorted_FQK.Qs[var][mode]
                a.sorted_FQK.kers[var][mode]
                a.sorted_FQK.chis[var][mode]
                        
        '''
        if freqs is None : 
            freqs = self.unsorted_FQK.freqs
        if Qs is None : 
            Qs = self.unsorted_FQK.Qs
        if chis is None : 
            chis = self.unsorted_FQK.chis
        
        
        self._sorted_freqs = {}
        self._sorted_Qs    = {}
        self._sorted_kers  = {}
        self._sorted_chis  = {}
        
        # For each variation, sorting the modes = matching the index of the mode
        # with the name of the mode given in self.target_freqs, sel.target_Qs...
        for var in self.var_list : 
            # Initialization of the dicts
            self._sorted_freqs[var] = {}
            self._sorted_Qs[var]    = {}
            self._sorted_kers[var]  = {}
            self._sorted_chis[var]  = {}
            
            ### Participation matrix of the resonators
            Pm_res = self.epr_hfss.results[var]['Pr']  
        
            ### TODO Ces if/else pourraient être évités si on avait récupéré les matrices
            # de participation dans 'getFQK'
            if self.target.chis or self.target.kers: # There is at least 1 junction
                ### Participation matrix of the transmons
                Pm_trm = self.epr_hfss.results[var]['Pm']
                # Sorting transmon and resonators  modes : finding the highest 
                # participation matrix coefficient to have the index of the mode.
                mode_idxs = Pm_res.idxmax(axis=0).append(Pm_trm.idxmax(axis=0))
            else : 
                mode_idxs = Pm_res.idxmax(axis=0)
        
            # Assign the sorted modes in the sorted dicts
            for mode_name in self.target.freqs : 
                idx = mode_idxs[mode_name]
                self._sorted_freqs[var][mode_name] = freqs[var][idx]
                
            for mode_name in self.target.Qs : 
                idx = mode_idxs[mode_name]
                self._sorted_Qs[var][mode_name] = Qs[var][idx]
            
            for mode_name in self.target.kers : 
                idx = mode_idxs[mode_name]
                self._sorted_kers[var][mode_name] = chis[var][idx][idx]
            
            for mode_name in self.target.chis : 
                idx1 = mode_idxs[mode_name[0]]
                idx2 = mode_idxs[mode_name[1]]
                self._sorted_chis[var][(mode_name[0],
                                mode_name[1])] = chis[var][idx1][idx2]
   
        Sorted = namedtuple('Sorted', ['freqs', 'Qs', 'kers', 'chis'])
        self.sorted_FQK = Sorted(freqs = self._sorted_freqs,
                                 Qs    = self._sorted_Qs,
                                 kers  = self._sorted_kers,
                                 chis  = self._sorted_chis)
        print('#####')
        print('The modes frequencies, Qs, Kers, and Chis are now sorted!')
        print('You will find them in the namedtuple : self.sorted_FQK')
        print('Units : MHz for everything')
        return self.sorted_FQK
    
    def get_sorted_modes(self, x) :
        '''
        This function just runs 3 functions : run_HFSS, getFQK, sort_modes
        Args : 
            see getFQK
        Returns : 
            see sort_modes
        '''
        self.run_HFSS(x)
        self.getFQK()
        return self.sort_modes()        
            
       
    def cost_function(self, x) : 
        """ This is the cost_function that you want to optimize, it should
            return 0 when the FQK from HFSS are the same as the FQK from 
            target. 
            It takes x as an arg 
            if your optimizer can only optimize a certain type of objects, 
            (arrays for instance), you should call it in the following way :

            def realcostfun(x) : 
                x0 = {}...
                return Autotune.cost_function(x0)
        """
        self.get_sorted_modes(x)

        cost = self.cost_fn(self.sorted_FQK, self.var_list, self.cost_fn_args)
        return(cost)
                   
    
    
    def optimize(self, bounds, cost_function = None, niter = 20,
                 method = 'Bayesian'):
        if cost_function == None : 
            cost_function = self.cost_function 

            
        
            
        if method == 'Bayesian' :
            
            
                        
            bounds_transformer = SequentialDomainReductionTransformer()
            optimizer = BayesianOptimization(f = cost_function, pbounds=xybounds, verbose=2, 
                        random_state=1, bounds_transformer = bounds_transformer)
            
            suggested_points = []
            targets = []
            
            for i in range(niter) : 
                optimizer = BayesianOptimization(f = cost_function, pbounds=xybounds, verbose=2, 
                        random_state=1, bounds_transformer = bounds_transformer)
                if suggested_points is not [] :
                    optimizer.register(suggested_points[-1], targets[-1])
                utility = UtilityFunction(kind="ucb", xi=0.0)
                suggested_points.append(optimizer.suggest(utility))
                targets.append(cost_function(suggested_points[-1]))
                
            
            
            
            
            
            bounds_transformer = SequentialDomainReductionTransformer()
            opti = BayesianOptimization(f = cost_function, pbounds=bounds, verbose=2, 
                        random_state=1, bounds_transformer = bounds_transformer)
    
            opti.maximize(int(niter/4), 3*int(niter/4)) 
            return(opti.max)
        
         
        
def mean_square(sorted_FQK, var_list, args):
    '''
    This functioun calculates the mean square betwenn the targets and the 
    result of sort_modes. It can be weighted.
    Args : 
        ##### sorted_FQK : namedtuple #####
        See sort_modes for description
        
        ##### var_list #####
        list of strings
        
        ##### args #####
        tuple of argumenents. This function only supports 2-tuples containg 
        target an weight (2 namedtuples)
        
    Returns : 
        The cost (mean square) should be 0 when the x parameters are tuned so 
        that the target matches the sorted_FQK/
      

    '''
    target = args[0]
    weight = args[1]
    cost = 0
    ### TODO C'est degeu ce if else 
    if weight.freqs == weight.Qs == weight.kers == weight.chis == None : 
        for var in var_list :
            for i in range(4) : 
                for key in target[i].keys() : 
                    cost += ((target[i][key] - sorted_FQK[i][var][key])/target[i][key])**2
    
    else : 
        for var in var_list :
            for i in range(4) : 
                for key in target[i].keys() : 
                    cost += weight[i][key]*((target[i][key] - sorted_FQK[i][var][key])/target[i][key])**2

                
    return cost
    
        
        
        
        
        
        
        

        
        
        