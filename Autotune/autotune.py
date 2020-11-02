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
    Class made to autotune a chip
    See autotune_test_optimizatoin for an example of use
    Steps : 
        1) Create a HFSS design. Make sure to add a line on each resonator you
           wish to tune. Make sure to run the HFSS setup once.
           
        2) Inquire the project_info but add the resonators too: 
            pinfo.resonators['resonator1'] = {'line' : 'line_res1'}
        
        If you just want to get HFSS modes, follow steps a.
        If you want to autotune your chip, follow steps b.
        
        3.a) Inquire target_Qs, target_freqs... to tell which hfss modes 
             you wish to get. (see init for syntax)
        
        4.a) Inquire var_name and x that will contain the names of the HFSS 
             variables that you wish to set to the value x.
             
        5.a) Launch auto = Autotune(pinfo, target_freqs, target_Qs, 
                 target_kers, target_chis, var_name) and do :
                                    
             auto.get_sorted_modes(x)
             
        3.b) Inquire target_Qs, target_freqs... to tell which hfss modes 
             you wish to optimize. (see init for syntax)
        
        4.b) Inquire var_name that will contain the names of the variables that
             will be tuned. Try to limit yourself to 4 variables max because
             optimizing more than 4 parameters if very long and unprecise. 
             You should probably optimize your parameters 3 by 3.
             
        5.b) Launch auto = Autotune(pinfo, target_freqs, target_Qs, 
                 target_kers, target_chis, var_name)
        
        6.b) Inquire bounds that will set the bounds for the variables and 
             make sure that bounds.keys are the same as var_name.
        
        7.b) Do auto.optimize(bounds)
                                    
    """
    
    
    
    
    ### TODO il y a toujours un bug qui fait que si on a pas lancé de variation
    ### avant sur HFSS directement et bah epr_hfss n'a pas de fields

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
            Pm_res = self.epr_hfss.resiults[var]['Pr']  
        
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
        """ This is the cost_function that you want to optimize
            Args : 
            x : see run_HFSS
            
            Returns : 
                the cost that should be 0 when your chip is optimimized.
                (when the FQK from HFSS are the same as the FQK from 
                target)
            
         
        """
        self.get_sorted_modes(x)

        cost = self.cost_fn(self.sorted_FQK, self.var_list, self.cost_fn_args)
        return(cost)
                   
    
    
    
        
    
    def optimize(self, bounds, cost_function = None, niter = 20,
                 method = 'Bayesian'):
        """
        This function should optimize your chip

        Parameters
        ----------
        bounds : dict
            The keys must be the same strings as the HFSS variables contained
            in var_name. The values are tuples (x,y) where x is the min bound
            associated to the key and y is the max bound.
            
        cost_function : function, optional
            The default is a mean square method. You should better change 
            cost_fn instead of cost_function (see __init__)
            
        niter : int, optional
            Number of passes. The default is 20. I should probably find another
            criteria such as the tolerance.
            
        method : string, optional
            Optimization method, for the moment, there is only the
            Bayesian Optimization.

        Returns
        -------
        A dictionnary containing as keys : 
            target : the best cost
            params : a dict containing the parameters that gave this cost
        """
    
    
        if method == 'Bayesian' :
            def bayesian_cost(**kwargs):
                x = np.array([kwargs[name] for name in self.var_name ])
                return -self.cost_function(x)
     
    
            bounds_transformer = SequentialDomainReductionTransformer()
            self.optimizer = BayesianOptimization(f = bayesian_cost, pbounds=bounds)
            
            
            self.optimizer.maximize(int(niter/3), int(2*niter/3))
            
            self.x0 = self.optimizer.max
            
            return self.x0
            
            # Work in progress for parrallelizing the optimization
         
            # guess = {}
            # for key, val in bounds.items() : 
            #     guess[key] = (val[1]-val[0])/3
            
            # self.suggested_points = []
            # self.targets = []
            # next_point = []
            
            
            # for i in range(niter) : 
            #     utility = UtilityFunction(kind="ucb", kappa = 2.576, xi=0.0)
                
            #     next_point = optimizer.suggest(utility)
            #     self.suggested_points.append(next_point)
            #     self.targets.append(bayesian_cost(**next_point))
                
            #     optimizer.register(next_point, self.targets[-1])
            #     if i == 0 : 
            #         optimizer.register(guess, bayesian_cost(**guess))
            
            #     next_point = []

            # return(targets, suggested_points)      
            
         
        
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
    costs = []
    ### TODO C'est degeu ce if else 
    if weight.freqs == weight.Qs == weight.kers == weight.chis == None : 
        for var in var_list :
            cost = 0
            for i in range(4) : 
                for key in target[i].keys() : 
                    cost += ((target[i][key] - sorted_FQK[i][var][key])/target[i][key])**2
            costs.append(cost)
    else : 
        for var in var_list :
            cost = 0
            for i in range(4) : 
                for key in target[i].keys() : 
                    cost += weight[i][key]*((target[i][key] - sorted_FQK[i][var][key])/target[i][key])**2
            costs.append(cost)
                
    if len(costs) == 1 :
        return costs[0]
    else : 
        return costs
    
        
        
        
        
        
        
        

        
        
        